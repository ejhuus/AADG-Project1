# Read Mapper — Assignment 1 (Algorithms for Genomic Data Analysis)

A seed-and-extend read mapper for long, noisy reads (~1 kbp, 5–10% errors) against a reference of up to 20 Mbp.

---

## Usage

```bash
python3 mapper.py reference.fasta reads.fasta output.txt
```

| Argument | Description |
|---|---|
| `reference.fasta` | Reference sequence in FASTA format (≤ 20 Mbp) |
| `reads.fasta` | Reads in FASTA format (~1 kbp each) |
| `output.txt` | Output file, one line per mapped read |

### Output format

One tab-separated line per **mapped** read (unmapped reads are omitted):

```
<read_id>    <start>    <end>
```

Progress and timing information is printed and saved in `.log` file.

---

## Pipeline

The pipeline schema for this mapper.


Index all reference k-mers in a sorted table, find exact k-mer seeds for each read, cluster the seeds by diagonal to locate the read, then verify and refine the location with semi-global edit-distance dynamic programming.

```
reference.fasta ─► [1] Read FASTA ─► [2] Encode ─► [3] Build k-mer index   (once)
                                                            │
reads.fasta ─► [1] Read FASTA ─► for each read:             ▼
                 [4] Encode read + reverse complement
                 [5] Seeding        → exact k-mer hits (ref_pos, read_pos)
                 [6] Diagonal voting → candidate locations
                 [7] DP verification → edit distance, start, end
                 [8] Accept / reject
                 [9] Write output line
```

### 1. Input — FASTA parsing
Parse the reference into a single sequence and stream the reads as `(id, sequence)` pairs.

| Function | Purpose |
|---|---|
| `read_reference(path)` | Read the reference sequence (concatenate lines) |
| `read_fasta(path)` | Generator yielding `(read_id, sequence)` |

### 2. Encoding
Convert DNA strings into NumPy `uint8` arrays (`A=0, C=1, G=2, T=3`, other characters such as `N` = 4) using a 256-entry lookup table, so all later steps can be vectorized.

| Function | Purpose |
|---|---|
| `encode(seq)` | String → `uint8` array |
| `revcomp(enc)` | Reverse complement in encoded space (`3 - x`, reversed; `N` stays `N`) |

### 3. Index construction (once)
Compute every k-mer of the reference as a 2-bit-packed `uint32` integer (k ≤ 16), mask out k-mers containing `N`, and sort them with `np.argsort`. All occurrences of the same k-mer are then contiguous, so any k-mer can be found by binary search, just as in a suffix array.

| Function | Purpose |
|---|---|
| `kmer_codes(enc, k)` | Vectorized rolling 2-bit k-mer codes + validity mask |
| `build_index(ref_enc, k)` | Returns `sorted_codes` and `positions` arrays |

Memory: ~160 MB for a 20 Mbp reference.

### 4. Read encoding
Encode the read and its reverse complement; both strands are mapped and the better result is kept.

*(Uses `encode` and `revcomp` from step 2.)*

### 5. Seeding
Compute all k-mers of the read and look them all up at once with `np.searchsorted`. Over-represented k-mers (more than `MAX_OCC` occurrences, typically repeats) are ignored. Each match becomes a hit `(ref_pos, read_pos)`.

| Function | Purpose |
|---|---|
| `find_hits(read_enc, index)` | Returns arrays `ref_pos`, `read_pos` of all exact seed hits |

With 10% errors and k = 15, a true location gets on the order of 100–200 hits, while random hits are scattered (~0.02 per k-mer).

### 6. Candidate selection — diagonal voting
Each hit votes for the diagonal `ref_pos - read_pos` (the estimated read start). Hits from the true location share nearly the same diagonal, shifted only slightly by indels. A sliding window of width `BAND` over the sorted diagonals finds the densest cluster and the best non-overlapping runner-up.

| Function | Purpose |
|---|---|
| `best_diagonals(ref_pos, read_pos)` | Returns the top candidate windows `(est_start, votes)` |

Reads whose best window has fewer than `MIN_VOTES` hits are reported as unmapped immediately.

### 7. Verification — semi-global edit-distance DP
For each surviving candidate, extract the reference window `[est_start - MARGIN, est_start + len(read) + MARGIN]` and run edit-distance DP with a zero first row (the read may start anywhere in the window). Each row is computed with NumPy: the up/diagonal terms come from the previous row, and the left-neighbour term is resolved with a running minimum:

```
D[i][j] - j = min(t[j] - j, D[i][j-1] - (j-1))   →   np.minimum.accumulate
```

The minimum of the last row gives the edit distance and the end position; a traceback gives the start position.

| Function | Purpose |
|---|---|
| `semiglobal_dp(read, window)` | Fills the DP matrix; returns `(D, edit_distance, end_col)` |
| `traceback(D, read, window, end_col)` | Returns the start column of the alignment |
| `verify(candidate, ref_enc)` | Builds the window, runs DP, returns `(edit_distance, start, end)` |

### 8. Decision
Accept the best candidate only if:
- its edit distance is at most `MAX_ERR_RATE × read_length`, and
- it is clearly better than the runner-up, by at least `MIN_ED_GAP` (avoids mapping repeats to the wrong copy).

Otherwise the read is left unmapped. Leaving a read unmapped is preferred over risking an incorrect mapping.

| Function | Purpose |
|---|---|
| `map_read(read_enc, index, ref_enc)` | Runs steps 4–8 for one read; returns `(start, end)` or `None` |

### 9. Output
| Function | Purpose |
|---|---|
| `main()` | Parses arguments, builds the index, loops over the reads, writes results, logs timings |

---

## Parameters

Defined as constants at the top of `mapper.py`.

| Parameter | Typical value | Meaning |
|---|---|---|
| `K` | 13–15 | k-mer (seed) length |
| `MAX_OCC` | 20–50 | Ignore k-mers occurring more often than this |
| `BAND` | 60–100 bp | Width of the diagonal voting window |
| `MIN_VOTES` | 8–15 | Minimum seed hits to consider a location |
| `MARGIN` | 100–150 bp | Extra reference on each side of the DP window |
| `MAX_ERR_RATE` | 0.15–0.17 | Maximum accepted edit distance / read length |
| `MIN_ED_GAP` | tuned | Required edit-distance gap between best and second-best candidate |

---

## Performance targets

| Stage | Expected cost |
|---|---|
| Index build (20 Mbp) | ~5–15 s |
| Seeding + voting per read | a few ms |
| DP verification per read | ~20–100 ms |
| Total per read | well under 0.2 s |
| Peak memory | < 500 MB |

Assignment limits: total time ≤ 100 + 3r seconds, memory < 1 GB, ≥ 80% of mappable reads mapped, ≤ 1% incorrectly mapped. A read counts as correctly mapped if its coordinates are within 20 bp of the true ones.

---

## Repository structure

```
mapper.py      # the mapper (all steps above, organized in sections)
evaluate.py    # helper: compares output with expected results (% mapped, % wrong, time per read)
slides/        # presentation describing the approach
AI_USAGE.md    # scope of AI usage and list of prompts
README.md
```

## AI usage

The scope of AI assistance and the list of prompts are documented in `AI_USAGE.md`, as required by the assignment.
