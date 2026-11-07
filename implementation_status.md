## DONE:
1. Parsing the FASTA files
---
2. Encoding strings into base-4 system.
---
3. Hashing and building the k-mers from an encoded sequence.
---
4. Building the reference index.
---
6. Finding hits in the reference of a given sequence
---

## TODO:
5. Reverse complement so reads from both strands can be mapped.
---
7. Diagonal voting: compute `ref_pos - read_pos` for each hit, find the densest `BAND`-wide window and the best non-overlapping runner-up.
---
8. Semi-global edit-distance DP fill the DP row by row with NumPy, get the edit distance and end from the last row, the start from the traceback.
---
9. Run mapping, accept only if edit distance ≤ `MAX_ERR_RATE` × read length and the best candidate beats the runner-up by `MIN_ED_GAP`.
---
10. Run `main.py` end to end on `files/example_files1` and check the output format and the `.log` file.
---
11. Evaluation script (`evaluate.py`): compare the output with the expected results (% mapped, % wrongly mapped within 20 bp, time per read).
---
12. Tune the parameters in `config.py` (`K`, `MAX_OCC`, `BAND`, `MIN_VOTES`, `MARGIN`, `MAX_ERR_RATE`, `MIN_ED_GAP`) against the targets: ≥ 80% mapped, ≤ 1% wrong, ≤ 100 + 3r seconds, < 1 GB memory.
---
13. Update `README.md` to match the code
---
14. Fill in `AI_USAGE.md` (scope of AI help and the list of prompts).
---
15. Prepare the slides describing the approach.

