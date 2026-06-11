# compactlib toy example

This toy dataset demonstrates the behavior of `compactlib` on three precursors.

Files:

- `toy_prosit_full.tsv` — toy Prosit-like full library.
- `toy_ms2pip_full.tsv` — toy MS2PIP-like full library.
- `04_toy_demo.sh` — runs `max-n`, `union`, `consensus`, `random-n`, and `reverse-intensity`.
- `inspect_toy_outputs.py` — prints compact outputs in a readable form.

Use top-n=3 in this example, because it is easier to inspect by eye.

Expected behavior:

1. `max-n`
   Keeps up to 3 most intense transitions per precursor.
   For `ACDEK/2`, the Prosit-like library has only 2 transitions, so both are retained.

2. `union`
   Computes top3 from library A and top3 from library B, merges them, and removes duplicate transitions
   by annotation-based transition key:
   ModifiedPeptide + PrecursorCharge + FragmentType + FragmentSeriesNumber + FragmentCharge.
   For top-n=3 this is Union6.

3. `consensus`
   Scores each transition as:
   score = 1/rank_a + 1/rank_b.
   If a transition is absent in one model, the missing contribution is 0.
   Rows are taken from library A when the transition exists in both libraries.

4. `random-n`
   Randomly selects up to 3 transitions per precursor with a fixed seed.

5. `reverse-intensity`
   Preserves the same transitions but reverses intensity values within each precursor.
