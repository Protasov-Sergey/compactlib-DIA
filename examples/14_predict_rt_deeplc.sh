#!/usr/bin/env bash
set -euo pipefail

# Example: add DeepLC RT predictions to a compactlib precursor table.
# Requires: pip install -e '.[deeplc]'

compactlib digest \
  --fasta examples/toy.fasta \
  --output toy_out/toy_precursors.tsv \
  --enzyme trypsin-p \
  --missed-cleavages 1 \
  --min-length 5 \
  --max-length 30 \
  --charges 2 \
  --quiet

compactlib predict-rt \
  --precursors toy_out/toy_precursors.tsv \
  --output toy_out/toy_precursors.with_rt.tsv \
  --backend deeplc \
  --rt-output-column NormalizedRetentionTime \
  --write-deeplc-input toy_out/toy_deeplc_input.tsv

compactlib predict \
  --precursors toy_out/toy_precursors.with_rt.tsv \
  --output toy_out/toy_mock_max3.with_rt.tsv \
  --backend mock \
  --fragment-types b,y \
  --fragment-charges 1 \
  --top-n 3 \
  --include-rt \
  --rt-column NormalizedRetentionTime
