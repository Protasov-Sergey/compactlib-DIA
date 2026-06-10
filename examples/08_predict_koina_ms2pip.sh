#!/usr/bin/env bash
set -euo pipefail

# Example: predict compact max7 library using Koina MS2PIP HCD model.
# Requires: pip install -e '.[koina]'
# Note: model input/output names can be adjusted if a Koina model exposes a different schema.

mkdir -p toy_out

compactlib digest \
  --fasta examples/toy.fasta \
  --output toy_out/toy_precursors_for_ms2pip.tsv \
  --enzyme trypsin-p \
  --missed-cleavages 1 \
  --min-length 5 \
  --max-length 30 \
  --charges 2

compactlib predict \
  --precursors toy_out/toy_precursors_for_ms2pip.tsv \
  --output toy_out/toy_ms2pip_koina_max7.tsv \
  --backend koina \
  --model ms2pip_HCD2021 \
  --server-url koina.wilhelmlab.org:443 \
  --ssl \
  --collision-energy 30 \
  --batch-size 256 \
  --sequence-source stripped \
  --top-n 7
