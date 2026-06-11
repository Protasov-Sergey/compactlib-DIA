#!/usr/bin/env bash
set -euo pipefail

# Example: FASTA digest -> Koina Prosit prediction -> compact max7 library.
# Requires: pip install -e '.[koina]'
# For local Koina server, use e.g. --server-url localhost:8500 --no-ssl

mkdir -p toy_out

compactlib digest \
  --fasta examples/toy.fasta \
  --output toy_out/toy_precursors_for_koina.tsv \
  --enzyme trypsin-p \
  --missed-cleavages 1 \
  --min-length 5 \
  --max-length 30 \
  --charges 2 \
  --variable-oxidation-m \
  --max-variable-mods 1

compactlib predict \
  --precursors toy_out/toy_precursors_for_koina.tsv \
  --output toy_out/toy_prosit_koina_max7.tsv \
  --backend koina \
  --model Prosit_2020_intensity_HCD \
  --server-url koina.wilhelmlab.org:443 \
  --ssl \
  --collision-energy 30 \
  --batch-size 256 \
  --sequence-source stripped \
  --top-n 7
