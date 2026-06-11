#!/usr/bin/env bash
set -euo pipefail

# Practical precursor candidate filtering example.
# This is useful for predicted DIA library generation, where exhaustive charge
# enumeration can create unrealistic candidates (e.g. very short peptides with
# high precursor charges).

mkdir -p toy_out

compactlib digest \
  --fasta examples/toy.fasta \
  --output toy_out/toy_precursors_filtered.tsv \
  --enzyme trypsin-p \
  --missed-cleavages 1 \
  --min-length 7 \
  --max-length 30 \
  --charges 2,3,4 \
  --min-precursor-mz 400 \
  --max-precursor-mz 1200 \
  --charge-length-rules 4:12 \
  --c-mod-format unimod

column -t -s $'\t' toy_out/toy_precursors_filtered.summary.tsv
