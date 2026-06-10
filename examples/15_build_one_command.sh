#!/usr/bin/env bash
set -euo pipefail

# One-command black-box smoke test with the toy FASTA and mock backends.
# Replace --prediction-backend mock with koina and --rt-backend mock with deeplc
# for real Prosit/MS2PIP/DeepLC prediction.

OUT=toy_out/build_one_command
rm -rf "$OUT"
mkdir -p "$OUT"

compactlib build \
  --fasta examples/toy.fasta \
  --output-dir "$OUT" \
  --prefix toy \
  --sample-size 20 \
  --preset generic-dia \
  --top-n 3 \
  --chunk-size 5 \
  --batch-size 5 \
  --prediction-backend mock \
  --rt-backend mock \
  --rt-output-column PredictedRetentionTime \
  --resume \
  --verbose

printf '\nCreated files:\n'
ls -lh "$OUT"/*.tsv
