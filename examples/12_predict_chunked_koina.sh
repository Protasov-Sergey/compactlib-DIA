#!/usr/bin/env bash
set -euo pipefail

PREC=${PREC:-human_precursors.tsv}
OUT=${OUT:-Prosit_max7.fullhuman.tsv}
WORK=${WORK:-Prosit_max7.chunked_work}

compactlib predict-chunked \
  --precursors "$PREC" \
  --output "$OUT" \
  --work-dir "$WORK" \
  --chunk-size 100000 \
  --backend koina \
  --model Prosit_2020_intensity_HCD \
  --server-url koina.wilhelmlab.org:443 \
  --ssl \
  --collision-energy 30 \
  --batch-size 512 \
  --sequence-source stripped \
  --preset generic-dia \
  --top-n 7 \
  --resume \
  --verbose
