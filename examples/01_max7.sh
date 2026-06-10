#!/usr/bin/env bash
set -euo pipefail

compactlib max-n \
  --input Full_prosit_library.csv.gz \
  --output Prosit_max7.tsv \
  --top-n 7
