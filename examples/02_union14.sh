#!/usr/bin/env bash
set -euo pipefail

compactlib union \
  --input-a Full_prosit_library.csv.gz \
  --input-b Full_ms2pip_library.csv.gz \
  --output Union14.tsv \
  --top-n 7
