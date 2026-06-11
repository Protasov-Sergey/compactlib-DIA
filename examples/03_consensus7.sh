#!/usr/bin/env bash
set -euo pipefail

compactlib consensus \
  --input-a Full_prosit_library.csv.gz \
  --input-b Full_ms2pip_library.csv.gz \
  --output Consensus_max7.tsv \
  --top-n 7
