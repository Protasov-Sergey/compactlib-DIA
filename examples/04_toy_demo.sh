#!/usr/bin/env bash
set -euo pipefail

# Run from repository root after: pip install -e '.[dev]'
# This toy example uses top-n=3 to make the algorithm easy to inspect by eye.

mkdir -p toy_out

compactlib max-n \
  --input examples/toy_prosit_full.tsv \
  --output toy_out/toy_prosit_max3.tsv \
  --top-n 3

compactlib union \
  --input-a examples/toy_prosit_full.tsv \
  --input-b examples/toy_ms2pip_full.tsv \
  --output toy_out/toy_union6.tsv \
  --top-n 3

compactlib consensus \
  --input-a examples/toy_prosit_full.tsv \
  --input-b examples/toy_ms2pip_full.tsv \
  --output toy_out/toy_consensus3.tsv \
  --top-n 3

compactlib random-n \
  --input examples/toy_prosit_full.tsv \
  --output toy_out/toy_prosit_random3.tsv \
  --top-n 3 \
  --seed 42

compactlib reverse-intensity \
  --input toy_out/toy_prosit_max3.tsv \
  --output toy_out/toy_prosit_max3_reverse.tsv

python examples/inspect_toy_outputs.py
