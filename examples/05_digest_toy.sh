#!/usr/bin/env bash
set -euo pipefail

mkdir -p toy_out

compactlib digest \
  --fasta examples/toy.fasta \
  --output toy_out/toy_precursors.tsv \
  --enzyme trypsin-p \
  --missed-cleavages 1 \
  --min-length 5 \
  --max-length 30 \
  --charges 2,3 \
  --c-mod-format unimod

python - <<'PY'
import pandas as pd
print(pd.read_csv('toy_out/toy_precursors.tsv', sep='\t').to_string(index=False))
print('\nSummary:')
print(pd.read_csv('toy_out/toy_precursors.summary.tsv', sep='\t').to_string(index=False))
PY
