#!/usr/bin/env bash
set -euo pipefail

mkdir -p toy_out

compactlib digest \
  --fasta examples/toy.fasta \
  --output toy_out/toy_precursors_mox.tsv \
  --enzyme trypsin-p \
  --missed-cleavages 1 \
  --min-length 5 \
  --max-length 30 \
  --charges 2 \
  --variable-oxidation-m \
  --max-variable-mods 1

compactlib predict \
  --precursors toy_out/toy_precursors_mox.tsv \
  --output toy_out/toy_mock_max3.tsv \
  --backend mock \
  --fragment-types b,y \
  --fragment-charges 1 \
  --top-n 3

python - <<'PY'
import pandas as pd
for p in ["toy_out/toy_precursors_mox.tsv", "toy_out/toy_mock_max3.tsv"]:
    print("\n", "="*80, "\n", p, sep="")
    df = pd.read_csv(p, sep="\t")
    print(df.head(20).to_string(index=False))
PY
