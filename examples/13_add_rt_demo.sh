#!/usr/bin/env bash
set -euo pipefail

# Demonstrate RT passthrough and post hoc RT merging using the mock backend.
# This example uses tiny toy files; real RT tables can come from external tools
# such as DeepLC.

OUTDIR="${1:-toy_out}"
mkdir -p "$OUTDIR"

compactlib digest \
  --fasta examples/toy.fasta \
  --output "$OUTDIR/toy_precursors.tsv" \
  --enzyme trypsin-p \
  --missed-cleavages 1 \
  --min-length 5 \
  --max-length 30 \
  --charges 2 \
  --quiet --no-progress

OUTDIR_FOR_PY="$OUTDIR" python - <<'PY'
import os
import pandas as pd
from pathlib import Path
outdir = Path(os.environ["OUTDIR_FOR_PY"])
prec = outdir / "toy_precursors.tsv"
df = pd.read_csv(prec, sep="\t")
rt = df[["ModifiedPeptide", "PrecursorCharge"]].drop_duplicates().copy()
rt["Tr_recalibrated"] = range(1, len(rt) + 1)
rt.to_csv(outdir / "toy_rt.tsv", sep="\t", index=False)
# Also create a precursor table containing RT for passthrough mode.
df = df.merge(rt, on=["ModifiedPeptide", "PrecursorCharge"], how="left")
df.to_csv(outdir / "toy_precursors_with_rt.tsv", sep="\t", index=False)
PY

compactlib predict \
  --precursors "$OUTDIR/toy_precursors_with_rt.tsv" \
  --output "$OUTDIR/toy_mock_max3_with_rt_passthrough.tsv" \
  --backend mock \
  --fragment-types b,y \
  --fragment-charges 1 \
  --top-n 3 \
  --include-rt \
  --rt-column Tr_recalibrated \
  --rt-output-column NormalizedRetentionTime \
  --quiet --no-progress

compactlib predict \
  --precursors "$OUTDIR/toy_precursors.tsv" \
  --output "$OUTDIR/toy_mock_max3.tsv" \
  --backend mock \
  --fragment-types b,y \
  --fragment-charges 1 \
  --top-n 3 \
  --quiet --no-progress

compactlib add-rt \
  --library "$OUTDIR/toy_mock_max3.tsv" \
  --rt-table "$OUTDIR/toy_rt.tsv" \
  --output "$OUTDIR/toy_mock_max3_with_rt_merged.tsv" \
  --key ModifiedPeptide,PrecursorCharge \
  --rt-column Tr_recalibrated \
  --rt-output-column NormalizedRetentionTime \
  --min-match-rate 1.0 \
  --quiet --no-progress

printf '\nCreated RT example outputs in %s\n' "$OUTDIR"
