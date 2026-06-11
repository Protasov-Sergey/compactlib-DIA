#!/usr/bin/env bash
set -euo pipefail

mkdir -p toy_out

cat > toy_out/toy_selenoprotein.fasta <<'EOF'
>sp|PSEL|PSEL_HUMAN toy selenoprotein
MPEPTIDEKAEENITESCQUR
EOF

compactlib digest \
  --fasta toy_out/toy_selenoprotein.fasta \
  --output toy_out/toy_precursors_with_u.tsv \
  --preset generic-dia \
  --allow-selenocysteine \
  --c-mod-format unimod

cat toy_out/toy_precursors_with_u.summary.tsv
