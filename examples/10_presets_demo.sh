#!/usr/bin/env bash
set -euo pipefail

# Demonstrate Skyline/DIA-NN-like practical filtering presets.
# Run from repository root after: pip install -e '.[dev]'

mkdir -p toy_out

echo "Available presets:"
compactlib presets

# Project preset: trypsin-p, 1 missed cleavage, length 7-30,
# charges 2,3,4, precursor m/z 400-1000.
compactlib digest \
  --fasta examples/toy.fasta \
  --output toy_out/toy_precursors_project_preset.tsv \
  --preset generic-dia

# The user can override any preset parameter transparently.
compactlib digest \
  --fasta examples/toy.fasta \
  --output toy_out/toy_precursors_project_preset_1200.tsv \
  --preset generic-dia \
  --max-precursor-mz 1200

# Prediction preset applies transition filters after prediction.
# Here mock backend is used only for demonstration/testing.
compactlib predict \
  --precursors toy_out/toy_precursors_project_preset.tsv \
  --output toy_out/toy_mock_project_preset_max7.tsv \
  --backend mock \
  --preset generic-dia \
  --top-n 7

printf "\nDigest summary:\n"
cat toy_out/toy_precursors_project_preset.summary.tsv

printf "\nPredict summary:\n"
cat toy_out/toy_mock_project_preset_max7.summary.tsv
