# Development notes

This document stores implementation notes that are not needed in the main README.

## Validated workflow components

The tool was tested on FASTA-derived human precursor tables and Koina-based Prosit/MS2PIP predictions.

Validated functionality includes:

- FASTA digest and precursor-space generation;
- Skyline-like precursor-space auditing;
- optional selenocysteine handling during digest;
- Koina long-format prediction output parsing;
- chunked prediction with resume;
- minimal DIA-library output by default;
- DeepLC RT integration;
- one-command `compactlib build` workflow;
- Prosit max7, MS2PIP max7 and Union14 generation.

## Release philosophy

The main README should describe the product-facing workflow. Development history, debugging details and benchmark-specific notes should stay in `docs/` or supplementary materials.
