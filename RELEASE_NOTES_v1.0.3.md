# compactlib-DIA v1.0.3

Publication-ready release of compactlib-DIA.

This release provides a one-command workflow for compact in silico DIA spectral library generation:

FASTA -> digest -> RT prediction -> Prosit maxN -> MS2PIP maxN -> Union2N

Main entry point:

compactlib build

Included functionality:

- FASTA-derived precursor generation
- optional DeepLC-based RT prediction
- Koina-based Prosit and MS2PIP intensity prediction
- chunked prediction with resume support
- compact top-N transition selection
- Union2N integrated library construction
- minimal DIA-library output format
- reproducible summary and parameter files

Publication-readiness updates in v1.0.3:

- explicit MIT LICENSE file
- GitHub Actions CI workflow
- polished README
- updated package metadata

The core compactlib build workflow is unchanged relative to v1.0.2.
