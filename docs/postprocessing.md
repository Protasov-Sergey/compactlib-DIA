# Post-processing existing libraries

In addition to the FASTA-to-library workflow, `compactlib-DIA` can post-process existing spectral libraries.

## Select top-N transitions

```bash
compactlib max-n \
  --input full_library.tsv \
  --output library_max7.tsv \
  --top-n 7
```

## Build Union2N library

```bash
compactlib union \
  --input-a Prosit_max7.tsv \
  --input-b MS2PIP_max7.tsv \
  --output Union14.tsv \
  --top-n 7
```

## Build consensus library

```bash
compactlib consensus \
  --input-a Prosit_max7.tsv \
  --input-b MS2PIP_max7.tsv \
  --output Consensus_max7.tsv \
  --top-n 7
```

## Control libraries

Random transition control:

```bash
compactlib random-n \
  --input library.tsv \
  --output library_random7.tsv \
  --top-n 7 \
  --seed 42
```

Reverse-intensity control:

```bash
compactlib reverse-intensity \
  --input library_max7.tsv \
  --output library_max7_reverse.tsv
```
