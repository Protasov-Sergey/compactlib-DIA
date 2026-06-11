# Full workflow

This document describes the complete FASTA-to-library workflow used by `compactlib-DIA`.

## Recommended one-command workflow

```bash
compactlib build \
  --fasta human.fasta \
  --output-dir compactlib_out \
  --prefix human \
  --preset generic-dia \
  --top-n 7 \
  --prediction-backend koina \
  --prosit-model Prosit_2020_intensity_HCD \
  --ms2pip-model ms2pip_HCD2021 \
  --server-url koina.wilhelmlab.org:443 \
  --ssl \
  --collision-energy 30 \
  --sequence-source stripped \
  --rt-backend deeplc \
  --rt-output-column PredictedRetentionTime \
  --chunk-size 100000 \
  --batch-size 512 \
  --resume \
  --verbose
```

## Main stages

```text
FASTA
  -> digest
  -> precursor table
  -> optional DeepLC RT prediction
  -> Prosit maxN prediction
  -> MS2PIP maxN prediction
  -> Union2N integrated library
```

## Output files

```text
<prefix>.precursors.all.tsv
<prefix>.precursors.with_rt.tsv
<prefix>.Prosit_max7.tsv
<prefix>.MS2PIP_max7.tsv
<prefix>.Union14.tsv
<prefix>.build.summary.tsv
<prefix>.build.params.json
```

## Step-by-step workflow

### 1. FASTA digest

```bash
compactlib digest \
  --fasta human.fasta \
  --output human.precursors.tsv \
  --preset generic-dia
```

### 2. RT prediction

```bash
compactlib predict-rt \
  --precursors human.precursors.tsv \
  --output human.precursors.with_rt.tsv \
  --backend deeplc \
  --rt-output-column PredictedRetentionTime
```

### 3. Prosit prediction

```bash
compactlib predict-chunked \
  --precursors human.precursors.with_rt.tsv \
  --output human.Prosit_max7.tsv \
  --work-dir human.Prosit_max7.work \
  --backend koina \
  --model Prosit_2020_intensity_HCD \
  --preset generic-dia \
  --top-n 7 \
  --include-rt \
  --rt-column PredictedRetentionTime \
  --resume
```

### 4. MS2PIP prediction

```bash
compactlib predict-chunked \
  --precursors human.precursors.with_rt.tsv \
  --output human.MS2PIP_max7.tsv \
  --work-dir human.MS2PIP_max7.work \
  --backend koina \
  --model ms2pip_HCD2021 \
  --preset generic-dia \
  --top-n 7 \
  --include-rt \
  --rt-column PredictedRetentionTime \
  --resume
```

### 5. Union library

```bash
compactlib union \
  --input-a human.Prosit_max7.tsv \
  --input-b human.MS2PIP_max7.tsv \
  --output human.Union14.tsv \
  --top-n 7
```
