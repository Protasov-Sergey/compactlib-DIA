# Chunked prediction and resume

Large FASTA-derived precursor tables should be processed using `predict-chunked` or `build` with `--resume`.

## Example

```bash
compactlib predict-chunked \
  --precursors human.precursors.with_rt.tsv \
  --output human.Prosit_max7.tsv \
  --work-dir human.Prosit_max7.work \
  --chunk-size 100000 \
  --backend koina \
  --model Prosit_2020_intensity_HCD \
  --preset generic-dia \
  --top-n 7 \
  --include-rt \
  --rt-column PredictedRetentionTime \
  --resume \
  --verbose
```

## Work directory

The work directory contains persistent chunk files, predicted chunk outputs, summaries and marker files.

```text
human.Prosit_max7.work/
├── precursor_chunks/
├── predicted_chunks/
├── chunk_summaries/
├── chunk_params/
├── chunks_manifest.tsv
└── predict_chunked.params.json
```

A chunk is considered complete only after a corresponding `.done.json` marker is written.

## Resume behavior

If the process stops, run the same command again with `--resume`. Completed chunks will be skipped and missing chunks will be predicted.

## Recommended chunk sizes

```text
small tests: 250-2500 precursors
large proteomes: 100000 precursors
```
