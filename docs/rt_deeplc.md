# Retention time prediction with DeepLC

`compactlib-DIA` supports RT prediction through DeepLC.

## Predict RT for a precursor table

```bash
compactlib predict-rt \
  --precursors human.precursors.tsv \
  --output human.precursors.with_rt.tsv \
  --backend deeplc \
  --rt-output-column PredictedRetentionTime
```

## Use RT in spectral library prediction

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

## Calibration

DeepLC predictions may require calibration for a specific LC setup. If no calibration table is provided, the output should be interpreted as a model-predicted RT scale rather than exact chromatographic minutes.

## Add externally predicted RT

If RT values were generated outside `compactlib-DIA`, use:

```bash
compactlib add-rt \
  --library human.Prosit_max7.tsv \
  --rt-table deeplc_predictions.tsv \
  --output human.Prosit_max7.with_rt.tsv \
  --key ModifiedPeptide,PrecursorCharge \
  --rt-column Tr_recalibrated \
  --rt-output-column PredictedRetentionTime
```
