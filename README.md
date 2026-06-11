# compactlib-DIA

<<<<<<< HEAD
**compactlib-DIA** is a utility for building compact in silico DIA spectral libraries from protein FASTA files using Prosit, MS2PIP, DeepLC and Koina-compatible prediction backends.
=======
[![CI](https://github.com/Protasov-Sergey/compactlib-DIA/actions/workflows/ci.yml/badge.svg)](https://github.com/Protasov-Sergey/compactlib-DIA/actions/workflows/ci.yml)

**compactlib-DIA** is a command-line utility for building compact in silico DIA spectral libraries from protein FASTA files using Prosit, MS2PIP, DeepLC and Koina-compatible prediction backends.
>>>>>>> 1234a03 (Prepare compactlib-DIA for publication release)

The main workflow produces:

- Prosit maxN spectral library;
- MS2PIP maxN spectral library;
- Union2N integrated spectral library;
- optional DeepLC-predicted retention times;
- reproducible summary and parameter sidecar files.

The recommended one-command entry point is:

```bash
compactlib build
```

---

## Quick start

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

For a small test run:

```bash
compactlib build \
  --fasta human.fasta \
  --output-dir compactlib_test_1k \
  --prefix human_1k \
  --sample-size 1000 \
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
  --chunk-size 250 \
  --batch-size 128 \
  --resume \
  --verbose
```

---

## Output files

A typical `compactlib build` run creates:

```text
human.Prosit_max7.tsv
human.MS2PIP_max7.tsv
human.Union14.tsv
human.precursors.all.tsv
human.precursors.with_rt.tsv
human.build.summary.tsv
human.build.params.json
```

The final spectral libraries are written in a minimal DIA-library format by default.

### Default output columns

```text
ModifiedPeptide
StrippedPeptide
PrecursorCharge
PrecursorMz
ProductMz
LibraryIntensity
FragmentType
FragmentSeriesNumber
FragmentCharge
ProteinId
PredictedRetentionTime
```

If RT prediction is not requested, the RT column is omitted.

---

## Installation

Clone the repository:

```bash
git clone git@github.com:Protasov-Sergey/compactlib-DIA.git
cd compactlib-DIA
```

Create an environment:

```bash
conda create -n compactlib-dia python=3.11 -y
conda activate compactlib-dia
```

Install with Koina and DeepLC support:

```bash
pip install -e '.[dev,koina,deeplc]'
```

For workflows reproducing older DeepLC-based analyses, a specific DeepLC version can be installed separately:

```bash
pip install deeplc==1.1.2
pip install -e '.[dev,koina]'
```

Run tests:

```bash
pytest -q
```

---

## Main commands

```text
compactlib build             # one-command FASTA-to-library workflow
compactlib digest            # FASTA digestion and precursor table generation
compactlib predict-rt        # DeepLC-based RT prediction for precursor tables
compactlib predict           # prediction for a precursor table
compactlib predict-chunked   # resumable chunked prediction for large precursor tables
compactlib add-rt            # merge externally predicted RT into an existing library
compactlib max-n             # select top-N transitions from an existing library
compactlib union             # build an integrated Union2N library from two libraries
compactlib consensus         # build a rank-based consensus library
compactlib random-n          # random transition selection control
compactlib reverse-intensity # reverse intensity-pattern control
compactlib presets           # show available filtering presets
```

---

## Presets

Recommended default:

```bash
--preset generic-dia
```

Available presets:

```text
generic-dia
generic-dia-windowed
generic-dia-strict
diann-like
skyline-dia-like
skyline-export-like
```

The default `generic-dia` preset is organism-neutral and can be used for human, yeast, bacteria and other FASTA-derived proteomes.

---

## Minimal step-by-step workflow

The one-command `build` workflow is recommended for most users. The same process can also be run step by step:

```bash
compactlib digest \
  --fasta human.fasta \
  --output human.precursors.tsv \
  --preset generic-dia
```

```bash
compactlib predict-rt \
  --precursors human.precursors.tsv \
  --output human.precursors.with_rt.tsv \
  --backend deeplc \
  --rt-output-column PredictedRetentionTime
```

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

---

## Limitations

- DeepLC RT values are model predictions and may require calibration for a specific LC method. You can use own RT values or RT predictor for build the spectral library.
- Selenocysteine-containing peptides can be retained during FASTA digestion for Skyline-compatible precursor-space auditing, but are generally not supported by Prosit/MS2PIP prediction backends.
- Variable modifications should be used carefully. The default workflow is designed for fixed carbamidomethylated cysteine and no variable modifications.
- Very large proteomes should be processed using `predict-chunked` or `build` with `--resume`.

---

## Documentation

Additional documentation is available in `docs/`:

```text
docs/full_workflow.md
docs/presets.md
docs/rt_deeplc.md
docs/chunked_prediction.md
docs/postprocessing.md
docs/troubleshooting.md
docs/development_notes.md
```

---

## Citation

If you use `compactlib-DIA`, please cite this repository and the corresponding publication when available.

```text
compactlib-DIA: compact in silico DIA spectral library generation using Prosit, MS2PIP, DeepLC and Koina.
```
