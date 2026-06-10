# compact-dia-library

`compactlib` is a lightweight Python CLI utility for reproducible construction of compact DIA spectral libraries from existing predicted libraries.

Version `1.0.0` supports post-processing of existing predicted libraries, FASTA digestion into precursor tables, optional variable M oxidation, prediction backends, resumable chunked Koina prediction, minimal DIA-library output by default, and organism-neutral Skyline/DIA-NN-like practical filtering presets.

## Why compact libraries?

Full in silico DIA libraries can contain many fragment transitions per precursor. `compactlib` implements deterministic rules for retaining a limited number of informative transitions per precursor, building model-integrated Union libraries, and generating control libraries.

## FASTA digestion (v0.2)

`compactlib digest` creates a precursor table from FASTA for downstream prediction backends.

```bash
compactlib digest \
  --fasta human.fasta \
  --output human_precursors.tsv \
  --preset generic-dia \
  --c-mod-format unimod
```

The output contains `ModifiedPeptide`, `StrippedPeptide`, `PrecursorCharge`, `PrecursorMz`, `ProteinId`, peptide length and missed-cleavage metadata. Fixed C carbamidomethylation is supported and encoded as `C(UniMod:4)` by default. Variable M oxidation can be enabled with `--variable-oxidation-m`.

`compactlib digest` supports optional precursor candidate filters. The most useful practical filters are precursor m/z range, for example `--min-precursor-mz 400 --max-precursor-mz 1200`, and charge-length rules, for example `--charge-length-rules 4:12`, which removes charge 4 candidates for peptides shorter than 12 amino acids. These filters are optional; omit them if you need exhaustive charge enumeration.

## Skyline/DIA-NN-like filtering presets

`compactlib` provides transparent convenience presets for practical DIA library construction. They are intended to mimic common Skyline/DIA-NN-style restrictions on peptide length, precursor charges, precursor m/z range, fragment ion types and fragment series number. They are not exact locked reproductions of any particular software release; every parameter can be overridden from the command line and is recorded in `params.json`.

List presets:

```bash
compactlib presets
```

Digest with the neutral broad DIA preset:

```bash
compactlib digest \
  --fasta human.fasta \
  --output human_precursors.tsv \
  --preset generic-dia
```

Override a preset parameter explicitly:

```bash
compactlib digest \
  --fasta human.fasta \
  --output human_precursors_windowed.tsv \
  --preset generic-dia \
  --min-precursor-mz 400 \
  --max-precursor-mz 1200 \
  --charges 2,3,4
```

Apply transition-filtering presets during prediction:

```bash
compactlib predict \
  --precursors human_precursors.tsv \
  --output Prosit_max7.tsv \
  --backend koina \
  --model Prosit_2020_intensity_HCD \
  --preset generic-dia \
  --top-n 7
```


Audit note: the Skyline-derived human precursor export used to calibrate `skyline-dia-like` contained 4,182,467 exported rows, 1,232,841 unique stripped peptides, 3,963,681 unique modified-peptide/charge precursors, charges 1-4, peptide length 7-30, and precursor m/z approximately 300-1800. The observed charge/length pattern is explained primarily by the m/z filter; no hard `4:12` charge-length rule is used in the audited Skyline-like preset.

Available presets:

- `diann-like`: trypsin/P, 1 missed cleavage, peptide length 7-30, charges 2/3/4, precursor m/z 400-1200, b/y fragments, fragment charges 1/2, fragment series >=3.
- `skyline-dia-like`: audited Skyline-like precursor export preset, peptide length 7-30, charges 1/2/3/4, precursor m/z 300-1800, recommended top-n 6.
- `skyline-export-like`: explicit alias for the audited Skyline precursor export preset.
- `generic-dia`: neutral broad DIA preset aligned to the audited Skyline-like precursor space: peptide length 7-30, charges 1/2/3/4, precursor m/z 300-1800, recommended top-n 7. Use this for human, yeast, E. coli and other FASTA-derived libraries unless a narrower instrument window is desired.
- `generic-dia-windowed`: narrower DIA-window-focused preset, charges 2/3/4, precursor m/z 400-1200. 
- `generic-dia-strict`: same as `generic-dia-windowed`, plus optional charge-length plausibility rule `4:12`. 




## Retention time prediction with DeepLC (v1.0.0)

`compactlib predict-rt` adds RT predictions to a precursor table before spectral prediction. The command uses DeepLC as an optional backend and converts compactlib peptide notation into DeepLC `seq`/`modifications` input. Fixed carbamidomethylated cysteines are indicated explicitly by default, because DeepLC expects fixed modifications to be provided in the input.

Install the optional dependency:

```bash
pip install -e '.[deeplc]'
```

The current DeepLC releases may require Python >=3.11. If you are using an older Python environment, install a compatible DeepLC release manually, for example the version used in a specific study, and then run the same `compactlib predict-rt` command.

Predict RT for a precursor table:

```bash
compactlib predict-rt \
  --precursors human_precursors.tsv \
  --output human_precursors.with_rt.tsv \
  --backend deeplc \
  --rt-output-column NormalizedRetentionTime
```

With a calibration table:

```bash
compactlib predict-rt \
  --precursors human_precursors.tsv \
  --output human_precursors.with_rt.tsv \
  --backend deeplc \
  --calibration-table calibration_peptides.tsv \
  --calibration-rt-column tr \
  --rt-output-column NormalizedRetentionTime
```

The resulting precursor table can be passed directly into `predict` or `predict-chunked`:

```bash
compactlib predict-chunked \
  --precursors human_precursors.with_rt.tsv \
  --output Prosit_max7.with_rt.tsv \
  --work-dir Prosit_max7.with_rt.work \
  --backend koina \
  --model Prosit_2020_intensity_HCD \
  --preset generic-dia \
  --top-n 7 \
  --include-rt \
  --rt-column NormalizedRetentionTime \
  --resume
```

This design keeps fragment intensity prediction and LC-dependent RT prediction as explicit, auditable steps while still providing an end-to-end FASTA-to-library workflow.

## Installation

```bash
cd compact-dia-library
pip install -e .
```

For parquet support:

```bash
pip install -e '.[parquet]'
```

For tests:

```bash
pip install -e '.[dev]'
pytest -q
```

## Required input columns

Input libraries must contain the following columns:

```text
ModifiedPeptide
PrecursorCharge
PrecursorMz
ProductMz
LibraryIntensity
FragmentType
FragmentSeriesNumber
FragmentCharge
```

All additional columns are preserved in the output library.

## Definitions

Precursor key:

```text
ModifiedPeptide + PrecursorCharge
```

Transition key:

```text
ModifiedPeptide + PrecursorCharge + FragmentType + FragmentSeriesNumber + FragmentCharge
```

Tie-breakers for deterministic sorting:

```text
FragmentType, FragmentSeriesNumber, FragmentCharge, ProductMz
```

## Commands

### max-n

Keep up to `N` most intense transitions per precursor:

```bash
compactlib max-n \
  --input Full_prosit_library.csv.gz \
  --output Prosit_max7.tsv \
  --top-n 7
```

### union

Build Union2N library by merging top-N transitions from two libraries:

```bash
compactlib union \
  --input-a Full_prosit_library.csv.gz \
  --input-b Full_ms2pip_library.csv.gz \
  --output Union14.tsv \
  --top-n 7
```

For `--top-n 7`, the output is a Union14-style library. Duplicate transitions are removed by transition annotation.

### consensus

Build consensus top-N library using rank-based score:

```text
score = 1/rank_a + 1/rank_b
```

```bash
compactlib consensus \
  --input-a Full_prosit_library.csv.gz \
  --input-b Full_ms2pip_library.csv.gz \
  --output Consensus_max7.tsv \
  --top-n 7
```

### random-n

Randomly select up to N transitions per precursor using a fixed seed:

```bash
compactlib random-n \
  --input Full_prosit_library.csv.gz \
  --output Prosit_random7.tsv \
  --top-n 7 \
  --seed 42
```

### reverse-intensity

Invert intensity pattern within each precursor while retaining the same transitions:

```bash
compactlib reverse-intensity \
  --input Prosit_max7.tsv \
  --output Prosit_reverse7.tsv
```

## Output files

For every output library, `compactlib` writes:

```text
output_library.tsv
output_library.summary.tsv
output_library.params.json
```

The summary contains the number of input/output precursors and transitions, transition counts per precursor, duplicate counts, and output file size. The JSON file records command-line parameters for reproducibility.

## Notes

- `max-n` keeps **up to N** transitions per precursor. Precursors with fewer than N transitions are retained.
- The implementation is deterministic except for `random-n`, which is reproducible with a fixed seed.
- Version 1.0.0 supports post-processing of existing libraries, FASTA digestion into precursor tables, optional precursor candidate filters, variable M oxidation expansion, DeepLC RT prediction integration, a predictor backend interface, and Koina-based Prosit/MS2PIP prediction adapters.

## Progress logging

`compactlib` prints progress logs for the major processing stages: reading input libraries, selecting or combining transitions, writing the output library, writing summary files, and writing parameters. When `rich` is available, a terminal spinner is shown during long stages.

Disable the spinner but keep logs:

```bash
compactlib max-n --input Full_prosit_library.csv.gz --output Prosit_max7.tsv --top-n 7 --no-progress
```

Run silently:

```bash
compactlib max-n --input Full_prosit_library.csv.gz --output Prosit_max7.tsv --top-n 7 --quiet
```

Each summary file also contains `compactlib_elapsed_sec`, the total wall-clock time spent inside the library construction utility.

## Toy example

A small toy example is provided in `examples/` to inspect algorithm behavior by eye:

```bash
bash examples/04_toy_demo.sh
```

It creates `toy_out/` with outputs for `max-n`, `union`, `consensus`, `random-n`, and `reverse-intensity`, then prints the selected transitions. The example uses `--top-n 3` so that the result is easy to inspect manually.

## Duplicate policy for two-input commands

For `union`, duplicate transitions are removed using the annotation-based transition key. If the same transition is present in both input libraries, the output row is taken from input A (`duplicate_policy = keep_input_a`).

For `consensus`, transition ranking is based on `1/rank_a + 1/rank_b`. If a selected transition is present in both input libraries, metadata and intensity values are taken from input A (`duplicate_policy = prefer_input_a`).

## Variable M oxidation in digest (v0.2.1)

`compactlib digest` can expand methionine oxidation as a variable modification:

```bash
compactlib digest \
  --fasta human.fasta \
  --output human_precursors_mox.tsv \
  --enzyme trypsin-p \
  --missed-cleavages 1 \
  --min-length 7 \
  --max-length 30 \
  --charges 2,3,4 \
  --variable-oxidation-m \
  --max-variable-mods 1 \
  --m-mod-format unimod
```

The unmodified peptidoform is always retained. With `--variable-oxidation-m`, all combinations with 1..`max-variable-mods` oxidized methionines are also generated. Oxidation is encoded as `M(UniMod:35)` by default. The precursor table includes:

```text
VariableOxidationM
OxidationMPositions
NVariableMods
```

`OxidationMPositions` uses 1-based semicolon-separated positions in the stripped peptide.

## Predictor interface (v0.2.1)

`compactlib predict` converts a precursor table into a transition library using a predictor backend interface. Version 0.2.1 includes a deterministic `mock` backend for tests and examples only:

```bash
compactlib predict \
  --precursors human_precursors_mox.tsv \
  --output human_mock_max7.tsv \
  --backend mock \
  --fragment-types b,y \
  --fragment-charges 1 \
  --top-n 7
```

The mock backend is not a scientific Prosit/MS2PIP replacement. It exists to validate the end-to-end pipeline and CLI behavior without external services. Koina/Prosit/MS2PIP adapters are planned for later versions.

## Koina prediction backend, v0.2.2

`compactlib predict` supports an optional Koina/koinapy backend for generating transition libraries from precursor tables produced by `compactlib digest`.

Install optional dependency:

```bash
pip install -e '.[koina]'
```

Example for Prosit:

```bash
compactlib predict \
  --precursors human_precursors.tsv \
  --output Prosit_max7.tsv \
  --backend koina \
  --model Prosit_2020_intensity_HCD \
  --server-url koina.wilhelmlab.org:443 \
  --ssl \
  --collision-energy 30 \
  --batch-size 1024 \
  --sequence-source stripped \
  --top-n 7
```

Example for MS2PIP:

```bash
compactlib predict \
  --precursors human_precursors.tsv \
  --output MS2PIP_max7.tsv \
  --backend koina \
  --model ms2pip_HCD2021 \
  --server-url koina.wilhelmlab.org:443 \
  --ssl \
  --collision-energy 30 \
  --batch-size 1024 \
  --sequence-source stripped \
  --top-n 7
```

For locally hosted Koina, use for example:

```bash
--server-url localhost:8500 --no-ssl
```

Koina model schemas can differ. Therefore input and output column names are configurable:

```bash
--koina-sequence-input peptide_sequences
--koina-charge-input precursor_charges
--koina-ce-input collision_energies
--koina-intensity-col intensities
--koina-annotation-col annotation
--koina-mz-col mz
```

If the model output does not contain fragment m/z values, `compactlib` calculates b/y ProductMz values from `StrippedPeptide`, fixed carbamidomethyl C, and optional M oxidation positions.

## Selenocysteine handling

By default, `compactlib digest` removes peptides containing non-standard amino acids.
To reproduce Skyline precursor exports that retain selenocysteine-containing peptides, use:

```bash
compactlib digest \
  --fasta human.fasta \
  --output human_precursors_with_u.tsv \
  --preset generic-dia \
  --allow-selenocysteine
```

This keeps `U`-containing peptides and calculates their precursor masses using the
selenocysteine residue mass. This mode is intended mainly for Skyline-compatible
precursor-space audits. Most Prosit/MS2PIP/Koina prediction models do not support
`U`-containing peptide sequences. Therefore, `compactlib predict --backend koina`
stops with a clear error if such precursors are present. To explicitly remove them
before prediction, pass:

```bash
compactlib predict \
  --precursors human_precursors_with_u.tsv \
  --output Prosit_max7.tsv \
  --backend koina \
  --model Prosit_2020_intensity_HCD \
  --drop-unsupported-aa \
  --top-n 7
```


## Predicted library output columns, v0.2.10

By default, `compactlib predict` and `compactlib predict-chunked` write a minimal DIA transition library instead of repeating all precursor/digest/prediction metadata on every transition row.

Default output columns:

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
```

This keeps full-proteome libraries much smaller. To retain all metadata columns for debugging or audit trails, add:

```bash
--include-metadata
```

For example:

```bash
compactlib predict-chunked \
  --precursors human_precursors.tsv \
  --output Prosit_max7.fullhuman.with_metadata.tsv \
  --work-dir Prosit_max7.with_metadata.chunked_work \
  --backend koina \
  --model Prosit_2020_intensity_HCD \
  --top-n 7 \
  --include-metadata
```

## Retention time support, v0.2.13

`compactlib` does not implement an internal RT predictor in v0.2.13. Instead, it supports RT values predicted externally, for example by DeepLC. This keeps LC-method-dependent RT prediction separate from fragment intensity prediction while preserving a reproducible DIA library construction workflow.

There are two supported RT workflows.

### 1. Pass RT through during prediction

If the precursor table already contains an RT column, add it to the predicted transition library using `--include-rt`:

```bash
compactlib predict-chunked \
  --precursors human_precursors_with_rt.tsv \
  --output Prosit_max7.with_rt.tsv \
  --work-dir Prosit_max7.with_rt.work \
  --backend koina \
  --model Prosit_2020_intensity_HCD \
  --preset generic-dia \
  --top-n 7 \
  --include-rt \
  --rt-column Tr_recalibrated \
  --rt-output-column NormalizedRetentionTime \
  --resume
```

The minimal output then contains the standard 10 transition-library columns plus `NormalizedRetentionTime`.

### 2. Merge RT into an existing library

Use `add-rt` to attach external precursor-level RT predictions to an already generated library:

```bash
compactlib add-rt \
  --library Prosit_max7.fullhuman.tsv \
  --rt-table deeplc_predictions.tsv \
  --output Prosit_max7.fullhuman.with_rt.tsv \
  --key ModifiedPeptide,PrecursorCharge \
  --rt-column Tr_recalibrated \
  --rt-output-column NormalizedRetentionTime \
  --min-match-rate 0.99
```

The command writes sidecar summary and parameter files reporting the number of library precursors, RT-table keys, matched/unmatched precursors, duplicate RT keys, conflicting RT values, and precursor-level match rate.

## Robust full-proteome prediction with persistent chunks

For large precursor tables, use `predict-chunked` instead of manually splitting files.
The command creates persistent precursor chunks, predicts each chunk, writes a `.done.json`
marker after successful completion, and merges the final library only when all chunks are complete.
Re-running the same command with `--resume` skips completed chunks and continues from the first missing chunk.

Example for full-human Prosit max7:

```bash
compactlib predict-chunked \
  --precursors human_precursors.tsv \
  --output Prosit_max7.fullhuman.tsv \
  --work-dir Prosit_max7.chunked_work \
  --chunk-size 100000 \
  --backend koina \
  --model Prosit_2020_intensity_HCD \
  --server-url koina.wilhelmlab.org:443 \
  --ssl \
  --collision-energy 30 \
  --batch-size 512 \
  --sequence-source stripped \
  --preset generic-dia \
  --top-n 7 \
  --resume \
  --verbose
```

Example for MS2PIP max7:

```bash
compactlib predict-chunked \
  --precursors human_precursors.tsv \
  --output MS2PIP_max7.fullhuman.tsv \
  --work-dir MS2PIP_max7.chunked_work \
  --chunk-size 100000 \
  --backend koina \
  --model ms2pip_HCD2021 \
  --server-url koina.wilhelmlab.org:443 \
  --ssl \
  --collision-energy 30 \
  --batch-size 512 \
  --sequence-source stripped \
  --preset generic-dia \
  --top-n 7 \
  --resume \
  --verbose
```

The work directory contains:

```text
precursor_chunks/      # persistent input chunks
predicted_chunks/      # completed predicted library chunks
chunk_summaries/       # per-chunk summary TSV files
chunk_params/          # per-chunk parameter JSON files
chunks_manifest.tsv    # chunk status table
predict_chunked.params.json
```

If a run fails because of a network or Koina issue, run the same command again. Completed chunks are skipped by default. Use `--force` only when you intentionally want to recompute all chunks.

## One-command build workflow

`compactlib build` runs the full FASTA-to-library workflow in one command:

```bash
compactlib build \
  --fasta proteome.fasta \
  --output-dir compactlib_out \
  --prefix proteome \
  --preset generic-dia \
  --top-n 7 \
  --chunk-size 100000 \
  --prediction-backend koina \
  --prosit-model Prosit_2020_intensity_HCD \
  --ms2pip-model ms2pip_HCD2021 \
  --server-url koina.wilhelmlab.org:443 \
  --ssl \
  --collision-energy 30 \
  --batch-size 512 \
  --rt-backend deeplc \
  --rt-output-column PredictedRetentionTime \
  --resume
```

The command creates minimal DIA-library outputs by default:

- `<prefix>.Prosit_max7.tsv`
- `<prefix>.MS2PIP_max7.tsv`
- `<prefix>.Union14.tsv`

For smoke tests, add `--sample-size 1000` and reduce `--chunk-size`, for example `--chunk-size 250`.
