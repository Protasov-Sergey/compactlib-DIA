# Changelog

## v1.0.3 - 2026-06-11

Publication-readiness update.

- Added explicit MIT `LICENSE` file.
- Replaced generic package author metadata with the project author.
- Added GitHub Actions CI workflow for package tests.
- Polished README wording for RT usage.
- No changes to the core `compactlib build` workflow.


## v1.0.2

- Added `compactlib build`, a one-command FASTA-to-library workflow.
- The build workflow performs digest, optional RT prediction, resumable Prosit/MS2PIP max-N prediction, and optional Union2N construction.
- Added optional `--sample-size` for small black-box smoke tests while preserving predictor-safe stratified precursor sampling.
- Build outputs use minimal DIA-library columns by default and can pass RT into Prosit/MS2PIP/Union outputs.

## v1.0.1

- Fixed DeepLC RT prediction without a calibration table by calling DeepLC with `calibrate=False` when needed.
- Added DeepLC RT prediction metadata field `deeplc_make_preds_calibrate_arg` to summary output.


## v1.0.0

- Added `compactlib predict-rt` for DeepLC-based retention time prediction.
- Added compactlib-to-DeepLC peptide/modification conversion (`seq`, `modifications`).
- Added optional DeepLC calibration table support.
- Added optional `.[deeplc]` dependency group.
- Kept external RT merging via `add-rt` and RT passthrough in `predict`/`predict-chunked`.

## v0.2.13

### Added
- Added `compactlib add-rt` to merge externally predicted precursor-level RT values into transition libraries.
- Added `--include-rt`, `--rt-column`, and `--rt-output-column` to `compactlib predict` and `compactlib predict-chunked` for RT passthrough from precursor tables into predicted libraries.
- Added RT merge summaries with matched/unmatched precursor counts, duplicate RT-key statistics, and match-rate validation.

### Notes
- v0.2.13 does not implement an internal RT predictor. It supports RT values produced externally, for example by DeepLC, and keeps LC-method-specific RT prediction separate from fragment intensity prediction.

## v0.2.12

### Changed
- Removed the legacy human-specific preset aliases `project-human-dia`, `project-human-dia-windowed`, and `project-human-dia-strict` from the final public preset list.
- The recommended organism-neutral presets are now `generic-dia`, `generic-dia-windowed`, and `generic-dia-strict`.
- CLI help, README, and tests were updated to avoid human-specific preset names.

## v0.2.11

### Changed
- Added organism-neutral presets `generic-dia`, `generic-dia-windowed`, and `generic-dia-strict` for cross-organism DIA library generation.

## v0.2.10

### Changed
- `compactlib predict` and `compactlib predict-chunked` now write minimal DIA library columns by default instead of repeating precursor/digest/prediction metadata on every transition row.
- Default predicted output columns are: `ModifiedPeptide`, `StrippedPeptide`, `PrecursorCharge`, `PrecursorMz`, `ProductMz`, `LibraryIntensity`, `FragmentType`, `FragmentSeriesNumber`, `FragmentCharge`, and `ProteinId`.

### Added
- Added `--include-metadata/--no-metadata` to `predict` and `predict-chunked`. Use `--include-metadata` to keep the previous full metadata output for debugging/auditing.
- `predict-chunked` final merge can now re-merge older full-metadata chunks into a minimal final output, so completed v0.2.9 chunks do not need to be recomputed just to reduce the final library columns.


## v0.2.9

- Added `compactlib predict-chunked` for resumable full-proteome prediction using persistent outer chunks.
- Added per-chunk `.done.json` markers, chunk summaries, chunk params, a manifest table, and final merge after all chunks are complete.
- Added atomic temporary output handling for predicted chunks and final chunk merge.

## v0.2.8

- Added Koina long-format output adapter for Prosit/MS2PIP-style predictions where each row corresponds to a fragment and input peptide/charge columns are repeated.
- Kept support for array-valued one-row-per-precursor Koina outputs.
- Added tests for long-format Koina conversion without requiring network access.


## v0.2.7

- Changed Koina default sequence source to `stripped`, which is appropriate for `Prosit_2020_intensity_HCD` and `ms2pip_HCD2021` when cysteine carbamidomethylation is fixed/implicit.
- Added a preflight error explaining why modified strings such as `C(UniMod:4)` fail for Koina inputs named `peptide_sequences`.
- Added a guard against silently dropping variable-modification information when `--sequence-source stripped` is used on modified precursor tables.


## v0.2.6

- Added optional selenocysteine support in FASTA digestion via `--allow-selenocysteine`.
- Added U residue mass for Skyline-compatible precursor m/z calculation.
- Added digest summary fields for U-containing peptides/precursors.
- Added predictor-side guardrails: real prediction backends stop on U-containing precursors unless `--drop-unsupported-aa` is provided.
- Documented that Prosit/MS2PIP/Koina models generally do not support selenocysteine-containing peptide sequences.

## v0.2.5

### Changed
- Updated `skyline-dia-like` and `project-human-dia` digest presets using the audited Skyline-derived human precursor export.
- `skyline-dia-like`, `skyline-export-like`, and `project-human-dia` now use peptide length 7-30, precursor charges 1,2,3,4, precursor m/z 300-1800, one missed cleavage, and no hard charge-length rule.
- Kept narrower DIA-window-focused settings as `project-human-dia-windowed`.
- Kept the optional charge-length plausibility heuristic as `project-human-dia-strict`.

### Added
- `skyline-export-like` preset as an explicit alias for the audited Skyline precursor export.
- Tests for the audited Skyline-like preset values.

### Notes
- The Skyline audit showed 4,182,467 exported rows, 1,232,841 unique stripped peptides, 3,963,681 unique modified-peptide/charge precursors, charges 1-4, peptide length 7-30, and precursor m/z approximately 300-1800.
- Charge 4 candidates started at peptide length 8; therefore the audited Skyline-like preset does not apply the earlier `4:12` hard rule.

## v0.2.4

### Added
- Skyline/DIA-NN-like practical filtering presets: `diann-like`, `skyline-dia-like`, `project-human-dia`, and `project-human-dia-strict`.
- `compactlib presets` command to print available preset definitions.
- `--preset` option for `compactlib digest` and `compactlib predict`.
- Post-prediction transition filtering by fragment type, fragment charge, fragment series number, and optional product m/z.
- Preset metadata in `summary.tsv` and `params.json`.
- Example `examples/10_presets_demo.sh`.

### Notes
- Presets are transparent convenience configurations, not exact locked reproductions of a specific Skyline or DIA-NN release. All preset parameters can be overridden from the command line.

## v0.2.3

### Added
- Precursor m/z filtering in `compactlib digest` via `--min-precursor-mz` and `--max-precursor-mz`.
- Optional charge-length plausibility rules via `--charge-length-rules`, e.g. `4:12,5:16`.
- Digest summary fields for precursor filtering: `n_precursors_before_precursor_filters`, `n_precursors_filtered_mz`, and `n_precursors_filtered_charge_length`.

### Notes
- These filters are optional and are intended to mimic practical Skyline/DIA-style precursor candidate constraints such as peptide length, charge range, and instrument m/z range.

## v0.2.2

### Added
- Optional Koina/koinapy prediction backend for `compactlib predict`.
- Configurable Koina model name, server URL, SSL, collision energy, batch size, input column names, and output column names.
- Examples for Prosit and MS2PIP prediction through Koina.

### Notes
- `mock` backend remains the default for tests and offline examples.
- Koina backend requires optional dependency: `pip install -e '.[koina]'`.


## v0.2.0

### Added
- `compactlib digest` command for FASTA digestion into precursor tables.
- Monoisotopic precursor m/z calculation with fixed C carbamidomethylation.
- Trypsin and trypsin/P digestion modes.
- Toy FASTA digestion example.


## v0.1.2

### Fixed
- Removed incomplete single-input summary fields from two-input commands (`union`, `consensus`).
- Added explicit `n_precursors_input_a/b` and `n_transitions_input_a/b` fields for two-input commands.
- Added precursor overlap statistics to `consensus` summaries.
- Added duplicate policies to `union` and `consensus` summaries.

### Added
- Toy example files and a small demo script for visual sanity checks.
- `.gitignore` for Python/Jupyter/example outputs.

## v0.1.1

### Added
- Progress/logging output for all CLI commands.
- `--verbose/--quiet` and `--progress/--no-progress` options.
- `compactlib_elapsed_sec` in summary files.

## v0.1.0

### Added
- `max-n`, `union`, `consensus`, `random-n`, and `reverse-intensity` commands.
- Summary and params sidecar files for every output library.

## v0.2.1

### Added
- Variable methionine oxidation expansion in `compactlib digest` via `--variable-oxidation-m`.
- `OxidationMPositions`, `VariableOxidationM`, and `NVariableMods` columns in precursor tables.
- Predictor backend interface under `compactlib.predictors`.
- Deterministic `mock` predictor backend for testing FASTA/precursor-to-library workflows without external services.
- `compactlib predict` command for precursor table -> transition library generation, with optional `--top-n` compaction.

### Notes
- The mock predictor is for testing and examples only. Scientific Prosit/MS2PIP/Koina adapters are planned for later versions.

## v0.2.9

- Added `compactlib predict-chunked` for robust full-proteome prediction workflows.
- The command splits a precursor table into persistent outer chunks, predicts each chunk, writes per-chunk completion markers, resumes completed chunks on repeated runs, and merges the final output only after all chunks are complete.
- Added chunk manifests, per-chunk summaries/parameters, and final chunked summary/parameter sidecars.
