# Troubleshooting

## Koina error with modified peptide strings

If Koina fails because peptide sequences contain modification syntax such as `C(UniMod:4)`, use:

```bash
--sequence-source stripped
```

This is the default recommended setting for fixed carbamidomethylated cysteine workflows.

## Selenocysteine-containing peptides

Prosit/MS2PIP prediction backends generally do not support `U`-containing peptide sequences. Do not use `--allow-selenocysteine` for routine Koina prediction unless unsupported amino acids are explicitly dropped before prediction.

## Large output files

By default, predicted libraries are written in a minimal DIA-library format. To include additional metadata, use:

```bash
--include-metadata
```

This can substantially increase file size.

## Interrupted prediction

Use `predict-chunked` or `build` with `--resume`. Re-running the same command will skip completed chunks.

## DeepLC RT scale

If no calibration table is provided, DeepLC output should be interpreted as a predicted RT scale rather than exact chromatographic minutes.
