# Presets

`compactlib-DIA` uses presets to define practical DIA-like precursor and transition filters.

## Recommended default

```bash
--preset generic-dia
```

`generic-dia` is organism-neutral and can be used for human, yeast, bacterial and other proteomes.

## Available presets

```text
generic-dia
generic-dia-windowed
generic-dia-strict
diann-like
skyline-dia-like
skyline-export-like
```

## generic-dia

Broad Skyline-like precursor-space preset.

```text
enzyme: trypsin-p
missed cleavages: 1
peptide length: 7-30
precursor charges: 1,2,3,4
precursor m/z: 300-1800
fragment types: b,y
fragment charges: 1,2,3
minimum fragment series number: 3
recommended top-n: 7
```

## generic-dia-windowed

More practical DIA-window preset.

```text
precursor charges: 2,3,4
precursor m/z: 400-1200
recommended top-n: 7
```

## generic-dia-strict

A stricter variant of `generic-dia-windowed` with an additional charge-length plausibility rule.

```text
charge 4 allowed only for peptide length >= 12
```

This is an optional heuristic, not a universal rule.

## Checking presets

```bash
compactlib presets
```
