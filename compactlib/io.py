from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd

REQUIRED_COLUMNS = [
    "ModifiedPeptide",
    "PrecursorCharge",
    "PrecursorMz",
    "ProductMz",
    "LibraryIntensity",
    "FragmentType",
    "FragmentSeriesNumber",
    "FragmentCharge",
]

PRECURSOR_KEY = ["ModifiedPeptide", "PrecursorCharge"]
TRANSITION_KEY = [
    "ModifiedPeptide",
    "PrecursorCharge",
    "FragmentType",
    "FragmentSeriesNumber",
    "FragmentCharge",
]

NUMERIC_COLUMNS = [
    "PrecursorCharge",
    "PrecursorMz",
    "ProductMz",
    "LibraryIntensity",
    "FragmentSeriesNumber",
    "FragmentCharge",
]


def validate_columns(df: pd.DataFrame, required: Iterable[str] = REQUIRED_COLUMNS) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            "Missing required columns: " + ", ".join(missing) +
            f". Available columns: {list(df.columns)}"
        )


def coerce_numeric_commas(df: pd.DataFrame, cols: Iterable[str] = NUMERIC_COLUMNS) -> pd.DataFrame:
    """Convert numeric columns, accepting decimal commas and spaces.

    DIA/Skyline-like tables sometimes contain values such as '1000,002558'.
    """
    out = df.copy()
    for c in cols:
        if c not in out.columns:
            continue
        if pd.api.types.is_numeric_dtype(out[c]):
            continue
        out[c] = (
            out[c]
            .astype(str)
            .str.replace("\u00A0", "", regex=False)
            .str.replace(" ", "", regex=False)
            .str.replace(",", ".", regex=False)
        )
        out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def normalize_required_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Validate required columns and coerce core numeric fields."""
    validate_columns(df)
    out = coerce_numeric_commas(df)

    # FragmentType is used as deterministic tie-breaker; keep original value in output.
    # Here we only ensure missing values do not break sorting/grouping.
    out["FragmentType"] = out["FragmentType"].astype(str)

    # Drop unusable rows. We keep rows with missing non-core extra columns.
    core_not_null = [
        "ModifiedPeptide", "PrecursorCharge", "ProductMz", "LibraryIntensity",
        "FragmentType", "FragmentSeriesNumber", "FragmentCharge"
    ]
    out = out.dropna(subset=core_not_null).copy()

    # Use integer-like charges/series where possible, but do not force object dtypes in output.
    for c in ["PrecursorCharge", "FragmentSeriesNumber", "FragmentCharge"]:
        out[c] = out[c].astype(int)

    return out


def read_table(path: str | Path) -> pd.DataFrame:
    """Read a generic .tsv, .csv, .csv.gz, .tsv.gz or .parquet table without spectral-library validation."""
    path = Path(path).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"Input file does not exist: {path}")

    name = path.name.lower()
    if name.endswith(".parquet"):
        try:
            return pd.read_parquet(path)
        except ImportError as e:
            raise ImportError(
                "Reading parquet requires pyarrow or fastparquet. Install with: "
                "pip install compact-dia-library[parquet]"
            ) from e
    if name.endswith(".tsv") or name.endswith(".tsv.gz"):
        return pd.read_csv(path, sep="\t")
    if name.endswith(".csv") or name.endswith(".csv.gz"):
        return pd.read_csv(path)
    raise ValueError(
        f"Unsupported input format for {path}. Supported: .tsv, .tsv.gz, .csv, .csv.gz, .parquet"
    )


def read_library(path: str | Path) -> pd.DataFrame:
    """Read .tsv, .csv, .csv.gz, .tsv.gz or .parquet spectral library."""
    path = Path(path).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"Input file does not exist: {path}")

    name = path.name.lower()
    if name.endswith(".parquet"):
        try:
            df = pd.read_parquet(path)
        except ImportError as e:
            raise ImportError(
                "Reading parquet requires pyarrow or fastparquet. Install with: "
                "pip install compact-dia-library[parquet]"
            ) from e
    elif name.endswith(".tsv") or name.endswith(".tsv.gz"):
        df = pd.read_csv(path, sep="\t")
    elif name.endswith(".csv") or name.endswith(".csv.gz"):
        df = pd.read_csv(path)
    else:
        raise ValueError(
            f"Unsupported input format for {path}. Supported: .tsv, .tsv.gz, .csv, .csv.gz, .parquet"
        )

    return normalize_required_columns(df)


def write_library(df: pd.DataFrame, path: str | Path) -> None:
    """Write library to .tsv, .tsv.gz, .csv, .csv.gz or .parquet."""
    path = Path(path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    name = path.name.lower()

    if name.endswith(".parquet"):
        try:
            df.to_parquet(path, index=False)
        except ImportError as e:
            raise ImportError(
                "Writing parquet requires pyarrow or fastparquet. Install with: "
                "pip install compact-dia-library[parquet]"
            ) from e
    elif name.endswith(".tsv") or name.endswith(".tsv.gz"):
        df.to_csv(path, sep="\t", index=False)
    elif name.endswith(".csv") or name.endswith(".csv.gz"):
        df.to_csv(path, index=False)
    else:
        raise ValueError(
            f"Unsupported output format for {path}. Supported: .tsv, .tsv.gz, .csv, .csv.gz, .parquet"
        )


def sidecar_paths(output_path: str | Path) -> tuple[Path, Path]:
    """Return (summary.tsv, params.json) paths derived from output library path."""
    path = Path(output_path).expanduser()
    name = path.name

    # Strip common compound suffixes to avoid output.tsv.gz.summary.tsv.
    for suffix in [".tsv.gz", ".csv.gz", ".tsv", ".csv", ".parquet"]:
        if name.lower().endswith(suffix):
            stem = name[: -len(suffix)]
            break
    else:
        stem = path.stem

    return path.with_name(f"{stem}.summary.tsv"), path.with_name(f"{stem}.params.json")


def write_table(df: pd.DataFrame, path: str | Path) -> None:
    """Write a generic table to .tsv, .tsv.gz, .csv, .csv.gz or .parquet.

    Unlike ``write_library``, this function does not validate spectral-library
    transition columns and is used for precursor tables produced by ``digest``.
    """
    path = Path(path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    name = path.name.lower()

    if name.endswith(".parquet"):
        try:
            df.to_parquet(path, index=False)
        except ImportError as e:
            raise ImportError(
                "Writing parquet requires pyarrow or fastparquet. Install with: "
                "pip install compact-dia-library[parquet]"
            ) from e
    elif name.endswith(".tsv") or name.endswith(".tsv.gz"):
        df.to_csv(path, sep="\t", index=False)
    elif name.endswith(".csv") or name.endswith(".csv.gz"):
        df.to_csv(path, index=False)
    else:
        raise ValueError(
            f"Unsupported output format for {path}. Supported: .tsv, .tsv.gz, .csv, .csv.gz, .parquet"
        )
