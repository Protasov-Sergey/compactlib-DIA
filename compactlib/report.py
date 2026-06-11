from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pandas as pd

from . import __version__
from .io import PRECURSOR_KEY


def _n_precursors(df: pd.DataFrame) -> int:
    return int(df[PRECURSOR_KEY].drop_duplicates().shape[0])


def library_summary(
    df_input: pd.DataFrame | None,
    df_output: pd.DataFrame,
    output_path: str | Path,
    mode: str,
    extra: dict[str, Any] | None = None,
) -> pd.DataFrame:
    """Create one-row summary for compact library construction.

    For one-input commands, the summary contains n_precursors_input and
    n_transitions_input. For two-input commands, these columns are intentionally
    omitted and input-specific fields should be provided via ``extra``:
    n_precursors_input_a, n_transitions_input_a, n_precursors_input_b,
    n_transitions_input_b, etc. This avoids NaN input fields in union/consensus
    reports.
    """
    counts = df_output.groupby(PRECURSOR_KEY).size() if len(df_output) else pd.Series(dtype=int)

    row: dict[str, Any] = {
        "compactlib_version": __version__,
        "mode": mode,
        "output_file": str(output_path),
    }

    if df_input is not None:
        row.update({
            "n_precursors_input": _n_precursors(df_input),
            "n_transitions_input": int(len(df_input)),
        })

    row.update({
        "n_precursors_output": _n_precursors(df_output),
        "n_transitions_output": int(len(df_output)),
        "mean_transitions_per_precursor": float(counts.mean()) if len(counts) else 0.0,
        "median_transitions_per_precursor": float(counts.median()) if len(counts) else 0.0,
        "min_transitions_per_precursor": int(counts.min()) if len(counts) else 0,
        "max_transitions_per_precursor": int(counts.max()) if len(counts) else 0,
    })

    if extra:
        row.update(extra)

    top_n = row.get("top_n")
    if top_n is not None and len(counts):
        row["n_precursors_with_less_than_n"] = int((counts < int(top_n)).sum())
    elif top_n is not None:
        row["n_precursors_with_less_than_n"] = 0

    output_path = Path(output_path).expanduser()
    if output_path.exists():
        row["library_size_mb"] = os.path.getsize(output_path) / 1024**2
    else:
        row["library_size_mb"] = None

    return pd.DataFrame([row])


def write_summary(summary: pd.DataFrame, path: str | Path) -> None:
    path = Path(path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(path, sep="\t", index=False)


def write_params(params: dict[str, Any], path: str | Path) -> None:
    path = Path(path).expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"compactlib_version": __version__, **params}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False, default=str)
