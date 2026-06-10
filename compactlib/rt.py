from __future__ import annotations

from typing import Iterable

import pandas as pd

from .io import PRECURSOR_KEY, coerce_numeric_commas


def parse_key_columns(key: str | Iterable[str]) -> list[str]:
    if isinstance(key, str):
        cols = [c.strip() for c in key.split(",") if c.strip()]
    else:
        cols = [str(c).strip() for c in key if str(c).strip()]
    if not cols:
        raise ValueError("At least one key column must be provided")
    return cols


def _normalise_key_column(series: pd.Series, column: str) -> pd.Series:
    """Return a stable string representation for join keys."""
    if column in {"PrecursorCharge", "FragmentCharge", "FragmentSeriesNumber"}:
        num = coerce_numeric_commas(pd.DataFrame({column: series}), [column])[column]
        return num.astype("Int64").astype(str)
    return series.astype(str).fillna("").str.strip()


def add_join_key(df: pd.DataFrame, key_cols: list[str], key_name: str = "__compactlib_join_key") -> pd.DataFrame:
    missing = [c for c in key_cols if c not in df.columns]
    if missing:
        raise ValueError(
            "Missing RT merge key columns: " + ", ".join(missing) +
            f". Available columns: {list(df.columns)}"
        )
    out = df.copy()
    parts = [_normalise_key_column(out[c], c) for c in key_cols]
    key = parts[0]
    for part in parts[1:]:
        key = key + "|" + part
    out[key_name] = key
    return out


def attach_rt(
    library: pd.DataFrame,
    rt_table: pd.DataFrame,
    *,
    key_cols: list[str] | None = None,
    rt_column: str,
    rt_output_column: str = "NormalizedRetentionTime",
    min_match_rate: float = 0.0,
) -> tuple[pd.DataFrame, dict]:
    """Attach an RT column from a precursor-level RT table to a transition library.

    The merge is performed on precursor keys, typically ModifiedPeptide +
    PrecursorCharge. The RT table is de-duplicated by key; if multiple rows for a
    key are present, the first non-null RT value is used and duplicate statistics
    are reported in the summary.
    """
    key_cols = list(key_cols or PRECURSOR_KEY)
    if rt_column not in rt_table.columns:
        raise ValueError(
            f"RT column '{rt_column}' is missing from RT table. "
            f"Available columns: {list(rt_table.columns)}"
        )

    lib = add_join_key(library, key_cols)
    rt = add_join_key(rt_table, key_cols)
    rt = rt.copy()
    rt[rt_output_column] = coerce_numeric_commas(pd.DataFrame({rt_output_column: rt[rt_column]}), [rt_output_column])[rt_output_column]

    n_rt_duplicate_key_rows = int(rt.duplicated("__compactlib_join_key").sum())
    n_rt_keys_with_conflicting_values = 0
    if len(rt):
        nunique_rt = rt.groupby("__compactlib_join_key")[rt_output_column].nunique(dropna=True)
        n_rt_keys_with_conflicting_values = int((nunique_rt > 1).sum())

    rt_small = (
        rt[["__compactlib_join_key", rt_output_column]]
        .dropna(subset=[rt_output_column])
        .drop_duplicates("__compactlib_join_key", keep="first")
    )

    merged = lib.merge(rt_small, on="__compactlib_join_key", how="left")
    matched_rows_mask = merged[rt_output_column].notna()

    precursor_rt = merged[["__compactlib_join_key", rt_output_column]].drop_duplicates("__compactlib_join_key")
    n_library_precursors = int(lib["__compactlib_join_key"].nunique())
    n_matched_precursors = int(precursor_rt[rt_output_column].notna().sum())
    n_unmatched_precursors = int(n_library_precursors - n_matched_precursors)
    pct_matched_precursors = float(n_matched_precursors / n_library_precursors * 100.0) if n_library_precursors else 0.0

    stats = {
        "rt_key": ",".join(key_cols),
        "rt_column": rt_column,
        "rt_output_column": rt_output_column,
        "n_library_rows": int(len(library)),
        "n_library_precursors": n_library_precursors,
        "n_rt_rows": int(len(rt_table)),
        "n_rt_unique_keys": int(rt_small["__compactlib_join_key"].nunique()),
        "n_rt_duplicate_key_rows": n_rt_duplicate_key_rows,
        "n_rt_keys_with_conflicting_values": n_rt_keys_with_conflicting_values,
        "n_matched_library_rows": int(matched_rows_mask.sum()),
        "n_unmatched_library_rows": int((~matched_rows_mask).sum()),
        "n_matched_precursors": n_matched_precursors,
        "n_unmatched_precursors": n_unmatched_precursors,
        "pct_matched_precursors": pct_matched_precursors,
        "min_match_rate": float(min_match_rate),
    }

    if n_library_precursors and (n_matched_precursors / n_library_precursors) < min_match_rate:
        raise ValueError(
            f"RT match rate is below threshold: {n_matched_precursors}/{n_library_precursors} "
            f"precursors matched ({pct_matched_precursors:.3f}%), threshold={min_match_rate * 100:.3f}%. "
            "Lower --min-match-rate if this is expected."
        )

    merged = merged.drop(columns=["__compactlib_join_key"])
    return merged, stats
