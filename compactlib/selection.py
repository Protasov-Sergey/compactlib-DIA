from __future__ import annotations

import pandas as pd

from .io import PRECURSOR_KEY, TRANSITION_KEY

TIE_BREAKERS = ["FragmentType", "FragmentSeriesNumber", "FragmentCharge", "ProductMz"]


def _with_original_order(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "_compactlib_original_order" not in out.columns:
        out["_compactlib_original_order"] = range(len(out))
    return out


def drop_helper_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.drop(columns=[c for c in df.columns if c.startswith("_compactlib_")], errors="ignore")


def sort_for_intensity_rank(df: pd.DataFrame) -> pd.DataFrame:
    """Deterministic precursor-wise intensity sorting.

    Primary key: LibraryIntensity descending.
    Tie-breakers: FragmentType, FragmentSeriesNumber, FragmentCharge, ProductMz.
    """
    out = _with_original_order(df)
    return out.sort_values(
        PRECURSOR_KEY + ["LibraryIntensity"] + TIE_BREAKERS + ["_compactlib_original_order"],
        ascending=[True, True, False, True, True, True, True, True],
        kind="mergesort",
    )


def deduplicate_transitions(df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Remove duplicate transition annotations within precursor.

    For duplicate transition keys, keep the most intense row. Ties are resolved
    deterministically by FragmentType, FragmentSeriesNumber, FragmentCharge, ProductMz.
    """
    before = len(df)
    ranked = sort_for_intensity_rank(df)
    deduped = ranked.drop_duplicates(subset=TRANSITION_KEY, keep="first").copy()
    removed = before - len(deduped)
    return deduped, removed


def add_transition_rank(df: pd.DataFrame) -> pd.DataFrame:
    """Add 1-based deterministic intensity rank within each precursor."""
    ranked = sort_for_intensity_rank(df).copy()
    ranked["_compactlib_rank"] = ranked.groupby(PRECURSOR_KEY, sort=False).cumcount() + 1
    return ranked


def select_top_n(df: pd.DataFrame, top_n: int, deduplicate: bool = True) -> tuple[pd.DataFrame, dict]:
    """Keep up to top_n most intense unique transitions per precursor."""
    if top_n <= 0:
        raise ValueError("top_n must be positive")

    work = df.copy()
    duplicates_removed = 0
    if deduplicate:
        work, duplicates_removed = deduplicate_transitions(work)

    ranked = sort_for_intensity_rank(work)
    selected = ranked.groupby(PRECURSOR_KEY, sort=False).head(top_n).copy()

    stats = {
        "top_n": top_n,
        "n_duplicate_transitions_removed": duplicates_removed,
    }
    return drop_helper_columns(selected), stats
