from __future__ import annotations

import numpy as np
import pandas as pd

from .io import PRECURSOR_KEY
from .selection import deduplicate_transitions, drop_helper_columns, sort_for_intensity_rank


def select_random_n(df: pd.DataFrame, top_n: int, seed: int = 42, deduplicate: bool = True) -> tuple[pd.DataFrame, dict]:
    """Randomly select up to N transitions per precursor with a fixed seed."""
    if top_n <= 0:
        raise ValueError("top_n must be positive")

    work = df.copy()
    duplicates_removed = 0
    if deduplicate:
        work, duplicates_removed = deduplicate_transitions(work)

    # Stable order before group iteration, then deterministic RNG.
    work = sort_for_intensity_rank(work)
    rng = np.random.default_rng(seed)
    chosen_parts = []

    for _, group in work.groupby(PRECURSOR_KEY, sort=False):
        idx = np.asarray(group.index)
        k = min(top_n, len(idx))
        chosen = rng.choice(idx, size=k, replace=False)
        # Deterministic output ordering independent of random order.
        chosen_parts.append(group.loc[sorted(chosen)])

    out = pd.concat(chosen_parts, axis=0) if chosen_parts else work.iloc[0:0].copy()
    stats = {
        "top_n": top_n,
        "seed": seed,
        "n_duplicate_transitions_removed": duplicates_removed,
    }
    return drop_helper_columns(out), stats


def reverse_intensity(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Reverse intensity pattern within each precursor while preserving transitions.

    The same rows/transitions are retained. Within each precursor, the smallest
    intensity receives the largest original intensity, the second-smallest receives
    the second-largest, etc. This matches the anti-profile control used in the paper.
    """
    out = df.copy()

    def _reverse(values: pd.Series) -> pd.Series:
        arr = values.to_numpy(copy=True)
        order_min = np.argsort(arr, kind="mergesort")
        order_max = order_min[::-1]
        new = np.empty_like(arr)
        new[order_min] = arr[order_max]
        return pd.Series(new, index=values.index)

    out["LibraryIntensity"] = out.groupby(PRECURSOR_KEY, group_keys=False)["LibraryIntensity"].apply(_reverse)
    stats = {"n_transitions_preserved": len(out)}
    return out, stats
