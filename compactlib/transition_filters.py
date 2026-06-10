from __future__ import annotations

from typing import Iterable

import pandas as pd


def _parse_csv_strs(value: str | None) -> list[str]:
    if value is None or str(value).strip() == "":
        return []
    return [x.strip().lower() for x in str(value).split(",") if x.strip()]


def _parse_csv_ints(value: str | None) -> list[int]:
    if value is None or str(value).strip() == "":
        return []
    return [int(x.strip()) for x in str(value).split(",") if x.strip()]


def filter_transitions(
    df: pd.DataFrame,
    fragment_types: str | None = None,
    fragment_charges: str | None = None,
    min_fragment_series: int | None = None,
    max_fragment_series: int | None = None,
    min_product_mz: float | None = None,
    max_product_mz: float | None = None,
) -> tuple[pd.DataFrame, dict]:
    """Filter predicted transitions using Skyline/DIA-NN-like practical rules.

    This is intentionally transparent and conservative. It filters by fragment
    ion type, fragment charge, fragment series number and optional product m/z.
    Missing columns are ignored except for columns required by active filters.
    """
    out = df.copy()
    before = len(out)
    stats = {
        "n_transitions_before_transition_filters": int(before),
        "transition_filter_fragment_types": fragment_types,
        "transition_filter_fragment_charges": fragment_charges,
        "transition_filter_min_fragment_series": min_fragment_series,
        "transition_filter_max_fragment_series": max_fragment_series,
        "transition_filter_min_product_mz": min_product_mz,
        "transition_filter_max_product_mz": max_product_mz,
    }

    if fragment_types:
        allowed = set(_parse_csv_strs(fragment_types))
        if allowed and "FragmentType" in out.columns:
            out = out[out["FragmentType"].astype(str).str.lower().isin(allowed)].copy()

    if fragment_charges:
        allowed_z = set(_parse_csv_ints(fragment_charges))
        if allowed_z and "FragmentCharge" in out.columns:
            charges = pd.to_numeric(out["FragmentCharge"], errors="coerce")
            out = out[charges.isin(allowed_z)].copy()

    if min_fragment_series is not None and "FragmentSeriesNumber" in out.columns:
        series = pd.to_numeric(out["FragmentSeriesNumber"], errors="coerce")
        out = out[series >= int(min_fragment_series)].copy()

    if max_fragment_series is not None and "FragmentSeriesNumber" in out.columns:
        series = pd.to_numeric(out["FragmentSeriesNumber"], errors="coerce")
        out = out[series <= int(max_fragment_series)].copy()

    if min_product_mz is not None and "ProductMz" in out.columns:
        mz = pd.to_numeric(out["ProductMz"], errors="coerce")
        out = out[mz >= float(min_product_mz)].copy()

    if max_product_mz is not None and "ProductMz" in out.columns:
        mz = pd.to_numeric(out["ProductMz"], errors="coerce")
        out = out[mz <= float(max_product_mz)].copy()

    stats["n_transitions_filtered_by_transition_filters"] = int(before - len(out))
    stats["n_transitions_after_transition_filters"] = int(len(out))
    return out, stats
