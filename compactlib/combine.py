from __future__ import annotations

import pandas as pd

from .io import PRECURSOR_KEY, TRANSITION_KEY
from .selection import (
    TIE_BREAKERS,
    add_transition_rank,
    drop_helper_columns,
    select_top_n,
    sort_for_intensity_rank,
)


def _precursor_set(df: pd.DataFrame) -> set[tuple]:
    return set(map(tuple, df[PRECURSOR_KEY].drop_duplicates().to_numpy()))


def _input_stats(df_a: pd.DataFrame, df_b: pd.DataFrame) -> dict:
    prec_a = _precursor_set(df_a)
    prec_b = _precursor_set(df_b)
    return {
        "n_precursors_input_a": len(prec_a),
        "n_transitions_input_a": int(len(df_a)),
        "n_precursors_input_b": len(prec_b),
        "n_transitions_input_b": int(len(df_b)),
        "n_common_precursors": len(prec_a & prec_b),
        "n_a_only_precursors": len(prec_a - prec_b),
        "n_b_only_precursors": len(prec_b - prec_a),
    }


def build_union(df_a: pd.DataFrame, df_b: pd.DataFrame, top_n: int) -> tuple[pd.DataFrame, dict]:
    """Build Union2N = topN(A) union topN(B), deduplicated by transition key."""
    top_a, stats_a = select_top_n(df_a, top_n=top_n, deduplicate=True)
    top_b, stats_b = select_top_n(df_b, top_n=top_n, deduplicate=True)

    a = top_a.copy()
    b = top_b.copy()
    a["_compactlib_source_order"] = 0
    b["_compactlib_source_order"] = 1
    a["_compactlib_original_order"] = range(len(a))
    b["_compactlib_original_order"] = range(len(b))

    combined = pd.concat([a, b], ignore_index=True, sort=False)

    # Prefer source A for duplicate transitions; deterministic tie-breakers afterwards.
    combined = combined.sort_values(
        TRANSITION_KEY + ["_compactlib_source_order"] + TIE_BREAKERS + ["_compactlib_original_order"],
        ascending=[True, True, True, True, True, True, True, True, True, True, True],
        kind="mergesort",
    )
    before = len(combined)
    out = combined.drop_duplicates(subset=TRANSITION_KEY, keep="first").copy()
    duplicates_removed = before - len(out)

    # Final deterministic library order.
    out = sort_for_intensity_rank(out)

    stats = {
        "top_n": top_n,
        "mode": "union",
        **_input_stats(df_a, df_b),
        "n_precursors_top_a": len(_precursor_set(top_a)),
        "n_precursors_top_b": len(_precursor_set(top_b)),
        "n_transitions_top_a": int(len(top_a)),
        "n_transitions_top_b": int(len(top_b)),
        "n_duplicate_transitions_removed_after_union": int(duplicates_removed),
        "n_duplicate_transitions_removed_a": int(stats_a.get("n_duplicate_transitions_removed", 0)),
        "n_duplicate_transitions_removed_b": int(stats_b.get("n_duplicate_transitions_removed", 0)),
        "duplicate_policy": "keep_input_a",
    }
    return drop_helper_columns(out), stats


def _ranked_table(df: pd.DataFrame, rank_col: str, source_col: str) -> pd.DataFrame:
    ranked = add_transition_rank(df)
    cols = TRANSITION_KEY + [rank_col, source_col]
    ranked = ranked.rename(columns={"_compactlib_rank": rank_col})
    ranked[source_col] = True
    return ranked[cols].copy()


def build_consensus(df_a: pd.DataFrame, df_b: pd.DataFrame, top_n: int) -> tuple[pd.DataFrame, dict]:
    """Build consensus top-N using score = 1/rank_a + 1/rank_b.

    Missing transition in one model contributes 0 to the score. Output rows are
    taken from source A when available, otherwise from source B. Helper rank/score
    columns are not written to the output library.
    """
    if top_n <= 0:
        raise ValueError("top_n must be positive")

    # Dedupe before ranking full profiles.
    from .selection import deduplicate_transitions

    da, removed_a = deduplicate_transitions(df_a)
    db, removed_b = deduplicate_transitions(df_b)

    ra = _ranked_table(da, "_compactlib_rank_a", "_compactlib_in_a")
    rb = _ranked_table(db, "_compactlib_rank_b", "_compactlib_in_b")

    ranks = ra.merge(rb, on=TRANSITION_KEY, how="outer")
    ranks["_compactlib_score"] = 0.0
    ranks.loc[ranks["_compactlib_rank_a"].notna(), "_compactlib_score"] += 1.0 / ranks.loc[
        ranks["_compactlib_rank_a"].notna(), "_compactlib_rank_a"
    ]
    ranks.loc[ranks["_compactlib_rank_b"].notna(), "_compactlib_score"] += 1.0 / ranks.loc[
        ranks["_compactlib_rank_b"].notna(), "_compactlib_rank_b"
    ]

    # Build source rows. Prefer A for transitions available in both.
    a_rows = da.copy()
    a_rows["_compactlib_prefer"] = 0
    b_rows = db.copy()
    b_rows["_compactlib_prefer"] = 1
    source_rows = pd.concat([a_rows, b_rows], ignore_index=True, sort=False)
    source_rows = source_rows.sort_values(
        TRANSITION_KEY + ["_compactlib_prefer"] + TIE_BREAKERS,
        ascending=[True, True, True, True, True, True, True, True, True, True],
        kind="mergesort",
    )
    source_rows = source_rows.drop_duplicates(subset=TRANSITION_KEY, keep="first").copy()

    scored = source_rows.merge(
        ranks[TRANSITION_KEY + ["_compactlib_rank_a", "_compactlib_rank_b", "_compactlib_score"]],
        on=TRANSITION_KEY,
        how="inner",
    )

    # Deterministic top-N by consensus score, then rank_a, rank_b, tie-breakers.
    scored["_compactlib_rank_a_sort"] = scored["_compactlib_rank_a"].fillna(10**12)
    scored["_compactlib_rank_b_sort"] = scored["_compactlib_rank_b"].fillna(10**12)
    scored = scored.sort_values(
        PRECURSOR_KEY
        + ["_compactlib_score", "_compactlib_rank_a_sort", "_compactlib_rank_b_sort"]
        + TIE_BREAKERS,
        ascending=[True, True, False, True, True, True, True, True, True],
        kind="mergesort",
    )
    out = scored.groupby(PRECURSOR_KEY, sort=False).head(top_n).copy()

    stats = {
        "top_n": top_n,
        "mode": "consensus",
        **_input_stats(df_a, df_b),
        "score": "1/rank_a + 1/rank_b",
        "n_duplicate_transitions_removed_a": int(removed_a),
        "n_duplicate_transitions_removed_b": int(removed_b),
        "n_candidate_transitions_after_merge": int(len(ranks)),
        "duplicate_policy": "prefer_input_a",
    }
    return drop_helper_columns(out), stats
