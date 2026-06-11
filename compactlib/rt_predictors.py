from __future__ import annotations

import inspect
import re
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

from .io import coerce_numeric_commas


_MOD_NAME_BY_TOKEN = {
    "unimod:4": "Carbamidomethyl",
    "carbamidomethyl": "Carbamidomethyl",
    "+57": "Carbamidomethyl",
    "57.021": "Carbamidomethyl",
    "57,021": "Carbamidomethyl",
    "unimod:35": "Oxidation",
    "oxidation": "Oxidation",
    "+15.99": "Oxidation",
    "+15,99": "Oxidation",
    "15.994": "Oxidation",
    "15,994": "Oxidation",
}


@dataclass
class DeepLCPredictionResult:
    precursor_table: pd.DataFrame
    deeplc_input: pd.DataFrame
    stats: dict


def _mod_name_from_token(token: str) -> Optional[str]:
    low = str(token).strip().lower()
    for needle, name in _MOD_NAME_BY_TOKEN.items():
        if needle in low:
            return name
    return None


def parse_modified_peptide_for_deeplc(
    modified_peptide: str,
    *,
    stripped_peptide: Optional[str] = None,
    fixed_carbamidomethyl_c: bool = True,
) -> tuple[str, str]:
    """Convert compactlib/Skyline-like modified peptide notation to DeepLC input.

    DeepLC expects two columns:
      - seq: unmodified peptide sequence
      - modifications: MS2PIP-style pipe-separated pairs, e.g.
        ``3|Carbamidomethyl|7|Oxidation``.

    This helper recognises the compactlib defaults ``C(UniMod:4)`` and
    ``M(UniMod:35)``, bracket mass variants such as ``C[+57.021464]`` and
    simple named modifications. Fixed carbamidomethylation can additionally be
    added to every cysteine because DeepLC documentation explicitly expects
    fixed modifications to be indicated in the input.
    """
    s = "" if pd.isna(modified_peptide) else str(modified_peptide).strip()
    seq_chars: list[str] = []
    mods: list[tuple[int, str]] = []

    i = 0
    current_pos = 0
    while i < len(s):
        ch = s[i]
        if ch.isalpha() and ch.upper() == ch and len(ch) == 1:
            current_pos += 1
            seq_chars.append(ch)
            i += 1

            # Parse one or more immediate modification tokens after the AA.
            while i < len(s) and s[i] in "([":
                open_ch = s[i]
                close_ch = ")" if open_ch == "(" else "]"
                j = s.find(close_ch, i + 1)
                if j == -1:
                    break
                token = s[i + 1 : j]
                mod_name = _mod_name_from_token(token)
                if mod_name is not None:
                    mods.append((current_pos, mod_name))
                i = j + 1
            continue
        i += 1

    seq = "".join(seq_chars)
    if not seq and stripped_peptide is not None and not pd.isna(stripped_peptide):
        seq = str(stripped_peptide).strip()
    elif stripped_peptide is not None and not pd.isna(stripped_peptide):
        # Prefer explicit stripped sequence if provided; this avoids edge cases
        # in unusual modification notations while preserving parsed mod positions.
        explicit = str(stripped_peptide).strip()
        if explicit:
            seq = explicit

    existing = set(mods)
    if fixed_carbamidomethyl_c:
        for pos, aa in enumerate(seq, start=1):
            if aa == "C" and (pos, "Carbamidomethyl") not in existing:
                mods.append((pos, "Carbamidomethyl"))
                existing.add((pos, "Carbamidomethyl"))

    mods = sorted(set(mods), key=lambda x: (x[0], x[1]))
    mod_string = "|".join([item for pos, name in mods for item in (str(pos), name)])
    return seq, mod_string


def prepare_deeplc_input(
    table: pd.DataFrame,
    *,
    modified_column: str = "ModifiedPeptide",
    stripped_column: str = "StrippedPeptide",
    fixed_carbamidomethyl_c: bool = True,
    rt_column: Optional[str] = None,
    rt_output_column: str = "tr",
    keep_key_columns: bool = True,
) -> pd.DataFrame:
    """Create a DeepLC-compatible input table from a compactlib precursor table."""
    if modified_column not in table.columns and stripped_column not in table.columns:
        raise ValueError(
            f"Need at least one peptide column for DeepLC input: '{modified_column}' or '{stripped_column}'. "
            f"Available columns: {list(table.columns)}"
        )

    rows = []
    for _, row in table.iterrows():
        modified = row[modified_column] if modified_column in table.columns else row[stripped_column]
        stripped = row[stripped_column] if stripped_column in table.columns else None
        seq, mods = parse_modified_peptide_for_deeplc(
            modified,
            stripped_peptide=stripped,
            fixed_carbamidomethyl_c=fixed_carbamidomethyl_c,
        )
        rec = {
            "seq": seq,
            "modifications": mods,
        }
        if keep_key_columns:
            if modified_column in table.columns:
                rec[modified_column] = row[modified_column]
            if stripped_column in table.columns:
                rec[stripped_column] = row[stripped_column]
        if rt_column is not None:
            if rt_column not in table.columns:
                raise ValueError(
                    f"RT column '{rt_column}' is missing from table. Available columns: {list(table.columns)}"
                )
            rec[rt_output_column] = row[rt_column]
        rows.append(rec)

    out = pd.DataFrame(rows)
    if rt_column is not None:
        out = coerce_numeric_commas(out, [rt_output_column])
    return out


def _normalise_deeplc_predictions(preds, n_expected: int) -> np.ndarray:
    if isinstance(preds, pd.DataFrame):
        for col in ["predicted_tr", "rt", "tr", "prediction", "predictions"]:
            if col in preds.columns:
                arr = preds[col].to_numpy()
                break
        else:
            if preds.shape[1] != 1:
                raise ValueError(
                    "DeepLC returned a DataFrame with multiple columns and no recognised RT column: "
                    f"{list(preds.columns)}"
                )
            arr = preds.iloc[:, 0].to_numpy()
    elif isinstance(preds, pd.Series):
        arr = preds.to_numpy()
    else:
        arr = np.asarray(preds)

    arr = np.asarray(arr, dtype=float).reshape(-1)
    if len(arr) != n_expected:
        raise ValueError(f"DeepLC returned {len(arr)} predictions for {n_expected} input peptides")
    return arr


def _instantiate_deeplc(model_path: Optional[str] = None):
    try:
        from deeplc import DeepLC  # type: ignore
    except ImportError as exc:
        raise ImportError(
            "DeepLC backend requested but the 'deeplc' package is not installed. "
            "Install it with: pip install compact-dia-library[deeplc]. "
            "The current DeepLC release may require Python >=3.11; for older Python versions, "
            "install a compatible DeepLC release manually."
        ) from exc

    if model_path:
        # DeepLC versions have used slightly different constructor names. Try the
        # most common ones, then fall back to the default constructor with a clear
        # warning in the stats if needed.
        attempts = [
            {"path_model": model_path},
            {"path_models": model_path},
            {"model_path": model_path},
        ]
        last_exc = None
        for kwargs in attempts:
            try:
                return DeepLC(**kwargs)
            except TypeError as exc:
                last_exc = exc
        raise TypeError(
            "Could not pass --deeplc-model-path to DeepLC constructor. "
            "Your DeepLC version may use a different API. Try omitting --deeplc-model-path."
        ) from last_exc
    return DeepLC()


def _call_deeplc_make_preds(dlc, deeplc_input: pd.DataFrame, *, calibrate: Optional[bool] = None):
    """Call DeepLC.make_preds across DeepLC API versions.

    Some DeepLC versions require ``calibrate=False`` when no calibration data
    have been supplied; otherwise they raise ``AssertionError: DeepLC instance
    is not yet calibrated``. Older versions may not expose the ``calibrate``
    keyword, so we try the modern signature first and fall back gracefully.
    """
    attempts = []
    if calibrate is not None:
        attempts.extend([
            lambda: dlc.make_preds(seq_df=deeplc_input, calibrate=calibrate),
            lambda: dlc.make_preds(deeplc_input, calibrate=calibrate),
        ])
    attempts.extend([
        lambda: dlc.make_preds(seq_df=deeplc_input),
        lambda: dlc.make_preds(deeplc_input),
    ])

    last_type_error = None
    for attempt in attempts:
        try:
            return attempt()
        except TypeError as exc:
            last_type_error = exc
            continue

    if last_type_error is not None:
        raise last_type_error
    raise RuntimeError("DeepLC.make_preds failed without returning predictions")


def predict_rt_with_deeplc(
    deeplc_input: pd.DataFrame,
    *,
    calibration_input: Optional[pd.DataFrame] = None,
    model_path: Optional[str] = None,
) -> tuple[np.ndarray, dict]:
    """Run DeepLC RT prediction through its Python API."""
    dlc = _instantiate_deeplc(model_path)

    has_calibration = calibration_input is not None and len(calibration_input) > 0
    if has_calibration:
        try:
            dlc.calibrate_preds(seq_df=calibration_input)
        except TypeError:
            dlc.calibrate_preds(calibration_input)

    # Without a calibration table, recent DeepLC versions require
    # make_preds(..., calibrate=False). Otherwise they assert that the instance
    # is not calibrated. When calibration was supplied, keep DeepLC defaults.
    preds = _call_deeplc_make_preds(
        dlc,
        deeplc_input,
        calibrate=False if not has_calibration else None,
    )

    arr = _normalise_deeplc_predictions(preds, len(deeplc_input))
    try:
        import deeplc as deeplc_module  # type: ignore
        version = getattr(deeplc_module, "__version__", "unknown")
    except Exception:
        version = "unknown"

    return arr, {
        "deeplc_version": version,
        "deeplc_model_path": model_path or "default",
        "deeplc_calibrated": has_calibration,
        "deeplc_make_preds_calibrate_arg": "False" if not has_calibration else "default_after_calibration",
    }


def predict_rt_mock(deeplc_input: pd.DataFrame) -> tuple[np.ndarray, dict]:
    """Deterministic lightweight RT predictor used for tests and examples."""
    seq_len = deeplc_input["seq"].astype(str).str.len().to_numpy(dtype=float)
    n_mods = deeplc_input["modifications"].fillna("").astype(str).map(
        lambda x: 0 if x == "" else max(1, len(x.split("|")) // 2)
    ).to_numpy(dtype=float)
    preds = seq_len * 0.5 + n_mods * 0.25
    return preds, {
        "deeplc_version": "mock",
        "deeplc_model_path": "mock",
        "deeplc_calibrated": False,
    }


def add_deeplc_rt_to_precursors(
    precursors: pd.DataFrame,
    *,
    backend: str = "deeplc",
    modified_column: str = "ModifiedPeptide",
    stripped_column: str = "StrippedPeptide",
    rt_output_column: str = "NormalizedRetentionTime",
    fixed_carbamidomethyl_c: bool = True,
    calibration_table: Optional[pd.DataFrame] = None,
    calibration_rt_column: str = "tr",
    deeplc_model_path: Optional[str] = None,
) -> DeepLCPredictionResult:
    """Predict RT for unique peptidoforms and merge it back to precursor rows."""
    if modified_column not in precursors.columns and stripped_column not in precursors.columns:
        raise ValueError(
            f"Need '{modified_column}' or '{stripped_column}' in precursor table. "
            f"Available columns: {list(precursors.columns)}"
        )

    key_cols = [c for c in [modified_column, stripped_column] if c in precursors.columns]
    unique_peptides = precursors[key_cols].drop_duplicates().reset_index(drop=True)
    deeplc_input = prepare_deeplc_input(
        unique_peptides,
        modified_column=modified_column,
        stripped_column=stripped_column,
        fixed_carbamidomethyl_c=fixed_carbamidomethyl_c,
        keep_key_columns=True,
    )

    calibration_input = None
    if calibration_table is not None:
        if {"seq", "modifications"}.issubset(calibration_table.columns):
            calibration_input = calibration_table.copy()
            if calibration_rt_column != "tr":
                if calibration_rt_column not in calibration_input.columns:
                    raise ValueError(
                        f"Calibration RT column '{calibration_rt_column}' is missing. "
                        f"Available columns: {list(calibration_input.columns)}"
                    )
                calibration_input["tr"] = calibration_input[calibration_rt_column]
            elif "tr" not in calibration_input.columns:
                raise ValueError("Calibration table with DeepLC-style seq/modifications must contain 'tr' or --calibration-rt-column")
            calibration_input = coerce_numeric_commas(calibration_input, ["tr"])
        else:
            calibration_input = prepare_deeplc_input(
                calibration_table,
                modified_column=modified_column,
                stripped_column=stripped_column,
                fixed_carbamidomethyl_c=fixed_carbamidomethyl_c,
                rt_column=calibration_rt_column,
                rt_output_column="tr",
                keep_key_columns=False,
            )

    backend_norm = backend.lower().strip()
    if backend_norm == "deeplc":
        preds, backend_stats = predict_rt_with_deeplc(
            deeplc_input[["seq", "modifications"]].copy(),
            calibration_input=calibration_input[["seq", "modifications", "tr"]].copy() if calibration_input is not None else None,
            model_path=deeplc_model_path,
        )
    elif backend_norm == "mock":
        preds, backend_stats = predict_rt_mock(deeplc_input)
    else:
        raise ValueError("Unknown RT backend: {backend}. Supported: deeplc, mock")

    rt_lookup = deeplc_input[key_cols].copy()
    rt_lookup[rt_output_column] = preds

    out = precursors.merge(rt_lookup, on=key_cols, how="left")
    n_missing = int(out[rt_output_column].isna().sum())
    if n_missing:
        raise ValueError(f"Internal RT merge failed: {n_missing} precursor rows did not receive RT")

    stats = {
        "rt_backend": backend_norm,
        "rt_output_column": rt_output_column,
        "n_precursor_rows_input": int(len(precursors)),
        "n_unique_peptidoforms_predicted": int(len(deeplc_input)),
        "n_precursor_rows_with_rt": int(out[rt_output_column].notna().sum()),
        "n_calibration_rows": int(len(calibration_input)) if calibration_input is not None else 0,
        "fixed_carbamidomethyl_c_for_deeplc": fixed_carbamidomethyl_c,
        "rt_min": float(np.nanmin(preds)) if len(preds) else np.nan,
        "rt_max": float(np.nanmax(preds)) if len(preds) else np.nan,
        "rt_median": float(np.nanmedian(preds)) if len(preds) else np.nan,
        **backend_stats,
    }
    return DeepLCPredictionResult(precursor_table=out, deeplc_input=deeplc_input, stats=stats)
