from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .base import PredictorBackend
from ..mass import fragment_mz, parse_positions_1based


_ANNOT_RE = re.compile(r"(?P<ion>[byBY])\s*(?P<num>\d+)(?:[^0-9]+(?P<charge>\d+))?")


def _is_listlike(value: Any) -> bool:
    if isinstance(value, (str, bytes)):
        return False
    if isinstance(value, (list, tuple, np.ndarray, pd.Series)):
        return True
    return hasattr(value, "__iter__")


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, float) and np.isnan(value):
        return []
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, pd.Series):
        return value.tolist()
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]




def _contains_modification_notation(values: pd.Series) -> bool:
    """Return True if peptide strings contain common modification syntax.

    Koina models such as Prosit_2020_intensity_HCD and ms2pip_HCD2021 use
    plain peptide sequence inputs (``peptide_sequences``). Passing compactlib
    modified notation such as ``C(UniMod:4)`` would reach the Koina peptide
    preprocessor as non-amino-acid characters and produce a remote model error.
    """
    text = values.astype(str)
    return bool(text.str.contains(r"\(|\[|UniMod|UNIMOD|\+\d", regex=True, na=False).any())


def _has_variable_modifications(batch: pd.DataFrame) -> bool:
    """Detect precursor tables where stripped sequence would lose PTM information."""
    if "NVariableMods" in batch.columns:
        vals = pd.to_numeric(batch["NVariableMods"], errors="coerce").fillna(0)
        if bool((vals > 0).any()):
            return True
    if "VariableOxidationM" in batch.columns:
        vals = batch["VariableOxidationM"].astype(str).str.lower()
        if bool(vals.isin(["true", "1", "yes"]).any()):
            return True
    if "OxidationMPositions" in batch.columns:
        vals = batch["OxidationMPositions"].fillna("").astype(str).str.strip()
        if bool((vals != "").any()):
            return True
    return False


def parse_fragment_annotation(annotation: Any, default_charge: int = 1) -> tuple[str, int, int]:
    """Parse common Koina/Prosit-style fragment annotations.

    Supported examples include: ``y7``, ``y7^2``, ``y7+2``, ``b3/1``.
    Returns ``(FragmentType, FragmentSeriesNumber, FragmentCharge)``.
    """
    s = str(annotation).strip()
    m = _ANNOT_RE.search(s)
    if not m:
        raise ValueError(f"Could not parse fragment annotation: {annotation!r}")
    ion = m.group("ion").lower()
    num = int(m.group("num"))
    charge = int(m.group("charge")) if m.group("charge") else int(default_charge)
    return ion, num, charge


def _find_column(df: pd.DataFrame, preferred: str | None, candidates: Iterable[str]) -> str | None:
    if preferred and preferred in df.columns:
        return preferred
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return None



def _first_non_missing(series: pd.Series) -> Any:
    """Return the first non-missing value without treating list/array values as missing."""
    for value in series:
        if value is None:
            continue
        if isinstance(value, float) and np.isnan(value):
            continue
        try:
            missing = pd.isna(value)
            if isinstance(missing, (bool, np.bool_)) and bool(missing):
                continue
        except Exception:
            pass
        return value
    return None


def _normalise_scalar(value: Any) -> str:
    """Normalise scalar values returned by Koina for key matching."""
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    return str(value)


def _to_int(value: Any, field_name: str) -> int:
    try:
        return int(float(value))
    except Exception as e:
        raise ValueError(f"Could not parse integer field {field_name}={value!r}") from e


@dataclass
class KoinaPredictor(PredictorBackend):
    """Koina/koinapy prediction backend.

    This backend delegates ML prediction to a Koina model and converts common
    array-style prediction outputs into compactlib's long transition table.

    Notes
    -----
    Koina model I/O can differ between models. Therefore column names for model
    inputs and prediction outputs are configurable from the CLI.
    """

    model_name: str
    server_url: str = "koina.wilhelmlab.org:443"
    ssl: bool = True
    collision_energy: float = 30.0
    batch_size: int = 1024
    sequence_source: str = "stripped"  # stripped or modified
    sequence_input: str = "peptide_sequences"
    charge_input: str = "precursor_charges"
    ce_input: str = "collision_energies"
    intensity_col: str = "intensities"
    annotation_col: str = "annotation"
    mz_col: str | None = None
    drop_zero_intensity: bool = True
    name: str = "koina"

    def __post_init__(self) -> None:
        try:
            from koinapy import Koina  # type: ignore
        except ImportError as e:
            raise ImportError(
                "Koina backend requires koinapy. Install with: "
                "pip install compact-dia-library[koina] or pip install koinapy"
            ) from e

        if self.batch_size <= 0:
            raise ValueError("batch_size must be positive")
        if self.sequence_source not in {"modified", "stripped"}:
            raise ValueError("sequence_source must be 'modified' or 'stripped'")

        self._client = Koina(self.model_name, self.server_url, ssl=self.ssl)

    def _make_inputs(self, batch: pd.DataFrame) -> pd.DataFrame:
        seq_col = "ModifiedPeptide" if self.sequence_source == "modified" else "StrippedPeptide"
        if seq_col not in batch.columns:
            raise ValueError(f"Precursor table is missing required sequence column: {seq_col}")
        if "PrecursorCharge" not in batch.columns:
            raise ValueError("Precursor table is missing required column: PrecursorCharge")

        sequences = batch[seq_col].astype(str)

        if self.sequence_source == "modified" and self.sequence_input == "peptide_sequences":
            if _contains_modification_notation(sequences):
                raise ValueError(
                    "Koina model input 'peptide_sequences' expects plain amino-acid sequences, "
                    "but --sequence-source modified produced strings with modification notation "
                    "such as C(UniMod:4) or M(UniMod:35). For fixed carbamidomethyl-C libraries, "
                    "rerun with --sequence-source stripped. For variable modifications, use a Koina "
                    "model/input schema that explicitly supports modified sequences."
                )

        if self.sequence_source == "stripped" and _has_variable_modifications(batch):
            raise ValueError(
                "The precursor table contains variable modifications, but --sequence-source stripped "
                "would remove PTM information before prediction. Filter variable-modified precursors "
                "or use a model/input schema that explicitly supports modified sequences."
            )

        inputs = pd.DataFrame()
        inputs[self.sequence_input] = sequences.values
        inputs[self.charge_input] = pd.to_numeric(batch["PrecursorCharge"], errors="raise").astype(int).values
        if self.ce_input:
            inputs[self.ce_input] = float(self.collision_energy)
        return inputs

    def _predict_batch(self, batch: pd.DataFrame) -> pd.DataFrame:
        inputs = self._make_inputs(batch)
        pred = self._client.predict(inputs)
        if not isinstance(pred, pd.DataFrame):
            pred = pd.DataFrame(pred)
        return pred

    def predict(self, precursors: pd.DataFrame) -> pd.DataFrame:
        required = ["ModifiedPeptide", "PrecursorCharge", "PrecursorMz"]
        missing = [c for c in required if c not in precursors.columns]
        if missing:
            raise ValueError(f"Precursor table is missing required columns for Koina prediction: {missing}")
        if self.sequence_source == "stripped" and "StrippedPeptide" not in precursors.columns:
            raise ValueError("sequence_source='stripped' requires StrippedPeptide column")
        if self.sequence_source == "modified" and "ModifiedPeptide" not in precursors.columns:
            raise ValueError("sequence_source='modified' requires ModifiedPeptide column")

        parts: list[pd.DataFrame] = []
        n = len(precursors)
        for start in range(0, n, self.batch_size):
            stop = min(start + self.batch_size, n)
            batch = precursors.iloc[start:stop].reset_index(drop=True)
            pred = self._predict_batch(batch)
            parts.append(self._convert_predictions(batch, pred))

        if parts:
            out = pd.concat(parts, ignore_index=True)
        else:
            out = pd.DataFrame()

        if len(out):
            out = out.sort_values(
                [
                    "ModifiedPeptide", "PrecursorCharge", "FragmentType",
                    "FragmentSeriesNumber", "FragmentCharge", "ProductMz",
                ],
                ascending=[True, True, True, True, True, True],
                kind="mergesort",
            ).reset_index(drop=True)
        return out

    def _convert_predictions(self, batch: pd.DataFrame, pred: pd.DataFrame) -> pd.DataFrame:
        # If the backend already returns compactlib-compatible long format,
        # accept it and only add backend metadata.
        long_required = {
            "ModifiedPeptide", "PrecursorCharge", "PrecursorMz", "ProductMz",
            "LibraryIntensity", "FragmentType", "FragmentSeriesNumber", "FragmentCharge",
        }
        if long_required.issubset(pred.columns):
            out = pred.copy()
            out["PredictionBackend"] = self.name
            out["PredictionModel"] = self.model_name
            return out

        intensity_col = _find_column(
            pred,
            self.intensity_col,
            ["intensities", "intensity", "LibraryIntensity", "fragment_intensities"],
        )
        annotation_col = _find_column(
            pred,
            self.annotation_col,
            ["annotation", "annotations", "ion", "ions", "fragment_annotations"],
        )
        mz_col = _find_column(
            pred,
            self.mz_col,
            ["mz", "mzs", "ProductMz", "product_mz", "fragment_mz"],
        )

        if intensity_col is None:
            raise ValueError(
                "Could not find intensity output column in Koina prediction. "
                f"Available columns: {list(pred.columns)}. Use --koina-intensity-col."
            )
        if annotation_col is None:
            raise ValueError(
                "Could not find annotation output column in Koina prediction. "
                f"Available columns: {list(pred.columns)}. Use --koina-annotation-col."
            )

        # Koina models can return predictions in two common layouts:
        #   1) one row per precursor with array-valued intensities/annotations;
        #   2) one row per fragment with scalar intensity/annotation and repeated
        #      input columns such as peptide_sequences and precursor_charges.
        sample_intensity = _first_non_missing(pred[intensity_col])
        sample_annotation = _first_non_missing(pred[annotation_col])
        array_style = _is_listlike(sample_intensity) or _is_listlike(sample_annotation)

        if array_style:
            return self._convert_array_predictions(batch, pred, intensity_col, annotation_col, mz_col)
        return self._convert_long_predictions(batch, pred, intensity_col, annotation_col, mz_col)

    def _convert_array_predictions(
        self,
        batch: pd.DataFrame,
        pred: pd.DataFrame,
        intensity_col: str,
        annotation_col: str,
        mz_col: str | None,
    ) -> pd.DataFrame:
        """Convert one-row-per-precursor Koina output with array-valued columns."""
        if len(pred) != len(batch):
            raise ValueError(
                "Koina prediction output looks array-valued but does not have one row per input precursor: "
                f"predictions={len(pred)}, input={len(batch)}. Available columns: {list(pred.columns)}."
            )

        rows: list[dict[str, Any]] = []
        for i in range(len(batch)):
            prow = batch.iloc[i]
            pred_row = pred.iloc[i]
            intensities = _as_list(pred_row[intensity_col])
            annotations = _as_list(pred_row[annotation_col])
            mzs = _as_list(pred_row[mz_col]) if mz_col is not None else []

            if len(annotations) != len(intensities):
                raise ValueError(
                    f"Annotation/intensity length mismatch for row {i}: "
                    f"{len(annotations)} vs {len(intensities)}"
                )
            if mzs and len(mzs) != len(intensities):
                raise ValueError(
                    f"m/z/intensity length mismatch for row {i}: {len(mzs)} vs {len(intensities)}"
                )

            rows.extend(self._rows_for_precursor(prow, intensities, annotations, mzs))

        return pd.DataFrame(rows)

    def _convert_long_predictions(
        self,
        batch: pd.DataFrame,
        pred: pd.DataFrame,
        intensity_col: str,
        annotation_col: str,
        mz_col: str | None,
    ) -> pd.DataFrame:
        """Convert one-row-per-fragment Koina output.

        Recent Koina/koinapy calls for Prosit-like models can return a long
        table: each row is a fragment, while input columns such as
        peptide_sequences and precursor_charges are repeated. This adapter maps
        each fragment row back to the corresponding compactlib precursor row.
        """
        seq_pred_col = _find_column(
            pred,
            self.sequence_input,
            [self.sequence_input, "peptide_sequences", "peptide_sequence", "sequence", "sequences"],
        )
        charge_pred_col = _find_column(
            pred,
            self.charge_input,
            [self.charge_input, "precursor_charges", "precursor_charge", "precursorCharge", "charge", "charges"],
        )
        index_pred_col = _find_column(
            pred,
            None,
            ["input_index", "input_id", "precursor_index", "row_index", "sample_index"],
        )

        batch_seq_col = "ModifiedPeptide" if self.sequence_source == "modified" else "StrippedPeptide"

        batch_by_key: dict[tuple[str, int], pd.Series] = {}
        duplicate_keys: set[tuple[str, int]] = set()
        for i in range(len(batch)):
            prow = batch.iloc[i]
            key = (_normalise_scalar(prow[batch_seq_col]), _to_int(prow["PrecursorCharge"], "PrecursorCharge"))
            if key in batch_by_key:
                duplicate_keys.add(key)
            batch_by_key[key] = prow

        if duplicate_keys and seq_pred_col is not None and charge_pred_col is not None:
            raise ValueError(
                "Cannot unambiguously map long Koina output to precursor rows because the input batch "
                f"contains duplicated ({batch_seq_col}, PrecursorCharge) keys. Example: {next(iter(duplicate_keys))!r}. "
                "Use a precursor table with unique peptide/charge rows or a Koina output that contains input_index."
            )

        rows: list[dict[str, Any]] = []
        missing_keys: set[tuple[str, int]] = set()

        for pred_i in range(len(pred)):
            pred_row = pred.iloc[pred_i]

            if seq_pred_col is not None and charge_pred_col is not None:
                key = (
                    _normalise_scalar(pred_row[seq_pred_col]),
                    _to_int(pred_row[charge_pred_col], charge_pred_col),
                )
                prow = batch_by_key.get(key)
                if prow is None:
                    missing_keys.add(key)
                    continue
            elif index_pred_col is not None:
                idx = _to_int(pred_row[index_pred_col], index_pred_col)
                # Support both zero-based and one-based indices if possible.
                if 0 <= idx < len(batch):
                    prow = batch.iloc[idx]
                elif 1 <= idx <= len(batch):
                    prow = batch.iloc[idx - 1]
                else:
                    raise ValueError(
                        f"Koina output index column {index_pred_col!r} contains out-of-range value {idx}; "
                        f"batch size is {len(batch)}."
                    )
            else:
                raise ValueError(
                    "Koina prediction output appears to be one-row-per-fragment, but compactlib cannot map "
                    "fragment rows back to input precursors. The output lacks peptide/charge input columns "
                    f"({self.sequence_input!r}, {self.charge_input!r}) and lacks an input index column. "
                    f"Available columns: {list(pred.columns)}."
                )

            inten = pred_row[intensity_col]
            ann = pred_row[annotation_col]
            mzs = [pred_row[mz_col]] if mz_col is not None else []
            rows.extend(self._rows_for_precursor(prow, [inten], [ann], mzs))

        if missing_keys:
            example = next(iter(missing_keys))
            raise ValueError(
                "Could not map some Koina long-format prediction rows back to input precursors. "
                f"Missing keys: {len(missing_keys)}; example={example!r}; "
                f"batch sequence column={batch_seq_col!r}; Koina columns={seq_pred_col!r}, {charge_pred_col!r}."
            )

        return pd.DataFrame(rows)

    def _rows_for_precursor(
        self,
        prow: pd.Series,
        intensities: list[Any],
        annotations: list[Any],
        mzs: list[Any],
    ) -> list[dict[str, Any]]:
        """Create compactlib transition rows for one precursor."""
        if len(annotations) != len(intensities):
            raise ValueError(
                "Annotation/intensity length mismatch: "
                f"{len(annotations)} vs {len(intensities)}"
            )
        if mzs and len(mzs) != len(intensities):
            raise ValueError(
                f"m/z/intensity length mismatch: {len(mzs)} vs {len(intensities)}"
            )

        stripped = str(prow.get("StrippedPeptide", "")).upper()
        ox_pos = parse_positions_1based(prow.get("OxidationMPositions", ""))
        carb_c = bool(prow.get("FixedCarbamidomethylC", True))

        rows: list[dict[str, Any]] = []
        for j, inten in enumerate(intensities):
            try:
                inten = float(inten)
            except Exception:
                continue
            if not np.isfinite(inten):
                continue
            if self.drop_zero_intensity and inten <= 0:
                continue

            # Some long-format outputs can provide a fragment charge column via the
            # annotation itself. If the annotation does not specify charge, default
            # to 1, which is correct for Prosit_2020_intensity_HCD y/b outputs.
            ion_type, series_number, frag_charge = parse_fragment_annotation(annotations[j], default_charge=1)
            if mzs:
                product_mz = float(mzs[j])
            else:
                if not stripped:
                    raise ValueError("Cannot calculate ProductMz without StrippedPeptide column or Koina m/z output")
                try:
                    product_mz = fragment_mz(
                        stripped,
                        ion_type=ion_type,
                        series_number=series_number,
                        charge=frag_charge,
                        carbamidomethyl_c=carb_c,
                        oxidation_m_positions=ox_pos,
                    )
                except Exception as e:
                    raise ValueError(
                        f"Could not calculate ProductMz for annotation={annotations[j]!r}, peptide={stripped!r}. "
                        "Provide Koina m/z output via --koina-mz-col."
                    ) from e

            out = prow.to_dict()
            out.update(
                {
                    "ProductMz": product_mz,
                    "LibraryIntensity": inten,
                    "FragmentType": ion_type,
                    "FragmentSeriesNumber": int(series_number),
                    "FragmentCharge": int(frag_charge),
                    "FragmentAnnotation": str(annotations[j]),
                    "PredictionBackend": self.name,
                    "PredictionModel": self.model_name,
                    "CollisionEnergy": float(self.collision_energy),
                }
            )
            rows.append(out)

        return rows
