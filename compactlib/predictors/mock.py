from __future__ import annotations

import hashlib
from dataclasses import dataclass

import pandas as pd

from .base import PredictorBackend
from ..mass import fragment_mz, parse_positions_1based


def _parse_csv_ints(value: str) -> list[int]:
    return [int(x.strip()) for x in str(value).split(",") if x.strip()]


def _parse_csv_strs(value: str) -> list[str]:
    return [x.strip().lower() for x in str(value).split(",") if x.strip()]


def _stable_noise(*parts: object) -> float:
    s = "|".join(map(str, parts)).encode("utf-8")
    h = hashlib.md5(s).hexdigest()[:8]
    return int(h, 16) / 0xFFFFFFFF


@dataclass
class MockPredictor(PredictorBackend):
    """Deterministic toy predictor for testing compactlib's prediction interface.

    This backend is not a scientific Prosit/MS2PIP replacement. It generates
    simple b/y transitions with deterministic pseudo-intensities so that the
    end-to-end FASTA -> precursors -> transitions -> compact library pipeline
    can be tested without external services.
    """

    fragment_types: str = "b,y"
    fragment_charges: str = "1"
    max_fragment_series: int | None = None
    name: str = "mock"

    def predict(self, precursors: pd.DataFrame) -> pd.DataFrame:
        required = ["ModifiedPeptide", "StrippedPeptide", "PrecursorCharge", "PrecursorMz"]
        missing = [c for c in required if c not in precursors.columns]
        if missing:
            raise ValueError(f"Precursor table is missing required columns for mock prediction: {missing}")

        ion_types = _parse_csv_strs(self.fragment_types)
        ion_charges = _parse_csv_ints(self.fragment_charges)
        bad_types = [t for t in ion_types if t not in {"b", "y"}]
        if bad_types:
            raise ValueError(f"Mock predictor supports only b/y fragment types, got: {bad_types}")

        rows = []
        for _, row in precursors.iterrows():
            seq = str(row["StrippedPeptide"]).upper()
            if len(seq) < 2:
                continue
            ox_pos = parse_positions_1based(row.get("OxidationMPositions", ""))
            max_series = len(seq) - 1
            if self.max_fragment_series is not None:
                max_series = min(max_series, int(self.max_fragment_series))

            for ion_type in ion_types:
                for n in range(1, max_series + 1):
                    for frag_z in ion_charges:
                        try:
                            pmz = fragment_mz(
                                seq,
                                ion_type=ion_type,
                                series_number=n,
                                charge=frag_z,
                                carbamidomethyl_c=bool(row.get("FixedCarbamidomethylC", True)),
                                oxidation_m_positions=ox_pos,
                            )
                        except ValueError:
                            continue

                        # Deterministic, smooth toy intensity profile: y ions and
                        # longer fragments tend to be more intense, with a small
                        # stable pseudo-random component for tie-breaking realism.
                        base = n / max(1, len(seq) - 1)
                        ion_bonus = 1.0 if ion_type == "y" else 0.78
                        charge_penalty = 1.0 / frag_z
                        noise = 0.05 * _stable_noise(row["ModifiedPeptide"], row["PrecursorCharge"], ion_type, n, frag_z)
                        intensity = 1000.0 * (ion_bonus * base * charge_penalty + noise)

                        out = row.to_dict()
                        out.update(
                            {
                                "ProductMz": pmz,
                                "LibraryIntensity": float(intensity),
                                "FragmentType": ion_type,
                                "FragmentSeriesNumber": int(n),
                                "FragmentCharge": int(frag_z),
                                "PredictionBackend": self.name,
                            }
                        )
                        rows.append(out)

        df = pd.DataFrame(rows)
        if len(df):
            df = df.sort_values(
                [
                    "ModifiedPeptide",
                    "PrecursorCharge",
                    "FragmentType",
                    "FragmentSeriesNumber",
                    "FragmentCharge",
                    "ProductMz",
                ],
                ascending=[True, True, True, True, True, True],
                kind="mergesort",
            ).reset_index(drop=True)
        return df
