from __future__ import annotations

from abc import ABC, abstractmethod

import pandas as pd


class PredictorBackend(ABC):
    """Base interface for fragment-intensity prediction backends.

    Backends receive a precursor table and return a DIA spectral library-like
    transition table with the required compactlib columns:
    ModifiedPeptide, PrecursorCharge, PrecursorMz, ProductMz, LibraryIntensity,
    FragmentType, FragmentSeriesNumber, FragmentCharge.
    """

    name: str = "base"

    @abstractmethod
    def predict(self, precursors: pd.DataFrame) -> pd.DataFrame:
        raise NotImplementedError


def get_predictor(
    backend: str,
    fragment_types: str = "b,y",
    fragment_charges: str = "1",
    max_fragment_series: int | None = None,
    model: str | None = None,
    server_url: str = "koina.wilhelmlab.org:443",
    ssl: bool = True,
    collision_energy: float = 30.0,
    batch_size: int = 1024,
    sequence_source: str = "stripped",
    sequence_input: str = "peptide_sequences",
    charge_input: str = "precursor_charges",
    ce_input: str = "collision_energies",
    intensity_col: str = "intensities",
    annotation_col: str = "annotation",
    mz_col: str | None = None,
    drop_zero_intensity: bool = True,
) -> PredictorBackend:
    backend = backend.lower().strip()
    if backend == "mock":
        from .mock import MockPredictor

        return MockPredictor(
            fragment_types=fragment_types,
            fragment_charges=fragment_charges,
            max_fragment_series=max_fragment_series,
        )
    if backend == "koina":
        from .koina import KoinaPredictor

        if not model:
            raise ValueError("Koina backend requires --model, e.g. Prosit_2020_intensity_HCD")
        return KoinaPredictor(
            model_name=model,
            server_url=server_url,
            ssl=ssl,
            collision_energy=collision_energy,
            batch_size=batch_size,
            sequence_source=sequence_source,
            sequence_input=sequence_input,
            charge_input=charge_input,
            ce_input=ce_input,
            intensity_col=intensity_col,
            annotation_col=annotation_col,
            mz_col=mz_col,
            drop_zero_intensity=drop_zero_intensity,
        )
    raise ValueError(
        f"Unsupported predictor backend {backend!r}. "
        "Supported backends: mock, koina."
    )
