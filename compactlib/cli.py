from __future__ import annotations

from pathlib import Path
from typing import Optional
import os
import json


import typer
import pandas as pd

from .combine import build_consensus, build_union
from .controls import reverse_intensity, select_random_n
from .digest import build_precursor_table_from_fasta
from .io import read_library, read_table, sidecar_paths, write_library, write_table
from .progress import ProgressLogger, format_shape
from .presets import preset_table, resolve_digest_preset, resolve_predict_preset
from .transition_filters import filter_transitions
from .report import library_summary, write_params, write_summary
from .selection import select_top_n
from .predictors import get_predictor
from .rt import attach_rt, parse_key_columns
from .rt_predictors import add_deeplc_rt_to_precursors, prepare_deeplc_input
from . import __version__

app = typer.Typer(
    help="Build compact DIA spectral libraries and FASTA-derived precursor tables.",
    no_args_is_help=True,
)


def _finish(
    df_input,
    df_output,
    output: Path,
    mode: str,
    params: dict,
    extra_summary: Optional[dict] = None,
    logger: Optional[ProgressLogger] = None,
) -> None:
    logger = logger or ProgressLogger(verbose=True, progress=False)
    elapsed = logger.total_elapsed()

    extra = dict(extra_summary or {})
    extra["compactlib_elapsed_sec"] = elapsed

    with logger.step(f"writing output library to {output}"):
        write_library(df_output, output)

    summary_path, params_path = sidecar_paths(output)

    with logger.step(f"writing summary to {summary_path}"):
        summary = library_summary(
            df_input=df_input,
            df_output=df_output,
            output_path=output,
            mode=mode,
            extra=extra,
        )
        write_summary(summary, summary_path)

    with logger.step(f"writing parameters to {params_path}"):
        write_params(params, params_path)

    logger.log(f"[compactlib] output rows: {len(df_output):,}")
    logger.log(f"[compactlib] total elapsed: {logger.total_elapsed():.2f} s")
    logger.log(f"[compactlib] wrote library: {output}")
    logger.log(f"[compactlib] wrote summary: {summary_path}")
    logger.log(f"[compactlib] wrote params: {params_path}")


def _read_one(path: Path, logger: ProgressLogger, label: str = "input"):
    with logger.step(f"reading {label} library from {path}"):
        df = read_library(path)
    logger.log(f"[compactlib] {label} library: {format_shape(df)}")
    return df


def _common_options(
    verbose: bool = typer.Option(True, "--verbose/--quiet", help="Print progress logs."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a terminal progress spinner when possible."),
):
    return verbose, progress


def _unsupported_aa_mask_for_prediction(df: pd.DataFrame) -> pd.Series:
    """Return mask for precursor rows containing amino acids unsupported by real predictors.

    Currently this checks selenocysteine (U). Mock prediction can handle it for
    testing because compactlib knows the U mass, but Koina/Prosit/MS2PIP models
    generally do not accept U-containing peptide sequences.
    """
    if "StrippedPeptide" in df.columns:
        seq = df["StrippedPeptide"].astype(str).str.upper()
    elif "ModifiedPeptide" in df.columns:
        seq = df["ModifiedPeptide"].astype(str).str.upper()
    else:
        return pd.Series(False, index=df.index)
    return seq.str.contains("U", regex=False)


MINIMAL_LIBRARY_COLUMNS = [
    "ModifiedPeptide",
    "StrippedPeptide",
    "PrecursorCharge",
    "PrecursorMz",
    "ProductMz",
    "LibraryIntensity",
    "FragmentType",
    "FragmentSeriesNumber",
    "FragmentCharge",
    "ProteinId",
]


def _minimal_library_columns(*, include_rt: bool = False, rt_output_column: str = "NormalizedRetentionTime") -> list[str]:
    cols = list(MINIMAL_LIBRARY_COLUMNS)
    if include_rt and rt_output_column not in cols:
        cols.append(rt_output_column)
    return cols


def _ensure_rt_output_column(
    df: pd.DataFrame,
    *,
    include_rt: bool,
    rt_column: str,
    rt_output_column: str,
) -> pd.DataFrame:
    if not include_rt:
        return df
    if rt_column not in df.columns and rt_output_column not in df.columns:
        raise typer.BadParameter(
            f"--include-rt was requested, but neither RT source column '{rt_column}' "
            f"nor output column '{rt_output_column}' is present. Available columns: {list(df.columns)}"
        )
    out = df.copy()
    source_col = rt_column if rt_column in out.columns else rt_output_column
    out[rt_output_column] = out[source_col]
    return out


def _format_predicted_library_output(
    df: pd.DataFrame,
    *,
    include_metadata: bool,
    include_rt: bool = False,
    rt_column: str = "NormalizedRetentionTime",
    rt_output_column: str = "NormalizedRetentionTime",
) -> pd.DataFrame:
    """Return the user-facing predicted library table.

    By default compactlib writes a minimal transition library suitable for DIA tools.
    Digest/protein/prediction metadata can be retained with --include-metadata.
    If --include-rt is set, an RT column is passed through from the precursor
    table into the predicted transition library.
    """
    df = _ensure_rt_output_column(
        df,
        include_rt=include_rt,
        rt_column=rt_column,
        rt_output_column=rt_output_column,
    )

    if include_metadata:
        return df

    output_columns = _minimal_library_columns(include_rt=include_rt, rt_output_column=rt_output_column)
    missing = [col for col in output_columns if col not in df.columns]
    if missing:
        raise typer.BadParameter(
            "Cannot write minimal predicted library because required columns are missing: "
            + ", ".join(missing)
            + ". Re-run with --include-metadata to preserve all available columns for debugging."
        )
    return df.loc[:, output_columns].copy()



@app.command("presets")
def cli_presets() -> None:
    """Print available Skyline/DIA-NN-like practical filtering presets."""
    typer.echo(preset_table())

@app.command("digest")
def cli_digest(
    fasta: Path = typer.Option(..., "--fasta", "-f", help="Input FASTA file"),
    output: Path = typer.Option(..., "--output", "-o", help="Output precursor table: .tsv, .csv, .parquet"),
    preset: str = typer.Option("none", "--preset", help="Filtering preset: none, generic-dia, generic-dia-windowed, generic-dia-strict, diann-like, skyline-dia-like, skyline-export-like"),
    enzyme: Optional[str] = typer.Option(None, "--enzyme", help="Override enzyme: trypsin-p or trypsin"),
    missed_cleavages: Optional[int] = typer.Option(None, "--missed-cleavages", min=0, help="Override maximum missed cleavages"),
    min_length: Optional[int] = typer.Option(None, "--min-length", min=1, help="Override minimum peptide length"),
    max_length: Optional[int] = typer.Option(None, "--max-length", min=1, help="Override maximum peptide length"),
    charges: Optional[str] = typer.Option(None, "--charges", help="Override comma-separated precursor charges, e.g. 2,3,4"),
    min_precursor_mz: Optional[float] = typer.Option(None, "--min-precursor-mz", help="Optional minimum precursor m/z filter, e.g. 400."),
    max_precursor_mz: Optional[float] = typer.Option(None, "--max-precursor-mz", help="Optional maximum precursor m/z filter, e.g. 1200."),
    charge_length_rules: Optional[str] = typer.Option(None, "--charge-length-rules", help="Override minimum peptide length per charge, e.g. '4:12,5:16'."),
    carbamidomethyl_c: bool = typer.Option(True, "--carbamidomethyl-c/--no-carbamidomethyl-c", help="Apply fixed C carbamidomethylation to mass and ModifiedPeptide"),
    c_mod_format: str = typer.Option("unimod", "--c-mod-format", help="C modification format: unimod, bracket-unimod, name, plain"),
    remove_invalid_aa: bool = typer.Option(True, "--remove-invalid-aa/--keep-invalid-aa", help="Remove peptides with non-canonical amino acids"),
    allow_selenocysteine: bool = typer.Option(False, "--allow-selenocysteine/--no-allow-selenocysteine", help="Retain U-containing peptides during digest. Most prediction backends do not support U."),
    variable_oxidation_m: bool = typer.Option(False, "--variable-oxidation-m/--no-variable-oxidation-m", help="Expand variable methionine oxidation peptidoforms, M(UniMod:35)."),
    max_variable_mods: int = typer.Option(1, "--max-variable-mods", min=0, help="Maximum number of variable M oxidations per peptide."),
    m_mod_format: str = typer.Option("unimod", "--m-mod-format", help="M oxidation format: unimod, bracket-unimod, name, plain"),
    verbose: bool = typer.Option(True, "--verbose/--quiet", help="Print progress logs."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a terminal progress spinner when possible."),
) -> None:
    """Digest FASTA and create a precursor table for later prediction."""
    digest_cfg = resolve_digest_preset(
        preset=preset,
        enzyme=enzyme,
        missed_cleavages=missed_cleavages,
        min_length=min_length,
        max_length=max_length,
        charges=charges,
        min_precursor_mz=min_precursor_mz,
        max_precursor_mz=max_precursor_mz,
        charge_length_rules=charge_length_rules,
    )
    logger = ProgressLogger(verbose=verbose, progress=progress)
    logger.log(f"[compactlib] digest preset: {digest_cfg['preset']}")
    with logger.step(f"digesting FASTA {fasta}"):
        df, stats = build_precursor_table_from_fasta(
            fasta=fasta,
            enzyme=digest_cfg["enzyme"],
            missed_cleavages=digest_cfg["missed_cleavages"],
            min_length=digest_cfg["min_length"],
            max_length=digest_cfg["max_length"],
            charges=digest_cfg["charges"],
            min_precursor_mz=digest_cfg["min_precursor_mz"],
            max_precursor_mz=digest_cfg["max_precursor_mz"],
            charge_length_rules=digest_cfg["charge_length_rules"],
            carbamidomethyl_c=carbamidomethyl_c,
            c_mod_format=c_mod_format,
            remove_invalid_aa=remove_invalid_aa,
            allow_selenocysteine=allow_selenocysteine,
            variable_oxidation_m=variable_oxidation_m,
            max_variable_mods=max_variable_mods,
            m_mod_format=m_mod_format,
        )
    logger.log(f"[compactlib] precursor table: {format_shape(df)}")
    if stats.get("n_selenocysteine_precursors_output", 0):
        logger.log("[compactlib] warning: retained selenocysteine (U) precursors. Koina/Prosit/MS2PIP prediction backends generally do not support U; use this mainly for Skyline-compatible precursor audits or filter/drop U before prediction.")

    with logger.step(f"writing precursor table to {output}"):
        write_table(df, output)

    summary_path, params_path = sidecar_paths(output)
    elapsed = logger.total_elapsed()
    stats = {**stats, "preset": digest_cfg["preset"], "compactlib_elapsed_sec": elapsed}

    with logger.step(f"writing summary to {summary_path}"):
        write_summary(pd.DataFrame([stats]), summary_path)

    params = {
        "command": "digest",
        "fasta": str(fasta),
        "output": str(output),
        "preset": digest_cfg["preset"],
        "enzyme": digest_cfg["enzyme"],
        "missed_cleavages": digest_cfg["missed_cleavages"],
        "min_length": digest_cfg["min_length"],
        "max_length": digest_cfg["max_length"],
        "charges": digest_cfg["charges"],
        "min_precursor_mz": digest_cfg["min_precursor_mz"],
        "max_precursor_mz": digest_cfg["max_precursor_mz"],
        "charge_length_rules": digest_cfg["charge_length_rules"],
        "carbamidomethyl_c": carbamidomethyl_c,
        "c_mod_format": c_mod_format,
        "remove_invalid_aa": remove_invalid_aa,
        "allow_selenocysteine": allow_selenocysteine,
        "variable_oxidation_m": variable_oxidation_m,
        "max_variable_mods": max_variable_mods,
        "m_mod_format": m_mod_format,
    }
    with logger.step(f"writing parameters to {params_path}"):
        write_params(params, params_path)

    logger.log(f"[compactlib] output rows: {len(df):,}")
    logger.log(f"[compactlib] total elapsed: {logger.total_elapsed():.2f} s")
    logger.log(f"[compactlib] wrote precursor table: {output}")
    logger.log(f"[compactlib] wrote summary: {summary_path}")
    logger.log(f"[compactlib] wrote params: {params_path}")


@app.command("predict")
def cli_predict(
    precursors: Path = typer.Option(..., "--precursors", "-p", help="Input precursor table produced by compactlib digest"),
    output: Path = typer.Option(..., "--output", "-o", help="Output predicted transition library"),
    backend: str = typer.Option("mock", "--backend", help="Prediction backend: mock or koina"),
    preset: str = typer.Option("none", "--preset", help="Transition filtering preset: none, generic-dia, generic-dia-windowed, generic-dia-strict, diann-like, skyline-dia-like, skyline-export-like"),
    top_n: Optional[int] = typer.Option(None, "--top-n", "-n", min=1, help="Optionally keep only top-N transitions per precursor after prediction"),
    fragment_types: Optional[str] = typer.Option(None, "--fragment-types", help="Override comma-separated fragment ion types, e.g. b,y"),
    fragment_charges: Optional[str] = typer.Option(None, "--fragment-charges", help="Override comma-separated fragment charges, e.g. 1,2"),
    min_fragment_series: Optional[int] = typer.Option(None, "--min-fragment-series", min=1, help="Optional minimum fragment series number applied after prediction, e.g. 3"),
    max_fragment_series: Optional[int] = typer.Option(None, "--max-fragment-series", min=1, help="Optional maximum fragment series number for mock backend and post-prediction filtering"),
    min_product_mz: Optional[float] = typer.Option(None, "--min-product-mz", help="Optional minimum product m/z transition filter"),
    max_product_mz: Optional[float] = typer.Option(None, "--max-product-mz", help="Optional maximum product m/z transition filter"),
    model: Optional[str] = typer.Option(None, "--model", help="Koina model name, e.g. Prosit_2020_intensity_HCD or ms2pip_HCD2021"),
    server_url: str = typer.Option("koina.wilhelmlab.org:443", "--server-url", help="Koina server URL, e.g. koina.wilhelmlab.org:443 or localhost:8500"),
    ssl: bool = typer.Option(True, "--ssl/--no-ssl", help="Use SSL for Koina connection"),
    collision_energy: float = typer.Option(30.0, "--collision-energy", help="Collision energy passed to Koina intensity models"),
    batch_size: int = typer.Option(1024, "--batch-size", min=1, help="Number of precursors per Koina prediction batch"),
    sequence_source: str = typer.Option("stripped", "--sequence-source", help="Koina sequence source: stripped or modified. Use stripped for Prosit_2020_intensity_HCD/ms2pip_HCD2021 unless a model explicitly supports modified notation."),
    koina_sequence_input: str = typer.Option("peptide_sequences", "--koina-sequence-input", help="Koina input column name for peptide sequence"),
    koina_charge_input: str = typer.Option("precursor_charges", "--koina-charge-input", help="Koina input column name for precursor charge"),
    koina_ce_input: str = typer.Option("collision_energies", "--koina-ce-input", help="Koina input column name for collision energy; set empty string to omit"),
    koina_intensity_col: str = typer.Option("intensities", "--koina-intensity-col", help="Koina output column containing intensities"),
    koina_annotation_col: str = typer.Option("annotation", "--koina-annotation-col", help="Koina output column containing fragment annotations"),
    koina_mz_col: Optional[str] = typer.Option(None, "--koina-mz-col", help="Optional Koina output column containing ProductMz values"),
    drop_zero_intensity: bool = typer.Option(True, "--drop-zero-intensity/--keep-zero-intensity", help="Drop non-positive predicted intensities"),
    drop_unsupported_aa: bool = typer.Option(False, "--drop-unsupported-aa/--keep-unsupported-aa", help="Drop U-containing precursors before real prediction backends. Without this, Koina backend stops with a clear error if U is present."),
    include_metadata: bool = typer.Option(False, "--include-metadata/--no-metadata", help="Include all precursor/digest/prediction metadata columns in predicted library output. Default writes a minimal DIA library."),
    include_rt: bool = typer.Option(False, "--include-rt/--no-rt", help="Pass an RT column from precursor table to predicted library output."),
    rt_column: str = typer.Option("NormalizedRetentionTime", "--rt-column", help="RT column name in precursor table when --include-rt is used."),
    rt_output_column: str = typer.Option("NormalizedRetentionTime", "--rt-output-column", help="RT column name written to output library when --include-rt is used."),
    verbose: bool = typer.Option(True, "--verbose/--quiet", help="Print progress logs."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a terminal progress spinner when possible."),
) -> None:
    """Predict a transition library from a precursor table using a backend interface.

    v0.2.2 supports a deterministic mock backend and an optional Koina/koinapy backend.
    Koina model input/output names are configurable because they may differ between models.
    """
    predict_cfg = resolve_predict_preset(
        preset=preset,
        fragment_types=fragment_types,
        fragment_charges=fragment_charges,
        min_fragment_series=min_fragment_series,
        max_fragment_series=max_fragment_series,
        min_product_mz=min_product_mz,
        max_product_mz=max_product_mz,
    )
    logger = ProgressLogger(verbose=verbose, progress=progress)
    logger.log(f"[compactlib] transition/prediction preset: {predict_cfg['preset']}")
    if top_n is None and predict_cfg.get("recommended_top_n") is not None:
        logger.log(f"[compactlib] preset recommends --top-n {predict_cfg['recommended_top_n']} for compact libraries")
    with logger.step(f"reading precursor table from {precursors}"):
        df_prec = read_table(precursors)
    logger.log(f"[compactlib] precursor table: {format_shape(df_prec)}")

    unsupported_mask = _unsupported_aa_mask_for_prediction(df_prec)
    n_unsupported_aa_input = int(unsupported_mask.sum())
    n_unsupported_aa_dropped = 0
    backend_norm = backend.lower().strip()
    if n_unsupported_aa_input:
        msg = (
            f"precursor table contains {n_unsupported_aa_input:,} selenocysteine (U) precursors. "
            "Koina/Prosit/MS2PIP prediction backends generally do not support U."
        )
        if backend_norm != "mock":
            if drop_unsupported_aa:
                logger.log(f"[compactlib] warning: {msg} Dropping them before prediction because --drop-unsupported-aa was provided.")
                df_prec = df_prec.loc[~unsupported_mask].reset_index(drop=True)
                n_unsupported_aa_dropped = n_unsupported_aa_input
                logger.log(f"[compactlib] precursor table after unsupported-AA filtering: {format_shape(df_prec)}")
            else:
                raise typer.BadParameter(
                    msg + " Re-run digest without --allow-selenocysteine, filter U-containing precursors, or pass --drop-unsupported-aa."
                )
        else:
            logger.log(f"[compactlib] warning: {msg} Mock backend can continue because compactlib has a U mass for testing.")

    with logger.step(f"initializing predictor backend: {backend}"):
        predictor = get_predictor(
            backend=backend,
            fragment_types=predict_cfg["fragment_types"],
            fragment_charges=predict_cfg["fragment_charges"],
            max_fragment_series=predict_cfg["max_fragment_series"],
            model=model,
            server_url=server_url,
            ssl=ssl,
            collision_energy=collision_energy,
            batch_size=batch_size,
            sequence_source=sequence_source,
            sequence_input=koina_sequence_input,
            charge_input=koina_charge_input,
            ce_input=koina_ce_input,
            intensity_col=koina_intensity_col,
            annotation_col=koina_annotation_col,
            mz_col=koina_mz_col,
            drop_zero_intensity=drop_zero_intensity,
        )

    with logger.step(f"predicting transitions using backend={backend}"):
        df_pred_raw = predictor.predict(df_prec)
    logger.log(f"[compactlib] raw predicted transition library: {format_shape(df_pred_raw)}")

    with logger.step("applying transition filters"):
        df_pred, filter_stats = filter_transitions(
            df_pred_raw,
            fragment_types=predict_cfg["fragment_types"],
            fragment_charges=predict_cfg["fragment_charges"],
            min_fragment_series=predict_cfg["min_fragment_series"],
            max_fragment_series=predict_cfg["max_fragment_series"],
            min_product_mz=predict_cfg["min_product_mz"],
            max_product_mz=predict_cfg["max_product_mz"],
        )
    logger.log(f"[compactlib] filtered predicted transition library: {format_shape(df_pred)}")

    stats = {
        "backend": backend,
        "preset": predict_cfg["preset"],
        "n_precursors_input": int(len(df_prec) + n_unsupported_aa_dropped),
        "n_precursors_with_unsupported_aa_input": n_unsupported_aa_input,
        "n_precursors_dropped_unsupported_aa": n_unsupported_aa_dropped,
        "n_precursors_predicted_input": int(len(df_prec)),
        "n_transitions_predicted_raw": int(len(df_pred_raw)),
        "n_transitions_predicted": int(len(df_pred)),
        **filter_stats,
    }
    params = {
        "command": "predict",
        "precursors": str(precursors),
        "output": str(output),
        "backend": backend,
        "preset": predict_cfg["preset"],
        "top_n": top_n,
        "fragment_types": predict_cfg["fragment_types"],
        "fragment_charges": predict_cfg["fragment_charges"],
        "min_fragment_series": predict_cfg["min_fragment_series"],
        "max_fragment_series": predict_cfg["max_fragment_series"],
        "min_product_mz": predict_cfg["min_product_mz"],
        "max_product_mz": predict_cfg["max_product_mz"],
        "model": model,
        "server_url": server_url,
        "ssl": ssl,
        "collision_energy": collision_energy,
        "batch_size": batch_size,
        "sequence_source": sequence_source,
        "koina_sequence_input": koina_sequence_input,
        "koina_charge_input": koina_charge_input,
        "koina_ce_input": koina_ce_input,
        "koina_intensity_col": koina_intensity_col,
        "koina_annotation_col": koina_annotation_col,
        "koina_mz_col": koina_mz_col,
        "drop_zero_intensity": drop_zero_intensity,
        "drop_unsupported_aa": drop_unsupported_aa,
        "include_metadata": include_metadata,
        "include_rt": include_rt,
        "rt_column": rt_column,
        "rt_output_column": rt_output_column,
        "output_columns": "all" if include_metadata else _minimal_library_columns(include_rt=include_rt, rt_output_column=rt_output_column),
    }

    if top_n is not None:
        with logger.step(f"selecting top-{top_n} transitions per precursor after prediction"):
            df_out, select_stats = select_top_n(df_pred, top_n=top_n, deduplicate=True)
        stats.update(select_stats)
        mode = f"predict+max-n"
        logger.log(f"[compactlib] compact predicted library: {format_shape(df_out)}")
    else:
        df_out = df_pred
        mode = "predict"

    df_out = _format_predicted_library_output(
        df_out,
        include_metadata=include_metadata,
        include_rt=include_rt,
        rt_column=rt_column,
        rt_output_column=rt_output_column,
    )
    stats["include_metadata"] = include_metadata
    stats["include_rt"] = include_rt
    stats["rt_column"] = rt_column if include_rt else ""
    stats["rt_output_column"] = rt_output_column if include_rt else ""
    stats["n_output_columns"] = len(df_out.columns)
    stats["output_columns"] = "all" if include_metadata else ",".join(_minimal_library_columns(include_rt=include_rt, rt_output_column=rt_output_column))

    _finish(df_pred, df_out, output, mode=mode, params=params, extra_summary=stats, logger=logger)


def _predict_core_from_dataframe(
    df_prec_input: pd.DataFrame,
    *,
    backend: str,
    predict_cfg: dict,
    top_n: Optional[int],
    model: Optional[str],
    server_url: str,
    ssl: bool,
    collision_energy: float,
    batch_size: int,
    sequence_source: str,
    koina_sequence_input: str,
    koina_charge_input: str,
    koina_ce_input: str,
    koina_intensity_col: str,
    koina_annotation_col: str,
    koina_mz_col: Optional[str],
    drop_zero_intensity: bool,
    drop_unsupported_aa: bool,
    include_metadata: bool,
    include_rt: bool,
    rt_column: str,
    rt_output_column: str,
    logger: ProgressLogger,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Run prediction/filtering/top-N selection on an in-memory precursor table.

    Returns (filtered_predicted_transitions, output_transitions, stats).
    This helper is used by chunked prediction; the historical ``predict`` command
    keeps its direct implementation for backwards compatibility.
    """
    df_prec = df_prec_input.copy()

    unsupported_mask = _unsupported_aa_mask_for_prediction(df_prec)
    n_unsupported_aa_input = int(unsupported_mask.sum())
    n_unsupported_aa_dropped = 0
    backend_norm = backend.lower().strip()
    if n_unsupported_aa_input:
        msg = (
            f"precursor table contains {n_unsupported_aa_input:,} selenocysteine (U) precursors. "
            "Koina/Prosit/MS2PIP prediction backends generally do not support U."
        )
        if backend_norm != "mock":
            if drop_unsupported_aa:
                logger.log(f"[compactlib] warning: {msg} Dropping them before prediction because --drop-unsupported-aa was provided.")
                df_prec = df_prec.loc[~unsupported_mask].reset_index(drop=True)
                n_unsupported_aa_dropped = n_unsupported_aa_input
                logger.log(f"[compactlib] precursor table after unsupported-AA filtering: {format_shape(df_prec)}")
            else:
                raise typer.BadParameter(
                    msg + " Re-run digest without --allow-selenocysteine, filter U-containing precursors, or pass --drop-unsupported-aa."
                )
        else:
            logger.log(f"[compactlib] warning: {msg} Mock backend can continue because compactlib has a U mass for testing.")

    with logger.step(f"initializing predictor backend: {backend}"):
        predictor = get_predictor(
            backend=backend,
            fragment_types=predict_cfg["fragment_types"],
            fragment_charges=predict_cfg["fragment_charges"],
            max_fragment_series=predict_cfg["max_fragment_series"],
            model=model,
            server_url=server_url,
            ssl=ssl,
            collision_energy=collision_energy,
            batch_size=batch_size,
            sequence_source=sequence_source,
            sequence_input=koina_sequence_input,
            charge_input=koina_charge_input,
            ce_input=koina_ce_input,
            intensity_col=koina_intensity_col,
            annotation_col=koina_annotation_col,
            mz_col=koina_mz_col,
            drop_zero_intensity=drop_zero_intensity,
        )

    with logger.step(f"predicting transitions using backend={backend}"):
        df_pred_raw = predictor.predict(df_prec)
    logger.log(f"[compactlib] raw predicted transition library: {format_shape(df_pred_raw)}")

    with logger.step("applying transition filters"):
        df_pred, filter_stats = filter_transitions(
            df_pred_raw,
            fragment_types=predict_cfg["fragment_types"],
            fragment_charges=predict_cfg["fragment_charges"],
            min_fragment_series=predict_cfg["min_fragment_series"],
            max_fragment_series=predict_cfg["max_fragment_series"],
            min_product_mz=predict_cfg["min_product_mz"],
            max_product_mz=predict_cfg["max_product_mz"],
        )
    logger.log(f"[compactlib] filtered predicted transition library: {format_shape(df_pred)}")

    stats = {
        "backend": backend,
        "preset": predict_cfg["preset"],
        "n_precursors_input": int(len(df_prec) + n_unsupported_aa_dropped),
        "n_precursors_with_unsupported_aa_input": n_unsupported_aa_input,
        "n_precursors_dropped_unsupported_aa": n_unsupported_aa_dropped,
        "n_precursors_predicted_input": int(len(df_prec)),
        "n_transitions_predicted_raw": int(len(df_pred_raw)),
        "n_transitions_predicted": int(len(df_pred)),
        **filter_stats,
    }

    if top_n is not None:
        with logger.step(f"selecting top-{top_n} transitions per precursor after prediction"):
            df_out, select_stats = select_top_n(df_pred, top_n=top_n, deduplicate=True)
        stats.update(select_stats)
        stats["top_n"] = top_n
        logger.log(f"[compactlib] compact predicted library: {format_shape(df_out)}")
    else:
        df_out = df_pred

    df_out = _format_predicted_library_output(
        df_out,
        include_metadata=include_metadata,
        include_rt=include_rt,
        rt_column=rt_column,
        rt_output_column=rt_output_column,
    )
    stats["include_metadata"] = include_metadata
    stats["include_rt"] = include_rt
    stats["rt_column"] = rt_column if include_rt else ""
    stats["rt_output_column"] = rt_output_column if include_rt else ""
    stats["n_output_columns"] = len(df_out.columns)
    stats["output_columns"] = "all" if include_metadata else ",".join(_minimal_library_columns(include_rt=include_rt, rt_output_column=rt_output_column))

    return df_pred, df_out, stats


def _predict_params_dict(
    *,
    command: str,
    precursors: Path,
    output: Path,
    backend: str,
    predict_cfg: dict,
    top_n: Optional[int],
    model: Optional[str],
    server_url: str,
    ssl: bool,
    collision_energy: float,
    batch_size: int,
    sequence_source: str,
    koina_sequence_input: str,
    koina_charge_input: str,
    koina_ce_input: str,
    koina_intensity_col: str,
    koina_annotation_col: str,
    koina_mz_col: Optional[str],
    drop_zero_intensity: bool,
    drop_unsupported_aa: bool,
    include_metadata: bool,
    include_rt: bool = False,
    rt_column: str = "NormalizedRetentionTime",
    rt_output_column: str = "NormalizedRetentionTime",
    extra: Optional[dict] = None,
) -> dict:
    params = {
        "command": command,
        "precursors": str(precursors),
        "output": str(output),
        "backend": backend,
        "preset": predict_cfg["preset"],
        "top_n": top_n,
        "fragment_types": predict_cfg["fragment_types"],
        "fragment_charges": predict_cfg["fragment_charges"],
        "min_fragment_series": predict_cfg["min_fragment_series"],
        "max_fragment_series": predict_cfg["max_fragment_series"],
        "min_product_mz": predict_cfg["min_product_mz"],
        "max_product_mz": predict_cfg["max_product_mz"],
        "model": model,
        "server_url": server_url,
        "ssl": ssl,
        "collision_energy": collision_energy,
        "batch_size": batch_size,
        "sequence_source": sequence_source,
        "koina_sequence_input": koina_sequence_input,
        "koina_charge_input": koina_charge_input,
        "koina_ce_input": koina_ce_input,
        "koina_intensity_col": koina_intensity_col,
        "koina_annotation_col": koina_annotation_col,
        "koina_mz_col": koina_mz_col,
        "drop_zero_intensity": drop_zero_intensity,
        "drop_unsupported_aa": drop_unsupported_aa,
        "include_metadata": include_metadata,
        "include_rt": include_rt,
        "rt_column": rt_column,
        "rt_output_column": rt_output_column,
        "output_columns": "all" if include_metadata else _minimal_library_columns(include_rt=include_rt, rt_output_column=rt_output_column),
    }
    if extra:
        params.update(extra)
    return params


def _output_stem_without_compound_suffix(path: Path) -> str:
    name = path.name
    for suffix in [".tsv.gz", ".csv.gz", ".tsv", ".csv", ".parquet"]:
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def _merge_table_chunks(chunk_paths: list[Path], output: Path, output_columns: Optional[list[str]] = None) -> None:
    """Stream-merge chunk tables into final output. Supports TSV/CSV plain/gz.

    When ``output_columns`` is provided, only these columns are written. This is
    useful for re-merging old full-metadata chunks into a minimal final library.
    """
    import gzip
    import shutil

    output = Path(output).expanduser()
    output.parent.mkdir(parents=True, exist_ok=True)
    name = output.name.lower()

    if name.endswith(".parquet"):
        raise typer.BadParameter(
            "Chunked final merge to parquet is intentionally not supported for large libraries. "
            "Use .tsv/.tsv.gz/.csv/.csv.gz as --output."
        )

    if name.endswith(".tsv") or name.endswith(".tsv.gz"):
        out_sep = "\t"
    elif name.endswith(".csv") or name.endswith(".csv.gz"):
        out_sep = ","
    else:
        raise typer.BadParameter("Chunked output must be .tsv, .tsv.gz, .csv or .csv.gz")

    tmp = output.with_name(
        _output_stem_without_compound_suffix(output)
        + ".tmp"
        + ''.join(output.suffixes[-2:] if output.name.lower().endswith((".tsv.gz", ".csv.gz")) else output.suffixes[-1:])
    )

    if output_columns is None and (name.endswith(".tsv") or name.endswith(".tsv.gz")):
        opener_out = gzip.open if name.endswith(".gz") else open
        with opener_out(tmp, "wt", encoding="utf-8") as fout:
            wrote_header = False
            for path in chunk_paths:
                opener_in = gzip.open if path.name.lower().endswith(".gz") else open
                with opener_in(path, "rt", encoding="utf-8") as fin:
                    header = fin.readline()
                    if not header:
                        continue
                    if not wrote_header:
                        fout.write(header)
                        wrote_header = True
                    shutil.copyfileobj(fin, fout)
    else:
        import csv

        opener_out = gzip.open if name.endswith(".gz") else open
        with opener_out(tmp, "wt", encoding="utf-8", newline="") as fout:
            writer = csv.writer(fout, delimiter=out_sep, lineterminator="\n")
            wrote_header = False
            selected_indices: Optional[list[int]] = None

            for path in chunk_paths:
                opener_in = gzip.open if path.name.lower().endswith(".gz") else open
                with opener_in(path, "rt", encoding="utf-8", newline="") as fin:
                    reader = csv.reader(fin, delimiter="\t")
                    try:
                        header = next(reader)
                    except StopIteration:
                        continue

                    if output_columns is None:
                        current_columns = header
                        indices = list(range(len(header)))
                    else:
                        missing = [col for col in output_columns if col not in header]
                        if missing:
                            raise typer.BadParameter(
                                "Cannot merge chunk into requested output profile because columns are missing in "
                                f"{path}: " + ", ".join(missing)
                            )
                        current_columns = output_columns
                        indices = [header.index(col) for col in output_columns]

                    if selected_indices is None:
                        selected_indices = indices
                    if not wrote_header:
                        writer.writerow(current_columns)
                        wrote_header = True

                    for row in reader:
                        writer.writerow([row[i] if i < len(row) else "" for i in indices])

    os.replace(tmp, output)


def _completed_chunk_marker(chunk_output: Path) -> Path:
    return chunk_output.with_name(chunk_output.name + ".done.json")


def _write_done_marker(marker: Path, payload: dict) -> None:
    marker.parent.mkdir(parents=True, exist_ok=True)
    with open(marker, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False, default=str)


def _all_chunks_done(chunk_outputs: list[Path]) -> bool:
    return all(path.exists() and _completed_chunk_marker(path).exists() for path in chunk_outputs)


def _final_chunked_summary(
    *,
    chunk_summaries: pd.DataFrame,
    output: Path,
    mode: str,
    params_extra: dict,
) -> pd.DataFrame:
    row = {
        "compactlib_version": __version__,
        "mode": mode,
        "output_file": str(output),
    }
    for col in [
        "n_precursors_input",
        "n_precursors_with_unsupported_aa_input",
        "n_precursors_dropped_unsupported_aa",
        "n_precursors_predicted_input",
        "n_transitions_predicted_raw",
        "n_transitions_predicted",
        "n_transitions_before_transition_filters",
        "n_transitions_filtered_by_transition_filters",
        "n_transitions_after_transition_filters",
        "n_duplicate_transitions_removed",
        "n_precursors_with_less_than_n",
        "n_precursors_output",
        "n_transitions_output",
    ]:
        if col in chunk_summaries.columns:
            row[col] = pd.to_numeric(chunk_summaries[col], errors="coerce").fillna(0).sum()

    if row.get("n_precursors_output", 0):
        row["mean_transitions_per_precursor"] = row.get("n_transitions_output", 0) / row["n_precursors_output"]
    else:
        row["mean_transitions_per_precursor"] = 0.0

    if "min_transitions_per_precursor" in chunk_summaries.columns:
        row["min_transitions_per_precursor"] = int(pd.to_numeric(chunk_summaries["min_transitions_per_precursor"], errors="coerce").min())
    if "max_transitions_per_precursor" in chunk_summaries.columns:
        row["max_transitions_per_precursor"] = int(pd.to_numeric(chunk_summaries["max_transitions_per_precursor"], errors="coerce").max())
    if "median_transitions_per_precursor" in chunk_summaries.columns:
        # Exact global median would require per-precursor counts. For top-N libraries
        # chunk medians are normally stable; keep the median of chunk medians.
        row["median_transitions_per_precursor"] = float(pd.to_numeric(chunk_summaries["median_transitions_per_precursor"], errors="coerce").median())

    row.update(params_extra)
    output = Path(output).expanduser()
    row["library_size_mb"] = os.path.getsize(output) / 1024**2 if output.exists() else None
    return pd.DataFrame([row])


@app.command("predict-chunked")
def cli_predict_chunked(
    precursors: Path = typer.Option(..., "--precursors", "-p", help="Input precursor table produced by compactlib digest"),
    output: Path = typer.Option(..., "--output", "-o", help="Final merged predicted transition library"),
    work_dir: Optional[Path] = typer.Option(None, "--work-dir", help="Directory for precursor chunks, predicted chunk outputs and manifests. Default: <output>.chunked_work"),
    chunk_size: int = typer.Option(100_000, "--chunk-size", min=1, help="Number of precursor rows per persistent outer chunk"),
    resume: bool = typer.Option(True, "--resume/--no-resume", help="Skip chunks with completed output and .done marker"),
    force: bool = typer.Option(False, "--force", help="Recompute chunks even if completed markers exist"),
    backend: str = typer.Option("mock", "--backend", help="Prediction backend: mock or koina"),
    preset: str = typer.Option("none", "--preset", help="Transition filtering preset: none, generic-dia, generic-dia-windowed, generic-dia-strict, diann-like, skyline-dia-like, skyline-export-like"),
    top_n: Optional[int] = typer.Option(None, "--top-n", "-n", min=1, help="Optionally keep only top-N transitions per precursor after prediction"),
    fragment_types: Optional[str] = typer.Option(None, "--fragment-types", help="Override comma-separated fragment ion types, e.g. b,y"),
    fragment_charges: Optional[str] = typer.Option(None, "--fragment-charges", help="Override comma-separated fragment charges, e.g. 1,2"),
    min_fragment_series: Optional[int] = typer.Option(None, "--min-fragment-series", min=1, help="Optional minimum fragment series number applied after prediction, e.g. 3"),
    max_fragment_series: Optional[int] = typer.Option(None, "--max-fragment-series", min=1, help="Optional maximum fragment series number for mock backend and post-prediction filtering"),
    min_product_mz: Optional[float] = typer.Option(None, "--min-product-mz", help="Optional minimum product m/z transition filter"),
    max_product_mz: Optional[float] = typer.Option(None, "--max-product-mz", help="Optional maximum product m/z transition filter"),
    model: Optional[str] = typer.Option(None, "--model", help="Koina model name, e.g. Prosit_2020_intensity_HCD or ms2pip_HCD2021"),
    server_url: str = typer.Option("koina.wilhelmlab.org:443", "--server-url", help="Koina server URL, e.g. koina.wilhelmlab.org:443 or localhost:8500"),
    ssl: bool = typer.Option(True, "--ssl/--no-ssl", help="Use SSL for Koina connection"),
    collision_energy: float = typer.Option(30.0, "--collision-energy", help="Collision energy passed to Koina intensity models"),
    batch_size: int = typer.Option(1024, "--batch-size", min=1, help="Number of precursors per Koina request inside each persistent chunk"),
    sequence_source: str = typer.Option("stripped", "--sequence-source", help="Koina sequence source: stripped or modified"),
    koina_sequence_input: str = typer.Option("peptide_sequences", "--koina-sequence-input", help="Koina input column name for peptide sequence"),
    koina_charge_input: str = typer.Option("precursor_charges", "--koina-charge-input", help="Koina input column name for precursor charge"),
    koina_ce_input: str = typer.Option("collision_energies", "--koina-ce-input", help="Koina input column name for collision energy; set empty string to omit"),
    koina_intensity_col: str = typer.Option("intensities", "--koina-intensity-col", help="Koina output column containing intensities"),
    koina_annotation_col: str = typer.Option("annotation", "--koina-annotation-col", help="Koina output column containing fragment annotations"),
    koina_mz_col: Optional[str] = typer.Option(None, "--koina-mz-col", help="Optional Koina output column containing ProductMz values"),
    drop_zero_intensity: bool = typer.Option(True, "--drop-zero-intensity/--keep-zero-intensity", help="Drop non-positive predicted intensities"),
    drop_unsupported_aa: bool = typer.Option(False, "--drop-unsupported-aa/--keep-unsupported-aa", help="Drop U-containing precursors before real prediction backends"),
    include_metadata: bool = typer.Option(False, "--include-metadata/--no-metadata", help="Include all precursor/digest/prediction metadata columns in predicted chunk/final outputs. Default writes a minimal DIA library."),
    include_rt: bool = typer.Option(False, "--include-rt/--no-rt", help="Pass an RT column from precursor chunks to predicted chunk/final outputs."),
    rt_column: str = typer.Option("NormalizedRetentionTime", "--rt-column", help="RT column name in precursor table when --include-rt is used."),
    rt_output_column: str = typer.Option("NormalizedRetentionTime", "--rt-output-column", help="RT column name written to output library when --include-rt is used."),
    verbose: bool = typer.Option(True, "--verbose/--quiet", help="Print progress logs."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a terminal progress spinner when possible."),
) -> None:
    """Predict a library using persistent outer chunks with automatic resume.

    The command writes per-chunk predicted libraries to a work directory. A chunk
    is considered complete only when both its output table and .done marker are
    present. Re-running the same command with --resume skips completed chunks.
    The final output file is written only after all chunks are complete.
    """
    predict_cfg = resolve_predict_preset(
        preset=preset,
        fragment_types=fragment_types,
        fragment_charges=fragment_charges,
        min_fragment_series=min_fragment_series,
        max_fragment_series=max_fragment_series,
        min_product_mz=min_product_mz,
        max_product_mz=max_product_mz,
    )

    logger = ProgressLogger(verbose=verbose, progress=progress)
    logger.log(f"[compactlib] chunked transition/prediction preset: {predict_cfg['preset']}")
    if top_n is None and predict_cfg.get("recommended_top_n") is not None:
        logger.log(f"[compactlib] preset recommends --top-n {predict_cfg['recommended_top_n']} for compact libraries")

    output = Path(output).expanduser()
    if work_dir is None:
        work_dir = output.parent / f"{_output_stem_without_compound_suffix(output)}.chunked_work"
    else:
        work_dir = Path(work_dir).expanduser()

    precursor_chunk_dir = work_dir / "precursor_chunks"
    predicted_chunk_dir = work_dir / "predicted_chunks"
    summary_chunk_dir = work_dir / "chunk_summaries"
    params_chunk_dir = work_dir / "chunk_params"
    for d in [precursor_chunk_dir, predicted_chunk_dir, summary_chunk_dir, params_chunk_dir]:
        d.mkdir(parents=True, exist_ok=True)

    with logger.step(f"reading precursor table from {precursors}"):
        df_prec_all = read_table(precursors)
    logger.log(f"[compactlib] precursor table: {format_shape(df_prec_all)}")

    n_total = len(df_prec_all)
    n_chunks = (n_total + chunk_size - 1) // chunk_size if n_total else 0
    logger.log(f"[compactlib] chunk size: {chunk_size:,}; chunks: {n_chunks:,}")
    logger.log(f"[compactlib] work dir: {work_dir}")

    command_params = _predict_params_dict(
        command="predict-chunked",
        precursors=precursors,
        output=output,
        backend=backend,
        predict_cfg=predict_cfg,
        top_n=top_n,
        model=model,
        server_url=server_url,
        ssl=ssl,
        collision_energy=collision_energy,
        batch_size=batch_size,
        sequence_source=sequence_source,
        koina_sequence_input=koina_sequence_input,
        koina_charge_input=koina_charge_input,
        koina_ce_input=koina_ce_input,
        koina_intensity_col=koina_intensity_col,
        koina_annotation_col=koina_annotation_col,
        koina_mz_col=koina_mz_col,
        drop_zero_intensity=drop_zero_intensity,
        drop_unsupported_aa=drop_unsupported_aa,
        include_metadata=include_metadata,
        include_rt=include_rt,
        rt_column=rt_column,
        rt_output_column=rt_output_column,
        extra={
            "chunk_size": chunk_size,
            "work_dir": str(work_dir),
            "resume": resume,
            "force": force,
            "n_input_rows": n_total,
            "n_chunks": n_chunks,
        },
    )
    write_params(command_params, work_dir / "predict_chunked.params.json")

    chunk_records = []
    chunk_outputs: list[Path] = []

    for i, start in enumerate(range(0, n_total, chunk_size), start=1):
        stop = min(start + chunk_size, n_total)
        chunk_id = f"chunk_{i:05d}"
        chunk_prec_path = precursor_chunk_dir / f"{chunk_id}.precursors.tsv"
        chunk_out = predicted_chunk_dir / f"{chunk_id}.predicted.tsv"
        chunk_summary_path = summary_chunk_dir / f"{chunk_id}.summary.tsv"
        chunk_params_path = params_chunk_dir / f"{chunk_id}.params.json"
        done_marker = _completed_chunk_marker(chunk_out)
        chunk_outputs.append(chunk_out)

        record = {
            "chunk_index": i,
            "chunk_id": chunk_id,
            "row_start_0based": start,
            "row_stop_0based_exclusive": stop,
            "n_precursors": stop - start,
            "precursor_chunk_file": str(chunk_prec_path),
            "predicted_chunk_file": str(chunk_out),
            "done_marker": str(done_marker),
        }
        chunk_records.append(record)

        if force:
            for pth in [chunk_out, chunk_summary_path, chunk_params_path, done_marker]:
                if pth.exists():
                    pth.unlink()

        if resume and chunk_out.exists() and done_marker.exists():
            logger.log(f"[compactlib] [skip] {chunk_id}: completed chunk exists")
            continue

        df_chunk = df_prec_all.iloc[start:stop].reset_index(drop=True)
        with logger.step(f"writing precursor chunk {chunk_id} to {chunk_prec_path}"):
            write_table(df_chunk, chunk_prec_path)

        logger.log(f"[compactlib] [run] {chunk_id}: rows {start:,}:{stop:,} ({len(df_chunk):,} precursors)")
        tmp_out = predicted_chunk_dir / f"{chunk_id}.predicted.tmp.tsv"
        if tmp_out.exists():
            tmp_out.unlink()

        chunk_logger = logger
        df_pred, df_out, stats = _predict_core_from_dataframe(
            df_chunk,
            backend=backend,
            predict_cfg=predict_cfg,
            top_n=top_n,
            model=model,
            server_url=server_url,
            ssl=ssl,
            collision_energy=collision_energy,
            batch_size=batch_size,
            sequence_source=sequence_source,
            koina_sequence_input=koina_sequence_input,
            koina_charge_input=koina_charge_input,
            koina_ce_input=koina_ce_input,
            koina_intensity_col=koina_intensity_col,
            koina_annotation_col=koina_annotation_col,
            koina_mz_col=koina_mz_col,
            drop_zero_intensity=drop_zero_intensity,
            drop_unsupported_aa=drop_unsupported_aa,
            include_metadata=include_metadata,
            include_rt=include_rt,
            rt_column=rt_column,
            rt_output_column=rt_output_column,
            logger=chunk_logger,
        )

        with logger.step(f"writing predicted chunk {chunk_id} to {chunk_out}"):
            write_library(df_out, tmp_out)
            os.replace(tmp_out, chunk_out)

        mode = "predict-chunked+max-n" if top_n is not None else "predict-chunked"
        summary = library_summary(
            df_input=df_pred,
            df_output=df_out,
            output_path=chunk_out,
            mode=mode,
            extra={
                **stats,
                "chunk_index": i,
                "chunk_id": chunk_id,
                "chunk_row_start_0based": start,
                "chunk_row_stop_0based_exclusive": stop,
            },
        )
        write_summary(summary, chunk_summary_path)

        chunk_params = {**command_params, "output": str(chunk_out), "chunk_index": i, "chunk_id": chunk_id}
        write_params(chunk_params, chunk_params_path)
        _write_done_marker(
            done_marker,
            {
                "chunk_id": chunk_id,
                "chunk_index": i,
                "output": str(chunk_out),
                "summary": str(chunk_summary_path),
                "n_precursors": len(df_chunk),
                "n_transitions_output": len(df_out),
            },
        )
        logger.log(f"[compactlib] [done] {chunk_id}: output rows={len(df_out):,}")

    manifest = pd.DataFrame(chunk_records)
    manifest["completed"] = [Path(r["predicted_chunk_file"]).exists() and Path(r["done_marker"]).exists() for r in chunk_records]
    write_summary(manifest, work_dir / "chunks_manifest.tsv")

    if not _all_chunks_done(chunk_outputs):
        n_done = sum(path.exists() and _completed_chunk_marker(path).exists() for path in chunk_outputs)
        raise RuntimeError(
            f"Only {n_done}/{len(chunk_outputs)} chunks are complete. Final output was not merged. "
            f"Re-run the same command with --resume to continue. Work dir: {work_dir}"
        )

    with logger.step(f"merging {len(chunk_outputs):,} completed chunks into final output {output}"):
        _merge_table_chunks(
            chunk_outputs,
            output,
            output_columns=None if include_metadata else _minimal_library_columns(include_rt=include_rt, rt_output_column=rt_output_column),
        )

    summary_paths = [summary_chunk_dir / f"chunk_{i:05d}.summary.tsv" for i in range(1, n_chunks + 1)]
    chunk_summaries = pd.concat([pd.read_csv(p, sep="\t") for p in summary_paths], ignore_index=True) if summary_paths else pd.DataFrame()
    final_summary = _final_chunked_summary(
        chunk_summaries=chunk_summaries,
        output=output,
        mode="predict-chunked+max-n" if top_n is not None else "predict-chunked",
        params_extra={
            "backend": backend,
            "preset": predict_cfg["preset"],
            "model": model,
            "top_n": top_n,
            "chunk_size": chunk_size,
            "n_chunks": n_chunks,
            "work_dir": str(work_dir),
            "include_metadata": include_metadata,
            "include_rt": include_rt,
            "rt_column": rt_column if include_rt else "",
            "rt_output_column": rt_output_column if include_rt else "",
            "output_columns": "all" if include_metadata else ",".join(_minimal_library_columns(include_rt=include_rt, rt_output_column=rt_output_column)),
            "compactlib_elapsed_sec": logger.total_elapsed(),
        },
    )

    final_summary_path, final_params_path = sidecar_paths(output)
    with logger.step(f"writing final summary to {final_summary_path}"):
        write_summary(final_summary, final_summary_path)
    with logger.step(f"writing final parameters to {final_params_path}"):
        write_params(command_params, final_params_path)

    logger.log(f"[compactlib] final output rows: {int(final_summary.loc[0, 'n_transitions_output']) if 'n_transitions_output' in final_summary.columns else 'NA'}")
    logger.log(f"[compactlib] total elapsed: {logger.total_elapsed():.2f} s")
    logger.log(f"[compactlib] wrote final library: {output}")
    logger.log(f"[compactlib] wrote final summary: {final_summary_path}")
    logger.log(f"[compactlib] wrote chunk manifest: {work_dir / 'chunks_manifest.tsv'}")



def _sample_precursors_stratified(
    df: pd.DataFrame,
    *,
    sample_size: Optional[int],
    charge_column: str = "PrecursorCharge",
    sequence_column: str = "StrippedPeptide",
    random_state: int = 42,
    drop_unsupported_aa: bool = True,
) -> tuple[pd.DataFrame, dict]:
    """Return predictor-safe optional subset without losing grouping columns.

    This helper is intentionally explicit rather than using groupby.apply,
    because pandas versions differ in whether grouping columns are retained.
    """
    out = df.copy()
    out.columns = [str(c).strip() for c in out.columns]

    required = [charge_column, sequence_column]
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise typer.BadParameter(
            f"Cannot sample precursors because required columns are missing: {missing}. "
            f"Available columns: {list(out.columns)}"
        )

    n_input = len(out)
    n_unsupported = int(_unsupported_aa_mask_for_prediction(out).sum())
    if drop_unsupported_aa and n_unsupported:
        out = out.loc[~_unsupported_aa_mask_for_prediction(out)].reset_index(drop=True)

    n_after_unsupported = len(out)
    if sample_size is None or sample_size <= 0 or sample_size >= len(out):
        return out.reset_index(drop=True), {
            "n_precursors_before_sampling": n_input,
            "n_precursors_with_unsupported_aa_before_sampling": n_unsupported,
            "n_precursors_after_unsupported_aa_filter": n_after_unsupported,
            "sample_size_requested": sample_size,
            "n_precursors_after_sampling": len(out),
            "sampling_mode": "none" if sample_size is None else "all_available",
            "sampling_random_state": random_state,
        }

    charges = list(pd.Series(out[charge_column].dropna().unique()).sort_values())
    if not charges:
        raise typer.BadParameter("Cannot sample precursor table: no precursor charges found.")

    base = sample_size // len(charges)
    remainder = sample_size % len(charges)
    selected_indices: list[int] = []
    selected_set: set[int] = set()

    for i, charge in enumerate(charges):
        sub = out.loc[out[charge_column] == charge]
        quota = base + (1 if i < remainder else 0)
        n = min(quota, len(sub))
        if n > 0:
            sampled = sub.sample(n=n, random_state=random_state)
            idxs = list(sampled.index)
            selected_indices.extend(idxs)
            selected_set.update(int(x) for x in idxs)

    # If some charge strata had fewer rows than their quota, fill the deficit
    # from remaining precursors while preserving reproducibility.
    deficit = sample_size - len(selected_indices)
    if deficit > 0:
        remaining = out.loc[[idx for idx in out.index if int(idx) not in selected_set]]
        if len(remaining):
            fill = remaining.sample(n=min(deficit, len(remaining)), random_state=random_state)
            selected_indices.extend(list(fill.index))

    sampled = out.loc[selected_indices].copy().reset_index(drop=True)
    sampled = sampled[out.columns]
    return sampled, {
        "n_precursors_before_sampling": n_input,
        "n_precursors_with_unsupported_aa_before_sampling": n_unsupported,
        "n_precursors_after_unsupported_aa_filter": n_after_unsupported,
        "sample_size_requested": sample_size,
        "n_precursors_after_sampling": len(sampled),
        "sampling_mode": "stratified_by_precursor_charge",
        "sampling_random_state": random_state,
    }


@app.command("build")
def cli_build(
    fasta: Path = typer.Option(..., "--fasta", "-f", help="Input FASTA file"),
    output_dir: Path = typer.Option(..., "--output-dir", "-o", help="Directory for all generated files"),
    prefix: str = typer.Option("compactlib", "--prefix", help="Output file prefix"),
    work_dir: Optional[Path] = typer.Option(None, "--work-dir", help="Working directory. Default: <output-dir>/<prefix>.build_work"),
    preset: str = typer.Option("generic-dia", "--preset", help="Digest and transition filtering preset"),
    top_n: int = typer.Option(7, "--top-n", "-n", min=1, help="Top-N transitions per model; union output has up to 2N transitions"),
    chunk_size: int = typer.Option(100_000, "--chunk-size", min=1, help="Outer chunk size for resumable prediction"),
    batch_size: int = typer.Option(512, "--batch-size", min=1, help="Koina request batch size inside each chunk"),
    sample_size: Optional[int] = typer.Option(None, "--sample-size", min=1, help="Optional predictor-safe stratified precursor subset size for smoke tests"),
    random_state: int = typer.Option(42, "--random-state", help="Random seed for --sample-size"),
    predict_rt: bool = typer.Option(True, "--rt/--no-rt", help="Predict/pass RT into final libraries"),
    rt_backend: str = typer.Option("deeplc", "--rt-backend", help="RT backend: deeplc or mock"),
    rt_output_column: str = typer.Option("PredictedRetentionTime", "--rt-output-column", help="RT column written to precursor and library tables"),
    calibration_table: Optional[Path] = typer.Option(None, "--calibration-table", help="Optional calibration table for DeepLC"),
    calibration_rt_column: str = typer.Option("tr", "--calibration-rt-column", help="RT column in calibration table"),
    deeplc_model_path: Optional[str] = typer.Option(None, "--deeplc-model-path", help="Optional DeepLC model path"),
    prediction_backend: str = typer.Option("koina", "--prediction-backend", help="Fragment intensity backend: koina or mock"),
    prosit_model: str = typer.Option("Prosit_2020_intensity_HCD", "--prosit-model", help="Koina Prosit model name"),
    ms2pip_model: str = typer.Option("ms2pip_HCD2021", "--ms2pip-model", help="Koina MS2PIP model name"),
    server_url: str = typer.Option("koina.wilhelmlab.org:443", "--server-url", help="Koina server URL"),
    ssl: bool = typer.Option(True, "--ssl/--no-ssl", help="Use SSL for Koina connection"),
    collision_energy: float = typer.Option(30.0, "--collision-energy", help="Collision energy passed to Koina models"),
    sequence_source: str = typer.Option("stripped", "--sequence-source", help="Koina sequence source: stripped or modified"),
    union: bool = typer.Option(True, "--union/--no-union", help="Build Union2N from Prosit and MS2PIP outputs"),
    resume: bool = typer.Option(True, "--resume/--no-resume", help="Reuse completed intermediate files/chunks"),
    force: bool = typer.Option(False, "--force", help="Recompute intermediates even if they already exist"),
    verbose: bool = typer.Option(True, "--verbose/--quiet", help="Print progress logs"),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show terminal progress spinner when possible"),
) -> None:
    """One-command FASTA-to-library workflow.

    Builds a digest precursor table, optionally predicts RT, predicts Prosit maxN
    and MS2PIP maxN libraries using resumable chunks, and optionally builds a
    Union2N library. The default outputs are minimal DIA-library tables.
    """
    logger = ProgressLogger(verbose=verbose, progress=progress)
    output_dir = Path(output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = Path(work_dir).expanduser() if work_dir is not None else output_dir / f"{prefix}.build_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    logger.log("[compactlib] running one-command build workflow")
    logger.log(f"[compactlib] output dir: {output_dir}")
    logger.log(f"[compactlib] work dir: {work_dir}")

    # ------------------------------------------------------------------
    # 1. Digest FASTA
    # ------------------------------------------------------------------
    precursors_all = output_dir / f"{prefix}.precursors.all.tsv"
    if resume and precursors_all.exists() and not force:
        logger.log(f"[compactlib] [skip] digest: {precursors_all} exists")
        df_prec_all = read_table(precursors_all)
        digest_stats = {"digest_skipped_existing": True, "n_precursors_digest_output": len(df_prec_all)}
    else:
        digest_cfg = resolve_digest_preset(
            preset=preset,
            enzyme=None,
            missed_cleavages=None,
            min_length=None,
            max_length=None,
            charges=None,
            min_precursor_mz=None,
            max_precursor_mz=None,
            charge_length_rules=None,
        )
        logger.log(f"[compactlib] digest preset: {digest_cfg['preset']}")
        with logger.step(f"digesting FASTA {fasta}"):
            df_prec_all, digest_stats = build_precursor_table_from_fasta(
                fasta=fasta,
                enzyme=digest_cfg["enzyme"],
                missed_cleavages=digest_cfg["missed_cleavages"],
                min_length=digest_cfg["min_length"],
                max_length=digest_cfg["max_length"],
                charges=digest_cfg["charges"],
                min_precursor_mz=digest_cfg["min_precursor_mz"],
                max_precursor_mz=digest_cfg["max_precursor_mz"],
                charge_length_rules=digest_cfg["charge_length_rules"],
                carbamidomethyl_c=True,
                c_mod_format="unimod",
                remove_invalid_aa=True,
                allow_selenocysteine=False,
                variable_oxidation_m=False,
                max_variable_mods=1,
                m_mod_format="unimod",
            )
        with logger.step(f"writing full precursor table to {precursors_all}"):
            write_table(df_prec_all, precursors_all)
        write_summary(pd.DataFrame([{**digest_stats, "preset": preset}]), precursors_all.with_suffix(".summary.tsv"))
        write_params({"command": "build:digest", "fasta": str(fasta), "output": str(precursors_all), "preset": preset}, precursors_all.with_suffix(".params.json"))

    # ------------------------------------------------------------------
    # 2. Optional predictor-safe subset for smoke tests
    # ------------------------------------------------------------------
    if sample_size is not None:
        active_precursors = output_dir / f"{prefix}.precursors.sample_{sample_size}.tsv"
        if resume and active_precursors.exists() and not force:
            logger.log(f"[compactlib] [skip] sampling: {active_precursors} exists")
            df_active = read_table(active_precursors)
            sample_stats = {"sampling_skipped_existing": True, "sample_size_requested": sample_size, "n_precursors_after_sampling": len(df_active)}
        else:
            with logger.step(f"creating predictor-safe stratified sample of {sample_size:,} precursors"):
                df_active, sample_stats = _sample_precursors_stratified(
                    df_prec_all,
                    sample_size=sample_size,
                    random_state=random_state,
                    drop_unsupported_aa=True,
                )
            with logger.step(f"writing active precursor table to {active_precursors}"):
                write_table(df_active, active_precursors)
            write_summary(pd.DataFrame([sample_stats]), active_precursors.with_suffix(".summary.tsv"))
            write_params({"command": "build:sample", "input": str(precursors_all), "output": str(active_precursors), **sample_stats}, active_precursors.with_suffix(".params.json"))
    else:
        active_precursors = output_dir / f"{prefix}.precursors.tsv"
        if resume and active_precursors.exists() and not force:
            logger.log(f"[compactlib] [skip] active precursor table: {active_precursors} exists")
            df_active = read_table(active_precursors)
            sample_stats = {"sampling_mode": "none", "n_precursors_after_sampling": len(df_active)}
        else:
            with logger.step("creating predictor-safe active precursor table"):
                df_active, sample_stats = _sample_precursors_stratified(
                    df_prec_all,
                    sample_size=None,
                    random_state=random_state,
                    drop_unsupported_aa=True,
                )
            with logger.step(f"writing active precursor table to {active_precursors}"):
                write_table(df_active, active_precursors)
            write_summary(pd.DataFrame([sample_stats]), active_precursors.with_suffix(".summary.tsv"))
            write_params({"command": "build:active-precursors", "input": str(precursors_all), "output": str(active_precursors), **sample_stats}, active_precursors.with_suffix(".params.json"))

    # ------------------------------------------------------------------
    # 3. Optional RT prediction
    # ------------------------------------------------------------------
    if predict_rt:
        precursors_for_prediction = output_dir / f"{prefix}.precursors.with_rt.tsv"
        if sample_size is not None:
            precursors_for_prediction = output_dir / f"{prefix}.precursors.sample_{sample_size}.with_rt.tsv"
        if resume and precursors_for_prediction.exists() and not force:
            logger.log(f"[compactlib] [skip] RT prediction: {precursors_for_prediction} exists")
        else:
            df_cal = None
            if calibration_table is not None:
                with logger.step(f"reading DeepLC calibration table from {calibration_table}"):
                    df_cal = read_table(calibration_table)
            with logger.step(f"predicting RT using backend={rt_backend}"):
                rt_result = add_deeplc_rt_to_precursors(
                    df_active,
                    backend=rt_backend,
                    modified_column="ModifiedPeptide",
                    stripped_column="StrippedPeptide",
                    rt_output_column=rt_output_column,
                    fixed_carbamidomethyl_c=True,
                    calibration_table=df_cal,
                    calibration_rt_column=calibration_rt_column,
                    deeplc_model_path=deeplc_model_path,
                )
            with logger.step(f"writing precursor table with RT to {precursors_for_prediction}"):
                write_table(rt_result.precursor_table, precursors_for_prediction)
            write_summary(pd.DataFrame([{**rt_result.stats, "compactlib_elapsed_sec": logger.total_elapsed()}]), precursors_for_prediction.with_suffix(".summary.tsv"))
            write_params({
                "command": "build:predict-rt",
                "precursors": str(active_precursors),
                "output": str(precursors_for_prediction),
                "backend": rt_backend,
                "rt_output_column": rt_output_column,
                "calibration_table": str(calibration_table) if calibration_table else None,
                "calibration_rt_column": calibration_rt_column,
                "deeplc_model_path": deeplc_model_path,
            }, precursors_for_prediction.with_suffix(".params.json"))
        include_rt_flag = True
        rt_column_for_prediction = rt_output_column
    else:
        precursors_for_prediction = active_precursors
        include_rt_flag = False
        rt_column_for_prediction = rt_output_column

    # ------------------------------------------------------------------
    # 4. Predict Prosit/MS2PIP maxN using resumable chunked prediction
    # ------------------------------------------------------------------
    prosit_out = output_dir / f"{prefix}.Prosit_max{top_n}.tsv"
    ms2pip_out = output_dir / f"{prefix}.MS2PIP_max{top_n}.tsv"

    logger.log("[compactlib] building Prosit compact library")
    cli_predict_chunked(
        precursors=precursors_for_prediction,
        output=prosit_out,
        work_dir=work_dir / f"Prosit_max{top_n}.work",
        chunk_size=chunk_size,
        resume=resume,
        force=force,
        backend=prediction_backend,
        preset=preset,
        top_n=top_n,
        fragment_types=None,
        fragment_charges=None,
        min_fragment_series=None,
        max_fragment_series=None,
        min_product_mz=None,
        max_product_mz=None,
        model=prosit_model if prediction_backend.lower() == "koina" else None,
        server_url=server_url,
        ssl=ssl,
        collision_energy=collision_energy,
        batch_size=batch_size,
        sequence_source=sequence_source,
        koina_sequence_input="peptide_sequences",
        koina_charge_input="precursor_charges",
        koina_ce_input="collision_energies",
        koina_intensity_col="intensities",
        koina_annotation_col="annotation",
        koina_mz_col=None,
        drop_zero_intensity=True,
        drop_unsupported_aa=True,
        include_metadata=False,
        include_rt=include_rt_flag,
        rt_column=rt_column_for_prediction,
        rt_output_column=rt_output_column,
        verbose=verbose,
        progress=progress,
    )

    logger.log("[compactlib] building MS2PIP compact library")
    cli_predict_chunked(
        precursors=precursors_for_prediction,
        output=ms2pip_out,
        work_dir=work_dir / f"MS2PIP_max{top_n}.work",
        chunk_size=chunk_size,
        resume=resume,
        force=force,
        backend=prediction_backend,
        preset=preset,
        top_n=top_n,
        fragment_types=None,
        fragment_charges=None,
        min_fragment_series=None,
        max_fragment_series=None,
        min_product_mz=None,
        max_product_mz=None,
        model=ms2pip_model if prediction_backend.lower() == "koina" else None,
        server_url=server_url,
        ssl=ssl,
        collision_energy=collision_energy,
        batch_size=batch_size,
        sequence_source=sequence_source,
        koina_sequence_input="peptide_sequences",
        koina_charge_input="precursor_charges",
        koina_ce_input="collision_energies",
        koina_intensity_col="intensities",
        koina_annotation_col="annotation",
        koina_mz_col=None,
        drop_zero_intensity=True,
        drop_unsupported_aa=True,
        include_metadata=False,
        include_rt=include_rt_flag,
        rt_column=rt_column_for_prediction,
        rt_output_column=rt_output_column,
        verbose=verbose,
        progress=progress,
    )

    # ------------------------------------------------------------------
    # 5. Union2N
    # ------------------------------------------------------------------
    union_out = None
    if union:
        union_out = output_dir / f"{prefix}.Union{2 * top_n}.tsv"
        if resume and union_out.exists() and not force:
            logger.log(f"[compactlib] [skip] union: {union_out} exists")
        else:
            with logger.step(f"reading Prosit and MS2PIP libraries for Union{2 * top_n}"):
                df_prosit = read_library(prosit_out)
                df_ms2pip = read_library(ms2pip_out)
            with logger.step(f"building Union{2 * top_n}"):
                df_union, union_stats = build_union(df_prosit, df_ms2pip, top_n=top_n)
            with logger.step(f"writing Union{2 * top_n} library to {union_out}"):
                write_library(df_union, union_out)
            union_summary_path, union_params_path = sidecar_paths(union_out)
            write_summary(library_summary(None, df_union, union_out, mode="union", extra=union_stats), union_summary_path)
            write_params({
                "command": "build:union",
                "input_a": str(prosit_out),
                "input_b": str(ms2pip_out),
                "output": str(union_out),
                "top_n": top_n,
            }, union_params_path)

    # ------------------------------------------------------------------
    # 6. Build-level manifest
    # ------------------------------------------------------------------
    build_manifest = {
        "command": "build",
        "fasta": str(fasta),
        "output_dir": str(output_dir),
        "prefix": prefix,
        "preset": preset,
        "top_n": top_n,
        "union_n": 2 * top_n if union else None,
        "sample_size": sample_size,
        "predict_rt": predict_rt,
        "rt_backend": rt_backend if predict_rt else None,
        "rt_output_column": rt_output_column if predict_rt else None,
        "prediction_backend": prediction_backend,
        "prosit_model": prosit_model,
        "ms2pip_model": ms2pip_model,
        "collision_energy": collision_energy,
        "chunk_size": chunk_size,
        "batch_size": batch_size,
        "resume": resume,
        "force": force,
        "precursors_all": str(precursors_all),
        "precursors_active": str(active_precursors),
        "precursors_for_prediction": str(precursors_for_prediction),
        "prosit_library": str(prosit_out),
        "ms2pip_library": str(ms2pip_out),
        "union_library": str(union_out) if union_out is not None else None,
        "compactlib_elapsed_sec": logger.total_elapsed(),
    }
    write_params(build_manifest, output_dir / f"{prefix}.build.params.json")
    write_summary(pd.DataFrame([build_manifest]), output_dir / f"{prefix}.build.summary.tsv")

    logger.log("[compactlib] build workflow complete")
    logger.log(f"[compactlib] Prosit: {prosit_out}")
    logger.log(f"[compactlib] MS2PIP: {ms2pip_out}")
    if union_out is not None:
        logger.log(f"[compactlib] Union{2 * top_n}: {union_out}")


@app.command("predict-rt")
def cli_predict_rt(
    precursors: Path = typer.Option(..., "--precursors", "-p", help="Input precursor table produced by compactlib digest"),
    output: Path = typer.Option(..., "--output", "-o", help="Output precursor table with predicted RT column"),
    backend: str = typer.Option("deeplc", "--backend", help="RT backend: deeplc or mock"),
    rt_output_column: str = typer.Option("NormalizedRetentionTime", "--rt-output-column", help="RT column name written to output precursor table"),
    modified_column: str = typer.Option("ModifiedPeptide", "--modified-column", help="Modified peptide column used for DeepLC input conversion"),
    stripped_column: str = typer.Option("StrippedPeptide", "--stripped-column", help="Stripped peptide column used for DeepLC seq column"),
    fixed_carbamidomethyl_c: bool = typer.Option(True, "--fixed-carbamidomethyl-c/--no-fixed-carbamidomethyl-c", help="Indicate fixed C carbamidomethylation in DeepLC modifications input"),
    calibration_table: Optional[Path] = typer.Option(None, "--calibration-table", help="Optional calibration table for DeepLC. Can be DeepLC-style seq/modifications/tr or compactlib-like peptide table with RT."),
    calibration_rt_column: str = typer.Option("tr", "--calibration-rt-column", help="RT column in calibration table"),
    deeplc_model_path: Optional[str] = typer.Option(None, "--deeplc-model-path", help="Optional DeepLC model path. Omit to use DeepLC defaults."),
    write_deeplc_input: Optional[Path] = typer.Option(None, "--write-deeplc-input", help="Optional path to write the generated DeepLC input table for audit/debugging"),
    verbose: bool = typer.Option(True, "--verbose/--quiet", help="Print progress logs."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a terminal progress spinner when possible."),
) -> None:
    """Predict RT for a precursor table using DeepLC and write a precursor table with RT.

    The command converts compactlib modified peptide notation into DeepLC's
    ``seq``/``modifications`` input format. By default, fixed carbamidomethyl C
    is explicitly indicated because DeepLC expects fixed modifications to be
    present in the input.
    """
    logger = ProgressLogger(verbose=verbose, progress=progress)

    with logger.step(f"reading precursor table from {precursors}"):
        df_prec = read_table(precursors)
    logger.log(f"[compactlib] precursor table: {format_shape(df_prec)}")

    df_cal = None
    if calibration_table is not None:
        with logger.step(f"reading DeepLC calibration table from {calibration_table}"):
            df_cal = read_table(calibration_table)
        logger.log(f"[compactlib] calibration table: {format_shape(df_cal)}")

    with logger.step(f"predicting RT using backend={backend}"):
        result = add_deeplc_rt_to_precursors(
            df_prec,
            backend=backend,
            modified_column=modified_column,
            stripped_column=stripped_column,
            rt_output_column=rt_output_column,
            fixed_carbamidomethyl_c=fixed_carbamidomethyl_c,
            calibration_table=df_cal,
            calibration_rt_column=calibration_rt_column,
            deeplc_model_path=deeplc_model_path,
        )
    logger.log(f"[compactlib] precursor table with RT: {format_shape(result.precursor_table)}")

    if write_deeplc_input is not None:
        with logger.step(f"writing generated DeepLC input to {write_deeplc_input}"):
            write_table(result.deeplc_input, write_deeplc_input)

    with logger.step(f"writing precursor table with RT to {output}"):
        write_table(result.precursor_table, output)

    summary_path, params_path = sidecar_paths(output)
    stats = {**result.stats, "compactlib_elapsed_sec": logger.total_elapsed()}
    with logger.step(f"writing summary to {summary_path}"):
        write_summary(pd.DataFrame([stats]), summary_path)

    params = {
        "command": "predict-rt",
        "precursors": str(precursors),
        "output": str(output),
        "backend": backend,
        "rt_output_column": rt_output_column,
        "modified_column": modified_column,
        "stripped_column": stripped_column,
        "fixed_carbamidomethyl_c": fixed_carbamidomethyl_c,
        "calibration_table": str(calibration_table) if calibration_table else None,
        "calibration_rt_column": calibration_rt_column,
        "deeplc_model_path": deeplc_model_path,
        "write_deeplc_input": str(write_deeplc_input) if write_deeplc_input else None,
    }
    with logger.step(f"writing parameters to {params_path}"):
        write_params(params, params_path)

    logger.log(f"[compactlib] wrote precursor table with RT: {output}")
    logger.log(f"[compactlib] wrote summary: {summary_path}")
    logger.log(f"[compactlib] wrote params: {params_path}")

@app.command("add-rt")
def cli_add_rt(
    library: Path = typer.Option(..., "--library", "-l", help="Input transition library to annotate with RT"),
    rt_table: Path = typer.Option(..., "--rt-table", help="Table containing precursor-level RT predictions"),
    output: Path = typer.Option(..., "--output", "-o", help="Output library with RT column"),
    key: str = typer.Option("ModifiedPeptide,PrecursorCharge", "--key", help="Comma-separated merge key columns, e.g. ModifiedPeptide,PrecursorCharge"),
    rt_column: str = typer.Option(..., "--rt-column", help="RT column name in --rt-table, e.g. Tr_recalibrated"),
    rt_output_column: str = typer.Option("NormalizedRetentionTime", "--rt-output-column", help="RT column name to write to output library"),
    min_match_rate: float = typer.Option(0.99, "--min-match-rate", min=0.0, max=1.0, help="Minimum required precursor-level RT match rate"),
    verbose: bool = typer.Option(True, "--verbose/--quiet", help="Print progress logs."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a terminal progress spinner when possible."),
) -> None:
    """Add externally predicted RT values to an existing transition library.

    This command intentionally does not predict RT itself. It merges an RT table
    produced by an external tool, such as DeepLC, into a compact transition
    library using a precursor-level key.
    """
    logger = ProgressLogger(verbose=verbose, progress=progress)
    key_cols = parse_key_columns(key)

    with logger.step(f"reading transition library from {library}"):
        df_lib = read_library(library)
    logger.log(f"[compactlib] library: {format_shape(df_lib)}")

    with logger.step(f"reading RT table from {rt_table}"):
        df_rt = read_table(rt_table)
    logger.log(f"[compactlib] RT table: {format_shape(df_rt)}")

    with logger.step(f"merging RT column '{rt_column}' using key={','.join(key_cols)}"):
        df_out, rt_stats = attach_rt(
            df_lib,
            df_rt,
            key_cols=key_cols,
            rt_column=rt_column,
            rt_output_column=rt_output_column,
            min_match_rate=min_match_rate,
        )
    logger.log(
        "[compactlib] RT precursor match rate: "
        f"{rt_stats['n_matched_precursors']:,}/{rt_stats['n_library_precursors']:,} "
        f"({rt_stats['pct_matched_precursors']:.3f}%)"
    )

    params = {
        "command": "add-rt",
        "library": str(library),
        "rt_table": str(rt_table),
        "output": str(output),
        "key": key_cols,
        "rt_column": rt_column,
        "rt_output_column": rt_output_column,
        "min_match_rate": min_match_rate,
    }
    _finish(df_lib, df_out, output, mode="add-rt", params=params, extra_summary=rt_stats, logger=logger)


@app.command("max-n")
def cli_max_n(
    input: Path = typer.Option(..., "--input", "-i", help="Input library: .tsv, .csv, .csv.gz, .parquet"),
    output: Path = typer.Option(..., "--output", "-o", help="Output compact library"),
    top_n: int = typer.Option(..., "--top-n", "-n", min=1, help="Maximum transitions per precursor"),
    verbose: bool = typer.Option(True, "--verbose/--quiet", help="Print progress logs."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a terminal progress spinner when possible."),
) -> None:
    """Keep up to N most intense transitions per precursor."""
    logger = ProgressLogger(verbose=verbose, progress=progress)
    df = _read_one(input, logger)
    with logger.step(f"selecting top-{top_n} transitions per precursor"):
        out, stats = select_top_n(df, top_n=top_n, deduplicate=True)
    logger.log(f"[compactlib] selected library: {format_shape(out)}")
    params = {"command": "max-n", "input": str(input), "output": str(output), "top_n": top_n}
    _finish(df, out, output, mode="max-n", params=params, extra_summary=stats, logger=logger)


@app.command("union")
def cli_union(
    input_a: Path = typer.Option(..., "--input-a", help="Input library A"),
    input_b: Path = typer.Option(..., "--input-b", help="Input library B"),
    output: Path = typer.Option(..., "--output", "-o", help="Output union library"),
    top_n: int = typer.Option(..., "--top-n", "-n", min=1, help="Top-N transitions from each input library"),
    verbose: bool = typer.Option(True, "--verbose/--quiet", help="Print progress logs."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a terminal progress spinner when possible."),
) -> None:
    """Build Union2N library: top-N(A) union top-N(B), with duplicate transitions removed."""
    logger = ProgressLogger(verbose=verbose, progress=progress)
    df_a = _read_one(input_a, logger, label="A")
    df_b = _read_one(input_b, logger, label="B")
    with logger.step(f"building Union{2 * top_n} library from top-{top_n} transitions of A and B"):
        out, stats = build_union(df_a, df_b, top_n=top_n)
    logger.log(f"[compactlib] union library: {format_shape(out)}")
    params = {
        "command": "union",
        "input_a": str(input_a),
        "input_b": str(input_b),
        "output": str(output),
        "top_n": top_n,
        "transition_key": "ModifiedPeptide+PrecursorCharge+FragmentType+FragmentSeriesNumber+FragmentCharge",
    }
    _finish(None, out, output, mode="union", params=params, extra_summary=stats, logger=logger)


@app.command("consensus")
def cli_consensus(
    input_a: Path = typer.Option(..., "--input-a", help="Input library A"),
    input_b: Path = typer.Option(..., "--input-b", help="Input library B"),
    output: Path = typer.Option(..., "--output", "-o", help="Output consensus library"),
    top_n: int = typer.Option(..., "--top-n", "-n", min=1, help="Consensus top-N transitions per precursor"),
    verbose: bool = typer.Option(True, "--verbose/--quiet", help="Print progress logs."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a terminal progress spinner when possible."),
) -> None:
    """Build consensus top-N using rank-based score 1/rank_a + 1/rank_b."""
    logger = ProgressLogger(verbose=verbose, progress=progress)
    df_a = _read_one(input_a, logger, label="A")
    df_b = _read_one(input_b, logger, label="B")
    with logger.step(f"building consensus top-{top_n} library using rank-based score"):
        out, stats = build_consensus(df_a, df_b, top_n=top_n)
    logger.log(f"[compactlib] consensus library: {format_shape(out)}")
    params = {
        "command": "consensus",
        "input_a": str(input_a),
        "input_b": str(input_b),
        "output": str(output),
        "top_n": top_n,
        "score": "1/rank_a + 1/rank_b",
        "transition_key": "ModifiedPeptide+PrecursorCharge+FragmentType+FragmentSeriesNumber+FragmentCharge",
    }
    _finish(None, out, output, mode="consensus", params=params, extra_summary=stats, logger=logger)


@app.command("random-n")
def cli_random_n(
    input: Path = typer.Option(..., "--input", "-i", help="Input library"),
    output: Path = typer.Option(..., "--output", "-o", help="Output random-N library"),
    top_n: int = typer.Option(..., "--top-n", "-n", min=1, help="Maximum transitions per precursor"),
    seed: int = typer.Option(42, "--seed", help="Random seed"),
    verbose: bool = typer.Option(True, "--verbose/--quiet", help="Print progress logs."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a terminal progress spinner when possible."),
) -> None:
    """Randomly keep up to N transitions per precursor using a fixed seed."""
    logger = ProgressLogger(verbose=verbose, progress=progress)
    df = _read_one(input, logger)
    with logger.step(f"randomly selecting up to {top_n} transitions per precursor with seed={seed}"):
        out, stats = select_random_n(df, top_n=top_n, seed=seed, deduplicate=True)
    logger.log(f"[compactlib] random library: {format_shape(out)}")
    params = {"command": "random-n", "input": str(input), "output": str(output), "top_n": top_n, "seed": seed}
    _finish(df, out, output, mode="random-n", params=params, extra_summary=stats, logger=logger)


@app.command("reverse-intensity")
def cli_reverse_intensity(
    input: Path = typer.Option(..., "--input", "-i", help="Input library"),
    output: Path = typer.Option(..., "--output", "-o", help="Output reverse-intensity library"),
    verbose: bool = typer.Option(True, "--verbose/--quiet", help="Print progress logs."),
    progress: bool = typer.Option(True, "--progress/--no-progress", help="Show a terminal progress spinner when possible."),
) -> None:
    """Invert intensity pattern within each precursor while keeping the same transitions."""
    logger = ProgressLogger(verbose=verbose, progress=progress)
    df = _read_one(input, logger)
    with logger.step("reversing intensity profiles within each precursor"):
        out, stats = reverse_intensity(df)
    logger.log(f"[compactlib] reverse-intensity library: {format_shape(out)}")
    params = {"command": "reverse-intensity", "input": str(input), "output": str(output)}
    _finish(df, out, output, mode="reverse-intensity", params=params, extra_summary=stats, logger=logger)


if __name__ == "__main__":
    app()
