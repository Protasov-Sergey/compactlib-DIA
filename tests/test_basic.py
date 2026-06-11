from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
from typer.testing import CliRunner

from compactlib.cli import app
from compactlib.io import PRECURSOR_KEY

runner = CliRunner()


def toy_library(scale: float = 1.0) -> pd.DataFrame:
    rows = []
    for pep, charge, pmz in [("PEPTIDE", 2, 500.2), ("ACDEFGK", 2, 600.3)]:
        for i, ion in enumerate(["b", "y", "b", "y", "b"], start=1):
            rows.append({
                "ModifiedPeptide": pep,
                "PrecursorCharge": charge,
                "PrecursorMz": pmz,
                "ProductMz": 100.0 + i + (0.1 if ion == "y" else 0.0),
                "LibraryIntensity": scale * (100 - i * 10),
                "FragmentType": ion,
                "FragmentSeriesNumber": i,
                "FragmentCharge": 1,
                "ExtraColumn": f"extra_{pep}_{i}",
            })
    return pd.DataFrame(rows)


def write_tsv(df: pd.DataFrame, path):
    df.to_csv(path, sep="\t", index=False)


def test_max_n_cli(tmp_path):
    inp = tmp_path / "lib.tsv"
    out = tmp_path / "max3.tsv"
    write_tsv(toy_library(), inp)

    res = runner.invoke(app, ["max-n", "--input", str(inp), "--output", str(out), "--top-n", "3"])
    assert res.exit_code == 0, res.output
    df = pd.read_csv(out, sep="\t")
    assert df.groupby(PRECURSOR_KEY).size().max() <= 3
    assert "ExtraColumn" in df.columns
    assert (tmp_path / "max3.summary.tsv").exists()
    assert (tmp_path / "max3.params.json").exists()


def test_union_cli(tmp_path):
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    out = tmp_path / "union.tsv"
    write_tsv(toy_library(scale=1.0), a)
    write_tsv(toy_library(scale=0.8), b)

    res = runner.invoke(app, ["union", "--input-a", str(a), "--input-b", str(b), "--output", str(out), "--top-n", "2"])
    assert res.exit_code == 0, res.output
    df = pd.read_csv(out, sep="\t")
    # Since a and b have identical transition annotations, union after dedup has at most 2 per precursor.
    assert df.groupby(PRECURSOR_KEY).size().max() <= 2


def test_consensus_cli(tmp_path):
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    out = tmp_path / "consensus.tsv"
    write_tsv(toy_library(scale=1.0), a)
    write_tsv(toy_library(scale=0.8), b)

    res = runner.invoke(app, ["consensus", "--input-a", str(a), "--input-b", str(b), "--output", str(out), "--top-n", "3"])
    assert res.exit_code == 0, res.output
    df = pd.read_csv(out, sep="\t")
    assert df.groupby(PRECURSOR_KEY).size().max() <= 3


def test_random_reproducible(tmp_path):
    inp = tmp_path / "lib.tsv"
    out1 = tmp_path / "random1.tsv"
    out2 = tmp_path / "random2.tsv"
    write_tsv(toy_library(), inp)

    args1 = ["random-n", "--input", str(inp), "--output", str(out1), "--top-n", "2", "--seed", "123"]
    args2 = ["random-n", "--input", str(inp), "--output", str(out2), "--top-n", "2", "--seed", "123"]
    assert runner.invoke(app, args1).exit_code == 0
    assert runner.invoke(app, args2).exit_code == 0

    df1 = pd.read_csv(out1, sep="\t")
    df2 = pd.read_csv(out2, sep="\t")
    pd.testing.assert_frame_equal(df1.reset_index(drop=True), df2.reset_index(drop=True))


def test_reverse_intensity_preserves_transitions(tmp_path):
    inp = tmp_path / "lib.tsv"
    out = tmp_path / "rev.tsv"
    df0 = toy_library()
    write_tsv(df0, inp)

    res = runner.invoke(app, ["reverse-intensity", "--input", str(inp), "--output", str(out)])
    assert res.exit_code == 0, res.output
    df = pd.read_csv(out, sep="\t")
    assert len(df) == len(df0)
    assert set(map(tuple, df[["ModifiedPeptide", "PrecursorCharge", "FragmentType", "FragmentSeriesNumber", "FragmentCharge"]].values)) == set(
        map(tuple, df0[["ModifiedPeptide", "PrecursorCharge", "FragmentType", "FragmentSeriesNumber", "FragmentCharge"]].values)
    )


def test_two_input_summaries_have_input_a_b_fields_and_no_nan_input_fields(tmp_path):
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    out_union = tmp_path / "union.tsv"
    out_cons = tmp_path / "consensus.tsv"
    write_tsv(toy_library(scale=1.0), a)
    write_tsv(toy_library(scale=0.8), b)

    res_u = runner.invoke(app, ["union", "--input-a", str(a), "--input-b", str(b), "--output", str(out_union), "--top-n", "2"])
    assert res_u.exit_code == 0, res_u.output
    sum_u = pd.read_csv(tmp_path / "union.summary.tsv", sep="\t")

    res_c = runner.invoke(app, ["consensus", "--input-a", str(a), "--input-b", str(b), "--output", str(out_cons), "--top-n", "2"])
    assert res_c.exit_code == 0, res_c.output
    sum_c = pd.read_csv(tmp_path / "consensus.summary.tsv", sep="\t")

    for s in [sum_u, sum_c]:
        assert "n_precursors_input" not in s.columns
        assert "n_transitions_input" not in s.columns
        for col in [
            "n_precursors_input_a", "n_transitions_input_a",
            "n_precursors_input_b", "n_transitions_input_b",
            "n_common_precursors", "n_a_only_precursors", "n_b_only_precursors",
            "duplicate_policy",
        ]:
            assert col in s.columns
        assert not s[["n_precursors_input_a", "n_transitions_input_a", "n_precursors_input_b", "n_transitions_input_b"]].isna().any().any()



def test_digest_builds_precursor_table(tmp_path):
    from compactlib.digest import build_precursor_table_from_fasta

    fasta = tmp_path / "toy.fasta"
    fasta.write_text(
        ">sp|P1|P1_HUMAN toy\nMPEPTIDEKACDEK\n"
        ">sp|P2|P2_HUMAN toy\nACDEKPEPTIDER\n",
        encoding="utf-8",
    )

    df, stats = build_precursor_table_from_fasta(
        fasta,
        enzyme="trypsin-p",
        missed_cleavages=1,
        min_length=5,
        max_length=30,
        charges="2,3",
    )

    assert {"ModifiedPeptide", "StrippedPeptide", "PrecursorCharge", "PrecursorMz", "ProteinId"}.issubset(df.columns)
    assert set(df["PrecursorCharge"]) == {2, 3}
    assert stats["n_proteins_input"] == 2
    assert stats["n_unique_stripped_peptides"] >= 2
    assert len(df) == stats["n_precursors_output"]
    assert df["PrecursorMz"].notna().all()


def test_digest_variable_oxidation_m(tmp_path):
    from compactlib.digest import build_precursor_table_from_fasta

    fasta = tmp_path / "toy_mox.fasta"
    fasta.write_text(
        ">sp|PMOX|PMOX_HUMAN toy\nMPEPTIDEMK\n",
        encoding="utf-8",
    )

    df, stats = build_precursor_table_from_fasta(
        fasta,
        enzyme="trypsin-p",
        missed_cleavages=0,
        min_length=5,
        max_length=30,
        charges="2",
        variable_oxidation_m=True,
        max_variable_mods=1,
    )

    # Peptide has two M residues, therefore unmodified + two singly oxidized peptidoforms.
    assert stats["n_unique_stripped_peptides"] == 1
    assert stats["n_unique_modified_peptides"] == 3
    assert stats["n_oxidized_modified_peptides"] == 2
    assert len(df) == 3
    assert df["ModifiedPeptide"].str.contains("M\\(UniMod:35\\)", regex=True).sum() == 2
    assert set(df["NVariableMods"]) == {0, 1}
    assert df["PrecursorMz"].round(6).nunique() == 2  # both singly oxidized forms have the same precursor mass


def test_predict_mock_backend_and_topn(tmp_path):
    fasta = tmp_path / "toy.fasta"
    prec = tmp_path / "prec.tsv"
    out = tmp_path / "pred_max3.tsv"
    fasta.write_text(
        ">sp|P1|P1_HUMAN toy\nMPEPTIDEK\n",
        encoding="utf-8",
    )

    res_d = runner.invoke(app, [
        "digest", "--fasta", str(fasta), "--output", str(prec),
        "--enzyme", "trypsin-p", "--missed-cleavages", "0",
        "--min-length", "5", "--max-length", "30", "--charges", "2",
        "--variable-oxidation-m", "--max-variable-mods", "1",
    ])
    assert res_d.exit_code == 0, res_d.output

    res_p = runner.invoke(app, [
        "predict", "--precursors", str(prec), "--output", str(out),
        "--backend", "mock", "--fragment-types", "b,y", "--fragment-charges", "1",
        "--top-n", "3",
    ])
    assert res_p.exit_code == 0, res_p.output
    df = pd.read_csv(out, sep="\t")
    expected_minimal = [
        "ModifiedPeptide", "StrippedPeptide", "PrecursorCharge", "PrecursorMz",
        "ProductMz", "LibraryIntensity", "FragmentType", "FragmentSeriesNumber",
        "FragmentCharge", "ProteinId",
    ]
    assert list(df.columns) == expected_minimal
    assert df.groupby(PRECURSOR_KEY).size().max() <= 3
    assert "PredictionBackend" not in df.columns
    assert (tmp_path / "pred_max3.summary.tsv").exists()
    assert (tmp_path / "pred_max3.params.json").exists()

    out_meta = tmp_path / "pred_max3_meta.tsv"
    res_p_meta = runner.invoke(app, [
        "predict", "--precursors", str(prec), "--output", str(out_meta),
        "--backend", "mock", "--fragment-types", "b,y", "--fragment-charges", "1",
        "--top-n", "3", "--include-metadata",
    ])
    assert res_p_meta.exit_code == 0, res_p_meta.output
    df_meta = pd.read_csv(out_meta, sep="\t")
    assert "PredictionBackend" in df_meta.columns
    assert "PeptideLength" in df_meta.columns


def test_digest_precursor_mz_filter_and_charge_length_rules(tmp_path):
    from compactlib.digest import build_precursor_table_from_fasta

    fasta = tmp_path / "toy_filters.fasta"
    fasta.write_text(
        ">sp|P1|P1_HUMAN toy\nPEPTIDEKVERYVERYLONGPEPTIDER\n",
        encoding="utf-8",
    )

    df_no_filter, stats_no_filter = build_precursor_table_from_fasta(
        fasta,
        enzyme="trypsin-p",
        missed_cleavages=0,
        min_length=7,
        max_length=40,
        charges="2,3,4",
    )
    assert 4 in set(df_no_filter["PrecursorCharge"])

    df, stats = build_precursor_table_from_fasta(
        fasta,
        enzyme="trypsin-p",
        missed_cleavages=0,
        min_length=7,
        max_length=40,
        charges="2,3,4",
        min_precursor_mz=400,
        max_precursor_mz=1200,
        charge_length_rules="4:12",
    )

    assert stats["n_precursors_before_precursor_filters"] >= len(df)
    assert stats["n_precursors_filtered_mz"] >= 0
    assert stats["n_precursors_filtered_charge_length"] >= 0
    assert df["PrecursorMz"].between(400, 1200).all()
    assert not ((df["PrecursorCharge"] == 4) & (df["PeptideLength"] < 12)).any()


def test_digest_cli_precursor_filters(tmp_path):
    fasta = tmp_path / "toy_cli_filters.fasta"
    out = tmp_path / "precursors.tsv"
    fasta.write_text(
        ">sp|P1|P1_HUMAN toy\nPEPTIDEKVERYVERYLONGPEPTIDER\n",
        encoding="utf-8",
    )

    res = runner.invoke(app, [
        "digest",
        "--fasta", str(fasta),
        "--output", str(out),
        "--enzyme", "trypsin-p",
        "--missed-cleavages", "0",
        "--min-length", "7",
        "--max-length", "40",
        "--charges", "2,3,4",
        "--min-precursor-mz", "400",
        "--max-precursor-mz", "1200",
        "--charge-length-rules", "4:12",
    ])
    assert res.exit_code == 0, res.output
    df = pd.read_csv(out, sep="\t")
    assert df["PrecursorMz"].between(400, 1200).all()
    assert not ((df["PrecursorCharge"] == 4) & (df["PeptideLength"] < 12)).any()
    summary = pd.read_csv(tmp_path / "precursors.summary.tsv", sep="\t")
    assert "n_precursors_filtered_mz" in summary.columns
    assert "n_precursors_filtered_charge_length" in summary.columns


def test_presets_command_lists_presets():
    res = runner.invoke(app, ["presets"])
    assert res.exit_code == 0, res.output
    assert "diann-like" in res.output
    assert "skyline-dia-like" in res.output
    assert "generic-dia" in res.output
    assert "project-human-dia" not in res.output


def test_digest_project_preset_records_params_and_filters(tmp_path):
    fasta = tmp_path / "toy_preset.fasta"
    out = tmp_path / "prec.tsv"
    fasta.write_text(
        ">sp|P1|P1_HUMAN toy\nPEPTIDEKVERYVERYLONGPEPTIDER\n",
        encoding="utf-8",
    )

    res = runner.invoke(app, [
        "digest", "--fasta", str(fasta), "--output", str(out),
        "--preset", "generic-dia-strict",
        "--max-precursor-mz", "1200",  # user override should replace preset 1000
    ])
    assert res.exit_code == 0, res.output
    df = pd.read_csv(out, sep="\t")
    summary = pd.read_csv(tmp_path / "prec.summary.tsv", sep="\t")
    assert summary.loc[0, "preset"] == "generic-dia-strict"
    assert float(summary.loc[0, "max_precursor_mz"]) == 1200.0
    assert summary.loc[0, "charge_length_rules"] == "4:12"
    assert "PrecursorMz" in df.columns
    assert len(df) == int(summary.loc[0, "n_precursors_output"])


def test_predict_preset_filters_mock_transitions(tmp_path):
    fasta = tmp_path / "toy_predict_preset.fasta"
    prec = tmp_path / "prec.tsv"
    out = tmp_path / "pred.tsv"
    fasta.write_text(
        ">sp|P1|P1_HUMAN toy\nMPEPTIDEK\n",
        encoding="utf-8",
    )
    res_d = runner.invoke(app, [
        "digest", "--fasta", str(fasta), "--output", str(prec),
        "--min-length", "5", "--max-length", "30", "--charges", "2",
        "--missed-cleavages", "0",
    ])
    assert res_d.exit_code == 0, res_d.output

    res_p = runner.invoke(app, [
        "predict", "--precursors", str(prec), "--output", str(out),
        "--backend", "mock", "--preset", "diann-like", "--top-n", "10",
    ])
    assert res_p.exit_code == 0, res_p.output
    df = pd.read_csv(out, sep="\t")
    assert set(df["FragmentType"].astype(str).str.lower()).issubset({"b", "y"})
    assert set(df["FragmentCharge"]).issubset({1, 2})
    assert df["FragmentSeriesNumber"].min() >= 3
    summary = pd.read_csv(tmp_path / "pred.summary.tsv", sep="\t")
    assert summary.loc[0, "preset"] == "diann-like"
    assert summary.loc[0, "transition_filter_min_fragment_series"] == 3


def test_skyline_audited_preset_matches_export_audit():
    from compactlib.presets import resolve_digest_preset

    cfg = resolve_digest_preset("skyline-dia-like")
    assert cfg["min_length"] == 7
    assert cfg["max_length"] == 30
    assert cfg["charges"] == "1,2,3,4"
    assert cfg["min_precursor_mz"] == 300.0
    assert cfg["max_precursor_mz"] == 1800.0
    assert cfg["charge_length_rules"] == ""

    generic_cfg = resolve_digest_preset("generic-dia")
    assert generic_cfg["charges"] == "1,2,3,4"
    assert generic_cfg["min_precursor_mz"] == 300.0
    assert generic_cfg["max_precursor_mz"] == 1800.0

    try:
        resolve_digest_preset("project-human-dia")
    except ValueError as exc:
        assert "Unknown preset" in str(exc)
    else:
        raise AssertionError("project-human-dia alias should have been removed")


def test_digest_selenocysteine_optional(tmp_path):
    from compactlib.digest import build_precursor_table_from_fasta

    fasta = tmp_path / "toy_u.fasta"
    fasta.write_text(
        ">sp|PU|PU_HUMAN toy\nAEENITESCQUR\n",
        encoding="utf-8",
    )

    df_default, stats_default = build_precursor_table_from_fasta(
        fasta,
        enzyme="trypsin-p",
        missed_cleavages=0,
        min_length=5,
        max_length=30,
        charges="2",
    )
    assert len(df_default) == 0
    assert stats_default["n_invalid_aa_filtered_peptides"] == 1

    df_u, stats_u = build_precursor_table_from_fasta(
        fasta,
        enzyme="trypsin-p",
        missed_cleavages=0,
        min_length=5,
        max_length=30,
        charges="2",
        allow_selenocysteine=True,
    )
    assert len(df_u) == 1
    assert df_u["ContainsSelenocysteine"].iloc[0] is True or bool(df_u["ContainsSelenocysteine"].iloc[0]) is True
    assert df_u["PrecursorMz"].notna().all()
    assert stats_u["n_unique_selenocysteine_peptides"] == 1
    assert stats_u["n_selenocysteine_precursors_output"] == 1


def test_koina_long_format_adapter_without_koinapy():
    from compactlib.predictors.koina import KoinaPredictor

    # Bypass __post_init__ because this unit test only exercises the converter
    # and should not require koinapy/network access.
    predictor = KoinaPredictor.__new__(KoinaPredictor)
    predictor.model_name = "Prosit_2020_intensity_HCD"
    predictor.name = "koina"
    predictor.sequence_source = "stripped"
    predictor.sequence_input = "peptide_sequences"
    predictor.charge_input = "precursor_charges"
    predictor.ce_input = "collision_energies"
    predictor.intensity_col = "intensities"
    predictor.annotation_col = "annotation"
    predictor.mz_col = None
    predictor.drop_zero_intensity = True
    predictor.collision_energy = 30.0

    batch = pd.DataFrame({
        "ModifiedPeptide": ["PEPTIDEK", "ACDEK"],
        "StrippedPeptide": ["PEPTIDEK", "ACDEK"],
        "PrecursorCharge": [2, 2],
        "PrecursorMz": [500.0, 400.0],
        "FixedCarbamidomethylC": [True, True],
        "OxidationMPositions": ["", ""],
    })
    pred = pd.DataFrame({
        "peptide_sequences": ["PEPTIDEK", "PEPTIDEK", "ACDEK"],
        "precursor_charges": [2, 2, 2],
        "collision_energies": [30, 30, 30],
        "annotation": ["y3", "b2", "y2"],
        "intensities": [0.9, 0.4, 0.8],
    })

    out = predictor._convert_predictions(batch, pred)
    assert len(out) == 3
    assert {"ProductMz", "LibraryIntensity", "FragmentType", "FragmentSeriesNumber", "FragmentCharge"}.issubset(out.columns)
    assert set(out["ModifiedPeptide"]) == {"PEPTIDEK", "ACDEK"}
    assert set(out["FragmentType"]) == {"b", "y"}
    assert out["ProductMz"].notna().all()


def test_predict_chunked_mock_backend_resume(tmp_path):
    fasta = tmp_path / "toy_chunked.fasta"
    prec = tmp_path / "prec.tsv"
    out = tmp_path / "chunked_pred.tsv"
    work = tmp_path / "chunked_work"
    fasta.write_text(
        ">sp|P1|P1_HUMAN toy\nMPEPTIDEKACDEK\n"
        ">sp|P2|P2_HUMAN toy\nACDEKPEPTIDER\n",
        encoding="utf-8",
    )

    res_d = runner.invoke(app, [
        "digest", "--fasta", str(fasta), "--output", str(prec),
        "--enzyme", "trypsin-p", "--missed-cleavages", "1",
        "--min-length", "5", "--max-length", "30", "--charges", "2",
    ])
    assert res_d.exit_code == 0, res_d.output

    args = [
        "predict-chunked", "--precursors", str(prec), "--output", str(out),
        "--work-dir", str(work), "--chunk-size", "2",
        "--backend", "mock", "--fragment-types", "b,y", "--fragment-charges", "1",
        "--top-n", "3", "--quiet", "--no-progress",
    ]
    res_p = runner.invoke(app, args)
    assert res_p.exit_code == 0, res_p.output
    assert out.exists()
    assert (tmp_path / "chunked_pred.summary.tsv").exists()
    assert (tmp_path / "chunked_pred.params.json").exists()
    assert (work / "chunks_manifest.tsv").exists()

    df = pd.read_csv(out, sep="\t")
    assert len(df) > 0
    assert df.groupby(PRECURSOR_KEY).size().max() <= 3

    manifest = pd.read_csv(work / "chunks_manifest.tsv", sep="\t")
    assert manifest["completed"].all()
    assert len(list((work / "predicted_chunks").glob("*.done.json"))) == len(manifest)

    # Re-running the same command should resume/skip completed chunks and still succeed.
    res_p2 = runner.invoke(app, args)
    assert res_p2.exit_code == 0, res_p2.output


def test_predict_chunked_can_remerge_full_metadata_chunks_as_minimal(tmp_path):
    fasta = tmp_path / "toy_chunked_meta.fasta"
    prec = tmp_path / "prec_meta.tsv"
    out_full = tmp_path / "chunked_full.tsv"
    out_min = tmp_path / "chunked_min.tsv"
    work = tmp_path / "chunked_meta_work"
    fasta.write_text(
        ">sp|P1|P1_HUMAN toy\nMPEPTIDEKACDEK\n",
        encoding="utf-8",
    )

    res_d = runner.invoke(app, [
        "digest", "--fasta", str(fasta), "--output", str(prec),
        "--enzyme", "trypsin-p", "--missed-cleavages", "1",
        "--min-length", "5", "--max-length", "30", "--charges", "2",
        "--quiet", "--no-progress",
    ])
    assert res_d.exit_code == 0, res_d.output

    args_full = [
        "predict-chunked", "--precursors", str(prec), "--output", str(out_full),
        "--work-dir", str(work), "--chunk-size", "2",
        "--backend", "mock", "--fragment-types", "b,y", "--fragment-charges", "1",
        "--top-n", "3", "--include-metadata", "--quiet", "--no-progress",
    ]
    res_full = runner.invoke(app, args_full)
    assert res_full.exit_code == 0, res_full.output
    df_full = pd.read_csv(out_full, sep="\t")
    assert "PredictionBackend" in df_full.columns

    # Same completed full-metadata chunks should be re-merged into a minimal final output
    # without recomputing chunks.
    args_min = [
        "predict-chunked", "--precursors", str(prec), "--output", str(out_min),
        "--work-dir", str(work), "--chunk-size", "2",
        "--backend", "mock", "--fragment-types", "b,y", "--fragment-charges", "1",
        "--top-n", "3", "--quiet", "--no-progress",
    ]
    res_min = runner.invoke(app, args_min)
    assert res_min.exit_code == 0, res_min.output
    df_min = pd.read_csv(out_min, sep="\t")
    assert list(df_min.columns) == [
        "ModifiedPeptide", "StrippedPeptide", "PrecursorCharge", "PrecursorMz",
        "ProductMz", "LibraryIntensity", "FragmentType", "FragmentSeriesNumber",
        "FragmentCharge", "ProteinId",
    ]
    assert "PredictionBackend" not in df_min.columns
    assert len(df_min) == len(df_full)


def test_add_rt_cli_merges_external_rt_table(tmp_path):
    lib = tmp_path / "lib.tsv"
    rt = tmp_path / "rt.tsv"
    out = tmp_path / "lib_rt.tsv"
    df = toy_library()
    # toy_library lacks StrippedPeptide/ProteinId but add-rt only requires transition columns.
    write_tsv(df, lib)
    rt_df = df[["ModifiedPeptide", "PrecursorCharge"]].drop_duplicates().copy()
    rt_df["Tr_recalibrated"] = [1.23, 4.56]
    write_tsv(rt_df, rt)

    res = runner.invoke(app, [
        "add-rt", "--library", str(lib), "--rt-table", str(rt), "--output", str(out),
        "--key", "ModifiedPeptide,PrecursorCharge", "--rt-column", "Tr_recalibrated",
        "--rt-output-column", "NormalizedRetentionTime", "--min-match-rate", "1.0",
    ])
    assert res.exit_code == 0, res.output
    merged = pd.read_csv(out, sep="\t")
    assert "NormalizedRetentionTime" in merged.columns
    assert merged["NormalizedRetentionTime"].notna().all()
    summary = pd.read_csv(tmp_path / "lib_rt.summary.tsv", sep="\t")
    assert summary.loc[0, "n_matched_precursors"] == 2
    assert summary.loc[0, "n_unmatched_precursors"] == 0
    assert summary.loc[0, "pct_matched_precursors"] == 100.0


def test_predict_include_rt_passthrough_minimal_output(tmp_path):
    fasta = tmp_path / "toy_rt.fasta"
    prec = tmp_path / "prec_rt.tsv"
    out = tmp_path / "pred_rt.tsv"
    fasta.write_text(
        ">sp|P1|P1_HUMAN toy\nMPEPTIDEK\n",
        encoding="utf-8",
    )

    res_d = runner.invoke(app, [
        "digest", "--fasta", str(fasta), "--output", str(prec),
        "--enzyme", "trypsin-p", "--missed-cleavages", "0",
        "--min-length", "5", "--max-length", "30", "--charges", "2",
        "--quiet", "--no-progress",
    ])
    assert res_d.exit_code == 0, res_d.output
    df_prec = pd.read_csv(prec, sep="\t")
    df_prec["Tr_recalibrated"] = 7.89
    df_prec.to_csv(prec, sep="\t", index=False)

    res_p = runner.invoke(app, [
        "predict", "--precursors", str(prec), "--output", str(out),
        "--backend", "mock", "--fragment-types", "b,y", "--fragment-charges", "1",
        "--top-n", "3", "--include-rt", "--rt-column", "Tr_recalibrated",
        "--rt-output-column", "NormalizedRetentionTime", "--quiet", "--no-progress",
    ])
    assert res_p.exit_code == 0, res_p.output
    pred = pd.read_csv(out, sep="\t")
    assert list(pred.columns) == [
        "ModifiedPeptide", "StrippedPeptide", "PrecursorCharge", "PrecursorMz",
        "ProductMz", "LibraryIntensity", "FragmentType", "FragmentSeriesNumber",
        "FragmentCharge", "ProteinId", "NormalizedRetentionTime",
    ]
    assert pred["NormalizedRetentionTime"].eq(7.89).all()
    summary = pd.read_csv(tmp_path / "pred_rt.summary.tsv", sep="\t")
    assert bool(summary.loc[0, "include_rt"]) is True or str(summary.loc[0, "include_rt"]).lower() == "true"
    assert summary.loc[0, "rt_column"] == "Tr_recalibrated"


def test_deeplc_input_conversion_fixed_c_and_oxidation():
    from compactlib.rt_predictors import parse_modified_peptide_for_deeplc

    seq, mods = parse_modified_peptide_for_deeplc(
        "AC(UniMod:4)DM(UniMod:35)K",
        stripped_peptide="ACDMK",
        fixed_carbamidomethyl_c=True,
    )
    assert seq == "ACDMK"
    assert "2|Carbamidomethyl" in mods
    assert "4|Oxidation" in mods

    seq2, mods2 = parse_modified_peptide_for_deeplc(
        "ACDMK",
        stripped_peptide="ACDMK",
        fixed_carbamidomethyl_c=True,
    )
    assert seq2 == "ACDMK"
    assert mods2 == "2|Carbamidomethyl"


def test_predict_rt_mock_backend_and_passthrough_to_predict(tmp_path):
    fasta = tmp_path / "toy_rt_deeplc.fasta"
    prec = tmp_path / "prec.tsv"
    prec_rt = tmp_path / "prec_with_rt.tsv"
    pred = tmp_path / "pred_with_rt.tsv"
    deeplc_input = tmp_path / "deeplc_input.tsv"
    fasta.write_text(
        ">sp|P1|P1_HUMAN toy\nACDMKPEPTIDER\n",
        encoding="utf-8",
    )

    res_d = runner.invoke(app, [
        "digest", "--fasta", str(fasta), "--output", str(prec),
        "--enzyme", "trypsin-p", "--missed-cleavages", "1",
        "--min-length", "5", "--max-length", "30", "--charges", "2",
        "--quiet", "--no-progress",
    ])
    assert res_d.exit_code == 0, res_d.output

    res_rt = runner.invoke(app, [
        "predict-rt", "--precursors", str(prec), "--output", str(prec_rt),
        "--backend", "mock", "--rt-output-column", "NormalizedRetentionTime",
        "--write-deeplc-input", str(deeplc_input), "--quiet", "--no-progress",
    ])
    assert res_rt.exit_code == 0, res_rt.output
    df_rt = pd.read_csv(prec_rt, sep="\t")
    assert "NormalizedRetentionTime" in df_rt.columns
    assert df_rt["NormalizedRetentionTime"].notna().all()
    assert deeplc_input.exists()
    dlc_in = pd.read_csv(deeplc_input, sep="\t")
    assert {"seq", "modifications"}.issubset(dlc_in.columns)
    assert dlc_in["modifications"].astype(str).str.contains("Carbamidomethyl").any()

    res_p = runner.invoke(app, [
        "predict", "--precursors", str(prec_rt), "--output", str(pred),
        "--backend", "mock", "--fragment-types", "b,y", "--fragment-charges", "1",
        "--top-n", "3", "--include-rt", "--rt-column", "NormalizedRetentionTime",
        "--quiet", "--no-progress",
    ])
    assert res_p.exit_code == 0, res_p.output
    pred_df = pd.read_csv(pred, sep="\t")
    assert "NormalizedRetentionTime" in pred_df.columns
    assert pred_df["NormalizedRetentionTime"].notna().all()
