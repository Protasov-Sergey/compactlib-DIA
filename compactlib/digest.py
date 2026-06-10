from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable

import pandas as pd

from .mass import AA_MASS, format_modified_peptide, format_positions_1based, precursor_mz

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
SELENOCYSTEINE_AA = "U"


@dataclass(frozen=True)
class ProteinRecord:
    protein_id: str
    description: str
    sequence: str


def parse_charges(charges: str | Iterable[int]) -> list[int]:
    if isinstance(charges, str):
        out = [int(x.strip()) for x in charges.split(",") if x.strip()]
    else:
        out = [int(x) for x in charges]
    if not out:
        raise ValueError("No precursor charges provided")
    if any(z <= 0 for z in out):
        raise ValueError("All charges must be positive integers")
    return sorted(set(out))




def parse_charge_length_rules(rules: str | dict[int, int] | None) -> dict[int, int]:
    """Parse minimum peptide length rules per precursor charge.

    Examples
    --------
    ``"4:12,5:16"`` means that charge 4 is allowed only for peptides with
    length >= 12 and charge 5 only for peptides with length >= 16.

    The rules are intentionally optional because different DIA workflows use
    different precursor charge ranges. Precursor m/z filtering is usually the
    first and most important practical constraint.
    """
    if rules is None:
        return {}
    if isinstance(rules, dict):
        out = {int(k): int(v) for k, v in rules.items()}
    else:
        text = str(rules).strip()
        if not text:
            return {}
        out = {}
        for part in text.split(","):
            part = part.strip()
            if not part:
                continue
            if ":" not in part:
                raise ValueError("charge_length_rules must look like '4:12,5:16'")
            z, min_len = part.split(":", 1)
            out[int(z.strip())] = int(min_len.strip())
    if any(z <= 0 for z in out):
        raise ValueError("charge_length_rules charges must be positive")
    if any(v <= 0 for v in out.values()):
        raise ValueError("charge_length_rules minimum lengths must be positive")
    return out

def _protein_id_from_header(header: str) -> str:
    first = header.split()[0]
    parts = first.split("|")
    if len(parts) >= 2 and parts[1]:
        return parts[1]
    return first


def parse_fasta(path: str | Path) -> list[ProteinRecord]:
    path = Path(path).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"FASTA file does not exist: {path}")

    records: list[ProteinRecord] = []
    header: str | None = None
    seq_parts: list[str] = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = "".join(seq_parts).replace("*", "").upper()
                    records.append(ProteinRecord(_protein_id_from_header(header), header, seq))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)

    if header is not None:
        seq = "".join(seq_parts).replace("*", "").upper()
        records.append(ProteinRecord(_protein_id_from_header(header), header, seq))

    return records


def cleavage_sites(sequence: str, enzyme: str = "trypsin-p") -> list[int]:
    enzyme = enzyme.lower()
    if enzyme not in {"trypsin", "trypsin-p", "trypsinp"}:
        raise ValueError("Only enzyme='trypsin' and enzyme='trypsin-p' are supported")

    seq = sequence.upper()
    sites = [0]
    for i, aa in enumerate(seq):
        if aa in {"K", "R"}:
            if enzyme in {"trypsin-p", "trypsinp"}:
                sites.append(i + 1)
            else:
                next_aa = seq[i + 1] if i + 1 < len(seq) else ""
                if next_aa != "P":
                    sites.append(i + 1)
    if sites[-1] != len(seq):
        sites.append(len(seq))
    return sorted(set(sites))


def digest_sequence(sequence: str, enzyme: str = "trypsin-p", missed_cleavages: int = 1) -> list[tuple[str, int]]:
    if missed_cleavages < 0:
        raise ValueError("missed_cleavages must be >= 0")

    sites = cleavage_sites(sequence, enzyme=enzyme)
    peptides: list[tuple[str, int]] = []
    for i in range(len(sites) - 1):
        max_j = min(i + missed_cleavages + 1, len(sites) - 1)
        for j in range(i + 1, max_j + 1):
            pep = sequence[sites[i]:sites[j]]
            mc = j - i - 1
            if pep:
                peptides.append((pep, mc))
    return peptides


def oxidation_m_position_sets(sequence: str, enabled: bool = False, max_variable_mods: int = 1) -> list[tuple[int, ...]]:
    """Return zero-based M-oxidation position combinations.

    The unmodified form is always included. If enabled, all combinations with
    1..max_variable_mods oxidized methionines are included.
    """
    if max_variable_mods < 0:
        raise ValueError("max_variable_mods must be >= 0")
    m_positions = [i for i, aa in enumerate(sequence.upper()) if aa == "M"]
    out: list[tuple[int, ...]] = [tuple()]
    if not enabled or max_variable_mods == 0 or not m_positions:
        return out
    for k in range(1, min(max_variable_mods, len(m_positions)) + 1):
        out.extend(tuple(c) for c in combinations(m_positions, k))
    return out


def build_precursor_table_from_fasta(
    fasta: str | Path,
    enzyme: str = "trypsin-p",
    missed_cleavages: int = 1,
    min_length: int = 7,
    max_length: int = 30,
    charges: str | Iterable[int] = "2,3,4",
    min_precursor_mz: float | None = None,
    max_precursor_mz: float | None = None,
    charge_length_rules: str | dict[int, int] | None = None,
    carbamidomethyl_c: bool = True,
    c_mod_format: str = "unimod",
    variable_oxidation_m: bool = False,
    max_variable_mods: int = 1,
    m_mod_format: str = "unimod",
    remove_invalid_aa: bool = True,
    allow_selenocysteine: bool = False,
) -> tuple[pd.DataFrame, dict]:
    """Create a precursor table from FASTA.

    Variable M oxidation is expanded at the peptidoform level. Modification
    positions are encoded as 1-based semicolon-separated positions in the
    ``OxidationMPositions`` output column.
    """
    if min_length <= 0:
        raise ValueError("min_length must be positive")
    if max_length < min_length:
        raise ValueError("max_length must be >= min_length")

    z_values = parse_charges(charges)
    charge_len_rules = parse_charge_length_rules(charge_length_rules)
    if min_precursor_mz is not None and max_precursor_mz is not None and max_precursor_mz < min_precursor_mz:
        raise ValueError("max_precursor_mz must be >= min_precursor_mz")
    records = parse_fasta(fasta)

    allowed_aa = set(CANONICAL_AA)
    if allow_selenocysteine:
        allowed_aa.add(SELENOCYSTEINE_AA)

    peptide_map: dict[str, dict] = {}
    n_raw_peptides = 0
    n_length_filtered = 0
    n_invalid_filtered = 0

    for rec in records:
        for pep, mc in digest_sequence(rec.sequence, enzyme=enzyme, missed_cleavages=missed_cleavages):
            n_raw_peptides += 1
            if len(pep) < min_length or len(pep) > max_length:
                n_length_filtered += 1
                continue
            invalid = sorted(set(pep) - allowed_aa)
            if invalid and remove_invalid_aa:
                n_invalid_filtered += 1
                continue

            item = peptide_map.setdefault(
                pep,
                {
                    "protein_ids": set(),
                    "protein_descriptions": set(),
                    "min_missed_cleavages": mc,
                    "max_missed_cleavages": mc,
                },
            )
            item["protein_ids"].add(rec.protein_id)
            item["protein_descriptions"].add(rec.description)
            item["min_missed_cleavages"] = min(item["min_missed_cleavages"], mc)
            item["max_missed_cleavages"] = max(item["max_missed_cleavages"], mc)

    rows = []
    n_modified_peptides = 0
    n_oxidized_modified_peptides = 0
    n_precursors_before_precursor_filters = 0
    n_precursors_filtered_charge_length = 0
    n_precursors_filtered_mz = 0
    for pep in sorted(peptide_map):
        item = peptide_map[pep]
        protein_ids = sorted(item["protein_ids"])
        protein_desc = sorted(item["protein_descriptions"])
        ox_sets = oxidation_m_position_sets(
            pep,
            enabled=variable_oxidation_m,
            max_variable_mods=max_variable_mods,
        )
        for ox_positions in ox_sets:
            n_modified_peptides += 1
            if ox_positions:
                n_oxidized_modified_peptides += 1
            modpep = format_modified_peptide(
                pep,
                carbamidomethyl_c=carbamidomethyl_c,
                c_mod_format=c_mod_format,
                oxidation_m_positions=ox_positions,
                m_mod_format=m_mod_format,
            )
            ox_pos_1based = format_positions_1based(ox_positions)
            for z in z_values:
                n_precursors_before_precursor_filters += 1
                min_len_for_charge = charge_len_rules.get(z)
                if min_len_for_charge is not None and len(pep) < min_len_for_charge:
                    n_precursors_filtered_charge_length += 1
                    continue

                mz = precursor_mz(
                    pep,
                    z,
                    carbamidomethyl_c=carbamidomethyl_c,
                    oxidation_m_positions=ox_positions,
                )
                if min_precursor_mz is not None and mz < float(min_precursor_mz):
                    n_precursors_filtered_mz += 1
                    continue
                if max_precursor_mz is not None and mz > float(max_precursor_mz):
                    n_precursors_filtered_mz += 1
                    continue

                rows.append(
                    {
                        "ModifiedPeptide": modpep,
                        "StrippedPeptide": pep,
                        "PrecursorCharge": z,
                        "PrecursorMz": mz,
                        "ProteinId": ";".join(protein_ids),
                        "ProteinCount": len(protein_ids),
                        "ProteinDescription": ";".join(protein_desc),
                        "PeptideLength": len(pep),
                        "MinMissedCleavages": item["min_missed_cleavages"],
                        "MaxMissedCleavages": item["max_missed_cleavages"],
                        "FixedCarbamidomethylC": bool(carbamidomethyl_c),
                        "VariableOxidationM": bool(ox_positions),
                        "OxidationMPositions": ox_pos_1based,
                        "NVariableMods": len(ox_positions),
                        "ContainsSelenocysteine": SELENOCYSTEINE_AA in pep,
                    }
                )

    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(
            ["ModifiedPeptide", "PrecursorCharge", "ProteinId"],
            ascending=[True, True, True],
            kind="mergesort",
        ).reset_index(drop=True)

    stats = {
        "mode": "digest",
        "fasta": str(Path(fasta).expanduser()),
        "enzyme": enzyme,
        "missed_cleavages": missed_cleavages,
        "min_length": min_length,
        "max_length": max_length,
        "charges": ",".join(map(str, z_values)),
        "min_precursor_mz": min_precursor_mz,
        "max_precursor_mz": max_precursor_mz,
        "charge_length_rules": ",".join(f"{k}:{v}" for k, v in sorted(charge_len_rules.items())),
        "carbamidomethyl_c": bool(carbamidomethyl_c),
        "c_mod_format": c_mod_format,
        "variable_oxidation_m": bool(variable_oxidation_m),
        "max_variable_mods": max_variable_mods,
        "m_mod_format": m_mod_format,
        "remove_invalid_aa": bool(remove_invalid_aa),
        "allow_selenocysteine": bool(allow_selenocysteine),
        "n_proteins_input": len(records),
        "n_raw_digested_peptides": n_raw_peptides,
        "n_length_filtered_peptides": n_length_filtered,
        "n_invalid_aa_filtered_peptides": n_invalid_filtered,
        "n_unique_stripped_peptides": len(peptide_map),
        "n_unique_selenocysteine_peptides": sum(1 for p in peptide_map if SELENOCYSTEINE_AA in p),
        "n_unique_modified_peptides": n_modified_peptides,
        "n_oxidized_modified_peptides": n_oxidized_modified_peptides,
        "n_precursors_before_precursor_filters": n_precursors_before_precursor_filters,
        "n_precursors_filtered_charge_length": n_precursors_filtered_charge_length,
        "n_precursors_filtered_mz": n_precursors_filtered_mz,
        "n_selenocysteine_precursors_output": int(df["ContainsSelenocysteine"].sum()) if len(df) and "ContainsSelenocysteine" in df.columns else 0,
        "n_precursors_output": int(len(df)),
    }
    return df, stats
