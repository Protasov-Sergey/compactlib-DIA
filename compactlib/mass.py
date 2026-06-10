from __future__ import annotations

from typing import Iterable, Sequence

AA_MASS = {
    "A": 71.037113805,
    "R": 156.101111050,
    "N": 114.042927470,
    "D": 115.026943065,
    "C": 103.009184505,
    "E": 129.042593135,
    "Q": 128.058577540,
    "G": 57.021463735,
    "H": 137.058911875,
    "I": 113.084064015,
    "L": 113.084064015,
    "K": 128.094963050,
    "M": 131.040484645,
    "F": 147.068413945,
    "P": 97.052763875,
    "S": 87.032028435,
    "T": 101.047678505,
    "W": 186.079312980,
    "Y": 163.063328575,
    "V": 99.068413945,
    # Selenocysteine residue mass. Kept for Skyline-compatible digest/mass
    # calculation. Most ML predictors (Prosit/MS2PIP/Koina models) do not
    # support U-containing peptide sequences.
    "U": 150.953633405,
}

H2O_MASS = 18.010564684
PROTON_MASS = 1.007276466621
CARBAMIDOMETHYL_MASS = 57.021463735
OXIDATION_M_MASS = 15.994914620


def _normalize_positions(positions: Iterable[int] | None) -> set[int]:
    """Return zero-based modification positions as a set."""
    if positions is None:
        return set()
    return {int(p) for p in positions}


def residue_mass(
    aa: str,
    position: int | None = None,
    carbamidomethyl_c: bool = True,
    oxidation_m_positions: Iterable[int] | None = None,
) -> float:
    aa = aa.upper()
    if aa not in AA_MASS:
        raise ValueError(f"Unsupported amino acid {aa!r}")
    ox_pos = _normalize_positions(oxidation_m_positions)
    mass = AA_MASS[aa]
    if carbamidomethyl_c and aa == "C":
        mass += CARBAMIDOMETHYL_MASS
    if aa == "M" and position is not None and position in ox_pos:
        mass += OXIDATION_M_MASS
    return mass


def peptide_neutral_mass(
    sequence: str,
    carbamidomethyl_c: bool = True,
    oxidation_m_positions: Iterable[int] | None = None,
) -> float:
    """Return monoisotopic neutral peptide mass including H2O.

    ``oxidation_m_positions`` are zero-based positions in the stripped peptide.
    """
    seq = sequence.strip().upper()
    ox_pos = _normalize_positions(oxidation_m_positions)
    mass = H2O_MASS
    for i, aa in enumerate(seq):
        mass += residue_mass(
            aa,
            position=i,
            carbamidomethyl_c=carbamidomethyl_c,
            oxidation_m_positions=ox_pos,
        )
    return mass


def precursor_mz(
    sequence: str,
    charge: int,
    carbamidomethyl_c: bool = True,
    oxidation_m_positions: Iterable[int] | None = None,
) -> float:
    if charge <= 0:
        raise ValueError("charge must be positive")
    neutral = peptide_neutral_mass(
        sequence,
        carbamidomethyl_c=carbamidomethyl_c,
        oxidation_m_positions=oxidation_m_positions,
    )
    return (neutral + charge * PROTON_MASS) / charge


def fragment_mz(
    sequence: str,
    ion_type: str,
    series_number: int,
    charge: int = 1,
    carbamidomethyl_c: bool = True,
    oxidation_m_positions: Iterable[int] | None = None,
) -> float:
    """Approximate monoisotopic b/y fragment m/z for a stripped sequence.

    This is sufficient for the mock predictor and toy tests. Production Prosit/MS2PIP
    backends should use model-provided ProductMz values.
    """
    if charge <= 0:
        raise ValueError("fragment charge must be positive")
    seq = sequence.strip().upper()
    if series_number <= 0 or series_number >= len(seq):
        raise ValueError("series_number must be in 1..len(sequence)-1")
    ion_type = ion_type.lower()
    ox_pos = _normalize_positions(oxidation_m_positions)

    if ion_type == "b":
        indices = range(0, series_number)
        neutral = sum(
            residue_mass(seq[i], i, carbamidomethyl_c=carbamidomethyl_c, oxidation_m_positions=ox_pos)
            for i in indices
        )
    elif ion_type == "y":
        indices = range(len(seq) - series_number, len(seq))
        neutral = H2O_MASS + sum(
            residue_mass(seq[i], i, carbamidomethyl_c=carbamidomethyl_c, oxidation_m_positions=ox_pos)
            for i in indices
        )
    else:
        raise ValueError("Only b and y ions are supported by fragment_mz")

    return (neutral + charge * PROTON_MASS) / charge


def _format_c(aa: str, c_mod_format: str) -> str:
    if c_mod_format == "unimod":
        return "C(UniMod:4)"
    if c_mod_format == "bracket-unimod":
        return "C[UNIMOD:4]"
    if c_mod_format == "name":
        return "C[Carbamidomethyl]"
    if c_mod_format == "plain":
        return "C"
    raise ValueError("c_mod_format must be one of: unimod, bracket-unimod, name, plain")


def _format_mox(m_mod_format: str) -> str:
    if m_mod_format == "unimod":
        return "M(UniMod:35)"
    if m_mod_format == "bracket-unimod":
        return "M[UNIMOD:35]"
    if m_mod_format == "name":
        return "M[Oxidation]"
    if m_mod_format == "plain":
        return "M"
    raise ValueError("m_mod_format must be one of: unimod, bracket-unimod, name, plain")


def format_modified_peptide(
    sequence: str,
    carbamidomethyl_c: bool = True,
    c_mod_format: str = "unimod",
    oxidation_m_positions: Iterable[int] | None = None,
    m_mod_format: str = "unimod",
) -> str:
    """Format fixed C carbamidomethylation and optional M oxidation.

    Modification position indices are zero-based in the stripped peptide.

    ``*_mod_format`` values:
      - ``unimod``: C(UniMod:4), M(UniMod:35)
      - ``bracket-unimod``: C[UNIMOD:4], M[UNIMOD:35]
      - ``name``: C[Carbamidomethyl], M[Oxidation]
      - ``plain``: do not annotate that modification type
    """
    seq = sequence.strip().upper()
    ox_pos = _normalize_positions(oxidation_m_positions)
    parts: list[str] = []
    for i, aa in enumerate(seq):
        if carbamidomethyl_c and aa == "C" and c_mod_format != "plain":
            parts.append(_format_c(aa, c_mod_format))
        elif aa == "M" and i in ox_pos and m_mod_format != "plain":
            parts.append(_format_mox(m_mod_format))
        else:
            parts.append(aa)
    return "".join(parts)


def parse_positions_1based(value: str | float | int | None) -> tuple[int, ...]:
    """Parse semicolon-separated 1-based positions and return zero-based tuple."""
    if value is None:
        return tuple()
    if isinstance(value, float) and value != value:  # NaN
        return tuple()
    s = str(value).strip()
    if not s or s.lower() in {"nan", "none", "-"}:
        return tuple()
    return tuple(int(float(x)) - 1 for x in s.split(";") if x.strip())


def format_positions_1based(positions: Sequence[int] | Iterable[int]) -> str:
    pos = sorted(int(p) for p in positions)
    return ";".join(str(p + 1) for p in pos)
