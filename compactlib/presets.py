from __future__ import annotations

from copy import deepcopy
from typing import Any

DIGEST_DEFAULTS: dict[str, Any] = {
    "enzyme": "trypsin-p",
    "missed_cleavages": 1,
    "min_length": 7,
    "max_length": 30,
    "charges": "2,3,4",
    "min_precursor_mz": None,
    "max_precursor_mz": None,
    "charge_length_rules": "",
}

PREDICT_DEFAULTS: dict[str, Any] = {
    "fragment_types": "b,y",
    "fragment_charges": "1",
    "min_fragment_series": None,
    "max_fragment_series": None,
    "min_product_mz": None,
    "max_product_mz": None,
    "recommended_top_n": None,
}

# Broad generic DIA precursor-space preset. It was calibrated against a Skyline
# precursor export produced from a Swiss-Prot FASTA in the compact DIA library
# manuscript, but the rules themselves are not organism-specific:
#   peptide length: 7-30
#   precursor charges: 1,2,3,4
#   precursor m/z: approximately 300-1800
#   maximum missed cleavages: 1
#   fixed C carbamidomethylation, no variable M oxidation by default
#
# Importantly, the audited Skyline export did not apply a hard
# "charge 4 only above length X" rule. Charge 4 candidates started at peptide
# length 8, and the observed charge/length pattern is mainly explained by the
# precursor m/z filter.
GENERIC_DIA_DIGEST_PRESET: dict[str, Any] = {
    "enzyme": "trypsin-p",
    "missed_cleavages": 1,
    "min_length": 7,
    "max_length": 30,
    "charges": "1,2,3,4",
    "min_precursor_mz": 300.0,
    "max_precursor_mz": 1800.0,
    "charge_length_rules": "",
}

GENERIC_DIA_WINDOWED_DIGEST_PRESET: dict[str, Any] = {
    "enzyme": "trypsin-p",
    "missed_cleavages": 1,
    "min_length": 7,
    "max_length": 30,
    "charges": "2,3,4",
    "min_precursor_mz": 400.0,
    "max_precursor_mz": 1200.0,
    "charge_length_rules": "",
}

DIGEST_PRESETS: dict[str, dict[str, Any]] = {
    "none": {},
    # Neutral default for broad FASTA-derived DIA library generation. This is
    # the recommended cross-organism preset for human/yeast/E. coli-style use.
    "generic-dia": deepcopy(GENERIC_DIA_DIGEST_PRESET),
    # Narrower DIA-window-focused preset for users who want to restrict the
    # candidate space to a common acquisition m/z range.
    "generic-dia-windowed": deepcopy(GENERIC_DIA_WINDOWED_DIGEST_PRESET),
    # Windowed preset plus optional plausibility rule: charge 4 only for peptide
    # length >=12. This is a heuristic, not an exact Skyline or DIA-NN setting.
    "generic-dia-strict": {
        **deepcopy(GENERIC_DIA_WINDOWED_DIGEST_PRESET),
        "charge_length_rules": "4:12",
    },
    # Practical DIA-NN-like candidate-space preset. It is not intended to be an
    # exact clone of any specific DIA-NN release, but captures the common idea:
    # tryptic DIA candidate space with restricted peptide length, precursor
    # charges and precursor m/z range.
    "diann-like": deepcopy(GENERIC_DIA_WINDOWED_DIGEST_PRESET),
    # Skyline-like precursor export preset calibrated from the user's Skyline
    # digest audit. This reproduces the broad Skyline candidate space rather than
    # a narrow DIA-window-specific precursor list.
    "skyline-dia-like": deepcopy(GENERIC_DIA_DIGEST_PRESET),
    # Explicit alias for the audited Skyline precursor export.
    "skyline-export-like": deepcopy(GENERIC_DIA_DIGEST_PRESET),
}

GENERIC_DIA_PREDICT_PRESET: dict[str, Any] = {
    "fragment_types": "b,y",
    "fragment_charges": "1,2,3",
    "min_fragment_series": 3,
    "max_fragment_series": None,
    "min_product_mz": None,
    "max_product_mz": None,
    "recommended_top_n": 7,
}

SKYLINE_PREDICT_PRESET: dict[str, Any] = {
    "fragment_types": "b,y",
    "fragment_charges": "1,2",
    "min_fragment_series": 3,
    "max_fragment_series": None,
    "min_product_mz": None,
    "max_product_mz": None,
    "recommended_top_n": 6,
}

PREDICT_PRESETS: dict[str, dict[str, Any]] = {
    "none": {},
    "generic-dia": deepcopy(GENERIC_DIA_PREDICT_PRESET),
    "generic-dia-windowed": deepcopy(GENERIC_DIA_PREDICT_PRESET),
    "generic-dia-strict": deepcopy(GENERIC_DIA_PREDICT_PRESET),
    "diann-like": {
        "fragment_types": "b,y",
        "fragment_charges": "1,2",
        "min_fragment_series": 3,
        "max_fragment_series": None,
        "min_product_mz": None,
        "max_product_mz": None,
        "recommended_top_n": 7,
    },
    "skyline-dia-like": deepcopy(SKYLINE_PREDICT_PRESET),
    "skyline-export-like": deepcopy(SKYLINE_PREDICT_PRESET),
}

def available_presets() -> list[str]:
    return sorted(DIGEST_PRESETS.keys())


def validate_preset_name(name: str) -> str:
    name = str(name or "none").strip().lower()
    if name not in DIGEST_PRESETS:
        allowed = ", ".join(available_presets())
        raise ValueError(f"Unknown preset '{name}'. Available presets: {allowed}")
    return name


def resolve_digest_preset(preset: str = "none", **overrides: Any) -> dict[str, Any]:
    """Resolve digest parameters from defaults, preset, and non-None overrides."""
    preset = validate_preset_name(preset)
    out = deepcopy(DIGEST_DEFAULTS)
    out.update(deepcopy(DIGEST_PRESETS[preset]))
    for key, value in overrides.items():
        if value is not None:
            out[key] = value
    out["preset"] = preset
    return out


def resolve_predict_preset(preset: str = "none", **overrides: Any) -> dict[str, Any]:
    """Resolve transition/prediction filter parameters from defaults, preset, and non-None overrides."""
    preset = validate_preset_name(preset)
    out = deepcopy(PREDICT_DEFAULTS)
    out.update(deepcopy(PREDICT_PRESETS[preset]))
    for key, value in overrides.items():
        if value is not None:
            out[key] = value
    out["preset"] = preset
    return out


def preset_table() -> str:
    lines: list[str] = []
    for name in available_presets():
        if name == "none":
            continue
        d = {**DIGEST_DEFAULTS, **DIGEST_PRESETS.get(name, {})}
        p = {**PREDICT_DEFAULTS, **PREDICT_PRESETS.get(name, {})}
        lines.append(f"{name}:")
        lines.append(f"  digest enzyme: {d['enzyme']}")
        lines.append(f"  missed cleavages: {d['missed_cleavages']}")
        lines.append(f"  peptide length: {d['min_length']}-{d['max_length']}")
        lines.append(f"  precursor charges: {d['charges']}")
        lines.append(f"  precursor m/z: {d['min_precursor_mz']}-{d['max_precursor_mz']}")
        lines.append(f"  charge-length rules: {d['charge_length_rules'] or 'none'}")
        lines.append(f"  fragment types: {p['fragment_types']}")
        lines.append(f"  fragment charges: {p['fragment_charges']}")
        lines.append(f"  min fragment series: {p['min_fragment_series'] or 'none'}")
        lines.append(f"  recommended top-n: {p['recommended_top_n'] or 'none'}")
        lines.append("")
    return "\n".join(lines).rstrip()
