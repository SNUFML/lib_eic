"""Analysis modules for LC-MS data processing."""

from .eic import DirectTarget, Target, extract_eic, extract_eics_multi
from .fitting import fit_gaussian_and_score, gaussian_func
from .ms2 import build_ms2_index, match_ms2

__all__ = [
    "DirectTarget",
    "Target",
    "extract_eic",
    "extract_eics_multi",
    "fit_gaussian_and_score",
    "gaussian_func",
    "build_ms2_index",
    "match_ms2",
]
