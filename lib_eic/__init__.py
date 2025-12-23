"""
LCMS Adduct Finder - Automated Targeted Feature Extraction & Adduct Verification Tool for LC-MS Data.
"""

__version__ = "0.2.2"

from .config import Config, load_config
from .chemistry.mass import get_exact_mass, calculate_target_mz
from .chemistry.adducts import ADDUCT_DEFINITIONS
from .validation import validate_mode, validate_formula

__all__ = [
    "Config",
    "load_config",
    "get_exact_mass",
    "calculate_target_mz",
    "ADDUCT_DEFINITIONS",
    "validate_mode",
    "validate_formula",
    "__version__",
]
