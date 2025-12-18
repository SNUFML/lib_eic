"""Mass calculation functions for LC-MS analysis."""

import math
import re
from functools import lru_cache
from typing import Dict, Any, Optional

from .adducts import MASS_E

# Regex for removing whitespace from formula strings
_FORMULA_WS_RE = re.compile(r"\s+")


def normalize_formula_str(value: str) -> str:
    """Normalize a chemical formula string by removing whitespace.

    Args:
        value: Chemical formula string (e.g., "C6 H12 O6" or "C6H12O6").

    Returns:
        Normalized formula string without whitespace.
    """
    if value is None:
        return ""
    value = str(value)
    value = _FORMULA_WS_RE.sub("", value).strip()
    return value


@lru_cache(maxsize=1024)
def get_exact_mass(formula_str: str) -> Optional[float]:
    """Calculate the exact monoisotopic mass of a chemical formula.

    Uses LRU cache to avoid recalculating masses for repeated formulas.

    Args:
        formula_str: Chemical formula string (e.g., "C6H12O6").

    Returns:
        Exact monoisotopic mass in Daltons, or None if formula is invalid.
    """
    formula_norm = normalize_formula_str(formula_str)
    if not formula_norm:
        return None

    try:
        from molmass import Formula

        mass = float(Formula(formula_norm).isotope.mass)
        return mass
    except Exception:
        return None


def calculate_target_mz_from_mass(
    exact_mass: float,
    adduct_info: Dict[str, Any],
) -> Optional[float]:
    """Calculate target m/z from exact mass and adduct information.

    Formula: ion_mass = (exact_mass * multiplier) + delta - (charge * MASS_E)
             target_mz = ion_mass / abs(charge)

    Args:
        exact_mass: Exact monoisotopic mass of the molecule.
        adduct_info: Dictionary with 'multiplier', 'delta', and 'net_charge' keys.

    Returns:
        Target m/z value, or None if calculation fails.
    """
    try:
        exact_mass = float(exact_mass)
        if not math.isfinite(exact_mass):
            return None

        multiplier = adduct_info["multiplier"]
        delta = adduct_info["delta"]
        charge = adduct_info["net_charge"]

        ion_mass = (exact_mass * multiplier) + delta - (charge * MASS_E)
        return float(ion_mass / abs(charge))
    except (KeyError, TypeError, ZeroDivisionError):
        return None


def calculate_target_mz(
    formula_str: str,
    adduct_info: Dict[str, Any],
) -> Optional[float]:
    """Calculate target m/z from chemical formula and adduct information.

    Convenience function that combines get_exact_mass and calculate_target_mz_from_mass.

    Args:
        formula_str: Chemical formula string (e.g., "C6H12O6").
        adduct_info: Dictionary with 'multiplier', 'delta', and 'net_charge' keys.

    Returns:
        Target m/z value, or None if calculation fails.
    """
    exact_mass = get_exact_mass(formula_str)
    if exact_mass is None:
        return None
    return calculate_target_mz_from_mass(exact_mass, adduct_info)


def ppm_to_da(mz: float, ppm: float) -> float:
    """Convert mass tolerance from ppm to Daltons.

    Args:
        mz: m/z value.
        ppm: Tolerance in parts per million.

    Returns:
        Tolerance in Daltons.
    """
    return float(mz) * (float(ppm) / 1e6)


def clear_mass_cache() -> None:
    """Clear the exact mass calculation cache."""
    get_exact_mass.cache_clear()
