"""EIC (Extracted Ion Chromatogram) extraction and target building."""

import logging
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

from ..chemistry.mass import (
    calculate_target_mz_from_mass,
    get_exact_mass,
)
from ..chemistry.adducts import get_enabled_adducts
from ..io.raw_file import RawFileReader

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Target:
    """Represents a target for EIC extraction.

    Attributes:
        formula: Chemical formula.
        adduct: Adduct name.
        mz: Target m/z value.
    """

    formula: str
    adduct: str
    mz: float

    @property
    def key(self) -> Tuple[str, str]:
        """Get a unique key for this target."""
        return (self.formula, self.adduct)


@dataclass(frozen=True)
class DirectTarget:
    """Represents a target for EIC extraction using a direct m/z value.

    Attributes:
        compound_name: Compound display name.
        polarity: Ionization polarity ("POS" or "NEG").
        mz: Target m/z value.
        mixture: Mixture identifier.
    """

    compound_name: str
    polarity: str
    mz: float
    mixture: str

    @property
    def key(self) -> Tuple[str, str, str]:
        """Get a unique key for this target: (compound_name, polarity, mixture)."""
        return (self.compound_name, self.polarity, self.mixture)


def build_targets(
    formulas: Iterable[str],
    mode: str,
    enabled_adducts: Optional[List[str]] = None,
) -> List[Target]:
    """Build target list for EIC extraction.

    Args:
        formulas: Iterable of chemical formulas.
        mode: Ionization mode ("POS" or "NEG").
        enabled_adducts: List of adduct names to use. If None, uses defaults.

    Returns:
        List of Target objects.
    """
    targets: List[Target] = []
    formulas_list = list(formulas)

    # Get adducts filtered by mode and enabled status
    adducts = get_enabled_adducts(enabled_names=enabled_adducts, mode=mode)

    for formula in formulas_list:
        exact_mass = get_exact_mass(formula)
        if exact_mass is None:
            logger.warning("Could not calculate mass for formula: %s", formula)
            continue

        for adduct_name, adduct_info in adducts.items():
            target_mz = calculate_target_mz_from_mass(exact_mass, adduct_info)
            if target_mz is None:
                continue

            targets.append(Target(formula=formula, adduct=adduct_name, mz=float(target_mz)))

    logger.debug("Built %d targets from %d formulas", len(targets), len(formulas_list))
    return targets


def extract_eic(
    reader: RawFileReader,
    target_mz: float,
    ppm_tolerance: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Extract a single EIC from a raw file.

    Args:
        reader: RawFileReader instance.
        target_mz: Target m/z value.
        ppm_tolerance: Mass tolerance in ppm.

    Returns:
        Tuple of (rt_minutes, intensity) arrays.
    """
    return reader.get_chromatogram(target_mz, ppm_tolerance)


def extract_eics_multi(
    reader: RawFileReader,
    targets: List[Target],
    ppm_tolerance: float,
    batch_size: int = 256,
) -> Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]]:
    """Extract multiple EICs in batch mode.

    Args:
        reader: RawFileReader instance.
        targets: List of Target objects.
        ppm_tolerance: Mass tolerance in ppm.
        batch_size: Number of chromatograms per batch.

    Returns:
        Dictionary mapping (formula, adduct) to (rt_minutes, intensity).

    Raises:
        RuntimeError: If batch extraction is not available.
    """
    if not reader.has_multi_chromatogram_api():
        raise RuntimeError("Multi-chromatogram API is not available")

    target_mzs = [t.mz for t in targets]
    results = reader.get_chromatograms_batch(target_mzs, ppm_tolerance, batch_size)

    return {t.key: result for t, result in zip(targets, results)}


def calculate_area(
    rt_min: np.ndarray,
    intensity: np.ndarray,
    method: str = "sum",
) -> float:
    """Calculate peak area from EIC data.

    Args:
        rt_min: Retention time array (minutes).
        intensity: Intensity array.
        method: Area calculation method ("sum" or "trapz").

    Returns:
        Calculated area value.
    """
    if not intensity.size:
        return 0.0

    if method == "trapz":
        if intensity.size < 2:
            return 0.0

        # NumPy 2.0+ renamed `trapz` -> `trapezoid` (and may deprecate `trapz`).
        # Implement a small fallback to avoid relying on optional attributes.
        if hasattr(np, "trapezoid"):
            return float(np.trapezoid(intensity, rt_min))

        dx = np.diff(rt_min)
        return float(np.sum((intensity[1:] + intensity[:-1]) * dx / 2.0))

    return float(np.sum(intensity))


def dedupe_preserve_order(items: Iterable[str]) -> List[str]:
    """Remove duplicates while preserving order.

    Args:
        items: Iterable of strings.

    Returns:
        List with duplicates removed, order preserved.
    """
    seen = set()
    out: List[str] = []
    for item in items:
        if item in seen:
            continue
        seen.add(item)
        out.append(item)
    return out
