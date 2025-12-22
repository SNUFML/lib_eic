"""Adduct definitions and mass constants for LC-MS analysis."""

from typing import Any, Dict, List, Optional

# Atomic/molecular mass constants (isotopic masses)
MASS_H = 1.0073  # Hydrogen
MASS_Na = 22.9892  # Sodium
MASS_ACN = 41.0265  # Acetonitrile
MASS_FA = 46.0055  # Formic acid
MASS_HCOO = 44.9983  # Formate
MASS_NH4 = 18.0338  # Ammonium
MASS_H2O = 18.0106  # Water
MASS_E = 0.00054858  # Electron mass


# Adduct definitions with multiplier, delta mass, and net charge
# Formula: ion_mass = (exact_mass * multiplier) + delta - (charge * MASS_E)
#          target_mz = ion_mass / abs(charge)

ADDUCT_DEFINITIONS: Dict[str, Dict[str, Any]] = {
    # ===== Dimer Adducts (n=2) =====
    "[2M-2H+Na]-": {
        "multiplier": 2,
        "delta": -2 * MASS_H + MASS_Na,
        "net_charge": -1,
        "enabled": False,
    },
    "[2M-H]-": {
        "multiplier": 2,
        "delta": -MASS_H,
        "net_charge": -1,
        "enabled": False,
    },
    "[2M+H]+": {
        "multiplier": 2,
        "delta": +MASS_H,
        "net_charge": +1,
        "enabled": False,
    },
    # ===== Monomer Adducts (n=1) =====
    "[M-2H2O+H]+": {
        "multiplier": 1,
        "delta": -(2 * MASS_H2O) + MASS_H,
        "net_charge": +1,
        "enabled": False,
    },
    "[M-3H2O+H]+": {
        "multiplier": 1,
        "delta": -(3 * MASS_H2O) + MASS_H,
        "net_charge": +1,
        "enabled": False,
    },
    "[M-H]-": {
        "multiplier": 1,
        "delta": -MASS_H,
        "net_charge": -1,
        "enabled": True,  # Default enabled
    },
    "[M-H2O-H]-": {
        "multiplier": 1,
        "delta": -(MASS_H2O + MASS_H),
        "net_charge": -1,
        "enabled": False,
    },
    "[M-H2O+H]+": {
        "multiplier": 1,
        "delta": -MASS_H2O + MASS_H,
        "net_charge": +1,
        "enabled": False,
    },
    "[M]-": {
        "multiplier": 1,
        "delta": 0,
        "net_charge": -1,
        "enabled": False,
    },
    "[M]+": {
        "multiplier": 1,
        "delta": 0,
        "net_charge": +1,
        "enabled": False,
    },
    "[M+ACN+H]+": {
        "multiplier": 1,
        "delta": +MASS_ACN + MASS_H,
        "net_charge": +1,
        "enabled": False,
    },
    "[M+FA-H]-": {
        "multiplier": 1,
        "delta": +(MASS_FA - MASS_H),
        "net_charge": -1,
        "enabled": False,
    },
    "[M+H]+": {
        "multiplier": 1,
        "delta": +MASS_H,
        "net_charge": +1,
        "enabled": True,  # Default enabled
    },
    "[M+Na]+": {
        "multiplier": 1,
        "delta": +MASS_Na,
        "net_charge": +1,
        "enabled": False,
    },
    "[M+NH4]+": {
        "multiplier": 1,
        "delta": +MASS_NH4,
        "net_charge": +1,
        "enabled": False,
    },
}


def get_enabled_adducts(
    enabled_names: Optional[List[str]] = None,
    mode: Optional[str] = None,
) -> Dict[str, Dict[str, Any]]:
    """Get adduct definitions filtered by enabled status and ionization mode.

    Args:
        enabled_names: List of adduct names to enable. If None, uses default 'enabled' flag.
        mode: Ionization mode ("POS" or "NEG"). If None, returns all.

    Returns:
        Dictionary of enabled adduct definitions.
    """
    result = {}

    for name, info in ADDUCT_DEFINITIONS.items():
        # Check if enabled
        if enabled_names is not None:
            if name not in enabled_names:
                continue
        elif not info.get("enabled", False):
            continue

        # Check ionization mode
        if mode is not None:
            charge = info["net_charge"]
            if mode == "POS" and charge < 0:
                continue
            if mode == "NEG" and charge > 0:
                continue

        # Return without 'enabled' key for compatibility
        result[name] = {
            "multiplier": info["multiplier"],
            "delta": info["delta"],
            "net_charge": info["net_charge"],
        }

    return result
