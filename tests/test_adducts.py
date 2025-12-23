"""Tests for chemistry/adducts.py module."""

from lib_eic.chemistry.adducts import (
    ADDUCT_DEFINITIONS,
    MASS_H,
    MASS_Na,
    MASS_NH4,
    MASS_H2O,
    MASS_E,
    get_enabled_adducts,
)


class TestMassConstants:
    """Tests for mass constants."""

    def test_hydrogen_mass(self):
        """Test hydrogen mass is approximately correct."""
        assert abs(MASS_H - 1.0073) < 0.001

    def test_sodium_mass(self):
        """Test sodium mass is approximately correct."""
        assert abs(MASS_Na - 22.9892) < 0.001

    def test_ammonium_mass(self):
        """Test ammonium mass is approximately correct."""
        assert abs(MASS_NH4 - 18.0338) < 0.001

    def test_water_mass(self):
        """Test water mass is approximately correct."""
        assert abs(MASS_H2O - 18.0106) < 0.001

    def test_electron_mass(self):
        """Test electron mass is approximately correct."""
        assert abs(MASS_E - 0.00054858) < 1e-7


class TestAdductDefinitions:
    """Tests for adduct definitions."""

    def test_contains_mh_plus(self):
        """Test [M+H]+ adduct is defined."""
        assert "[M+H]+" in ADDUCT_DEFINITIONS
        adduct = ADDUCT_DEFINITIONS["[M+H]+"]
        assert adduct["multiplier"] == 1
        assert adduct["net_charge"] == 1
        assert adduct["delta"] == MASS_H

    def test_contains_mh_minus(self):
        """Test [M-H]- adduct is defined."""
        assert "[M-H]-" in ADDUCT_DEFINITIONS
        adduct = ADDUCT_DEFINITIONS["[M-H]-"]
        assert adduct["multiplier"] == 1
        assert adduct["net_charge"] == -1
        assert adduct["delta"] == -MASS_H

    def test_all_adducts_have_required_keys(self):
        """Test all adducts have required keys."""
        required_keys = {"multiplier", "delta", "net_charge"}
        for name, adduct in ADDUCT_DEFINITIONS.items():
            assert required_keys.issubset(adduct.keys()), f"Adduct {name} missing keys"

    def test_charges_are_valid(self):
        """Test all adduct charges are +1 or -1."""
        for name, adduct in ADDUCT_DEFINITIONS.items():
            assert adduct["net_charge"] in (-1, 1), f"Invalid charge for {name}"


class TestGetEnabledAdducts:
    """Tests for get_enabled_adducts function."""

    def test_default_enabled_adducts(self):
        """Test default enabled adducts."""
        enabled = get_enabled_adducts()
        # By default, [M+H]+ and [M-H]- should be enabled
        assert "[M+H]+" in enabled
        assert "[M-H]-" in enabled

    def test_filter_by_positive_mode(self):
        """Test filtering by positive ionization mode."""
        enabled = get_enabled_adducts(mode="POS")
        for name, info in enabled.items():
            assert info["net_charge"] > 0, f"{name} should be positive"

    def test_filter_by_negative_mode(self):
        """Test filtering by negative ionization mode."""
        enabled = get_enabled_adducts(mode="NEG")
        for name, info in enabled.items():
            assert info["net_charge"] < 0, f"{name} should be negative"

    def test_explicit_adduct_selection(self):
        """Test explicit adduct selection."""
        selected = ["[M+H]+", "[M+Na]+"]
        enabled = get_enabled_adducts(enabled_names=selected)
        assert set(enabled.keys()) == set(selected)

    def test_explicit_selection_with_mode_filter(self):
        """Test explicit selection combined with mode filter."""
        selected = ["[M+H]+", "[M-H]-"]  # One positive, one negative
        enabled = get_enabled_adducts(enabled_names=selected, mode="POS")
        assert "[M+H]+" in enabled
        assert "[M-H]-" not in enabled

    def test_returns_without_enabled_key(self):
        """Test returned adducts don't have 'enabled' key."""
        enabled = get_enabled_adducts()
        for name, info in enabled.items():
            assert "enabled" not in info
