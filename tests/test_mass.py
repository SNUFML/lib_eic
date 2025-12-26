"""Tests for chemistry/mass.py module."""

from lib_eic.chemistry.mass import (
    normalize_formula_str,
    get_exact_mass,
    calculate_target_mz_from_mass,
    calculate_target_mz,
    ppm_to_da,
    clear_mass_cache,
)
from lib_eic.chemistry.adducts import MASS_H, MASS_E


class TestNormalizeFormulaStr:
    """Tests for normalize_formula_str function."""

    def test_removes_whitespace(self):
        assert normalize_formula_str("C6 H12 O6") == "C6H12O6"

    def test_removes_tabs_and_newlines(self):
        assert normalize_formula_str("C6\tH12\nO6") == "C6H12O6"

    def test_handles_none(self):
        assert normalize_formula_str(None) == ""

    def test_handles_empty_string(self):
        assert normalize_formula_str("") == ""

    def test_preserves_valid_formula(self):
        assert normalize_formula_str("C6H12O6") == "C6H12O6"


class TestGetExactMass:
    """Tests for get_exact_mass function."""

    def setup_method(self):
        """Clear cache before each test."""
        clear_mass_cache()

    def test_glucose_mass(self):
        """Test mass calculation for glucose (C6H12O6)."""
        mass = get_exact_mass("C6H12O6")
        assert mass is not None
        # Glucose monoisotopic mass is approximately 180.0634
        assert abs(mass - 180.0634) < 0.001

    def test_water_mass(self):
        """Test mass calculation for water (H2O)."""
        mass = get_exact_mass("H2O")
        assert mass is not None
        # Water monoisotopic mass is approximately 18.0106
        assert abs(mass - 18.0106) < 0.001

    def test_caffeine_mass(self):
        """Test mass calculation for caffeine (C8H10N4O2)."""
        mass = get_exact_mass("C8H10N4O2")
        assert mass is not None
        # Caffeine monoisotopic mass is approximately 194.0804
        assert abs(mass - 194.0804) < 0.001

    def test_empty_formula_returns_none(self):
        assert get_exact_mass("") is None

    def test_none_formula_returns_none(self):
        assert get_exact_mass(None) is None

    def test_invalid_formula_returns_none(self):
        assert get_exact_mass("InvalidFormula123XYZ") is None

    def test_caching_works(self):
        """Test that caching works correctly."""
        # First call
        mass1 = get_exact_mass("C6H12O6")
        # Second call should return cached value
        mass2 = get_exact_mass("C6H12O6")
        assert mass1 == mass2


class TestCalculateTargetMz:
    """Tests for m/z calculation functions."""

    def test_mh_plus_adduct(self):
        """Test [M+H]+ adduct calculation."""
        exact_mass = 180.0634  # Glucose
        adduct_info = {"multiplier": 1, "delta": MASS_H, "net_charge": 1}

        target_mz = calculate_target_mz_from_mass(exact_mass, adduct_info)

        assert target_mz is not None
        # Expected: (180.0634 + 1.0073 - 0.00054858) / 1 = ~181.0702
        expected = (exact_mass + MASS_H - MASS_E) / 1
        assert abs(target_mz - expected) < 0.0001

    def test_mh_minus_adduct(self):
        """Test [M-H]- adduct calculation."""
        exact_mass = 180.0634  # Glucose
        adduct_info = {"multiplier": 1, "delta": -MASS_H, "net_charge": -1}

        target_mz = calculate_target_mz_from_mass(exact_mass, adduct_info)

        assert target_mz is not None
        # Expected: (180.0634 - 1.0073 - (-1) * 0.00054858) / 1 = ~179.0566
        expected = (exact_mass - MASS_H + MASS_E) / 1
        assert abs(target_mz - expected) < 0.0001

    def test_dimer_adduct(self):
        """Test [2M+H]+ dimer adduct calculation."""
        exact_mass = 180.0634  # Glucose
        adduct_info = {"multiplier": 2, "delta": MASS_H, "net_charge": 1}

        target_mz = calculate_target_mz_from_mass(exact_mass, adduct_info)

        assert target_mz is not None
        expected = (exact_mass * 2 + MASS_H - MASS_E) / 1
        assert abs(target_mz - expected) < 0.0001

    def test_invalid_mass_returns_none(self):
        adduct_info = {"multiplier": 1, "delta": MASS_H, "net_charge": 1}
        assert calculate_target_mz_from_mass(float("nan"), adduct_info) is None
        assert calculate_target_mz_from_mass(float("inf"), adduct_info) is None

    def test_calculate_target_mz_from_formula(self):
        """Test calculate_target_mz convenience function."""
        clear_mass_cache()
        adduct_info = {"multiplier": 1, "delta": MASS_H, "net_charge": 1}

        target_mz = calculate_target_mz("C6H12O6", adduct_info)

        assert target_mz is not None
        assert target_mz > 180  # Should be slightly higher due to H+


class TestPpmToDa:
    """Tests for ppm_to_da function."""

    def test_10ppm_at_500mz(self):
        """Test 10 ppm tolerance at 500 m/z."""
        tolerance = ppm_to_da(500.0, 10.0)
        expected = 500.0 * (10.0 / 1e6)  # 0.005 Da
        assert abs(tolerance - expected) < 1e-10

    def test_5ppm_at_1000mz(self):
        """Test 5 ppm tolerance at 1000 m/z."""
        tolerance = ppm_to_da(1000.0, 5.0)
        expected = 1000.0 * (5.0 / 1e6)  # 0.005 Da
        assert abs(tolerance - expected) < 1e-10
