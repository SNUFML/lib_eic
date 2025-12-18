"""Tests for validation.py module."""

import pytest
from lib_eic.validation import (
    ValidationError,
    validate_mode,
    validate_formula,
    validate_formulas,
    validate_ppm_tolerance,
)


class TestValidateMode:
    """Tests for validate_mode function."""

    def test_pos_mode(self):
        assert validate_mode("POS") == "POS"
        assert validate_mode("pos") == "POS"
        assert validate_mode("POSITIVE") == "POS"
        assert validate_mode("positive") == "POS"
        assert validate_mode("+") == "POS"

    def test_neg_mode(self):
        assert validate_mode("NEG") == "NEG"
        assert validate_mode("neg") == "NEG"
        assert validate_mode("NEGATIVE") == "NEG"
        assert validate_mode("negative") == "NEG"
        assert validate_mode("-") == "NEG"

    def test_whitespace_handling(self):
        assert validate_mode("  POS  ") == "POS"
        assert validate_mode("\tNEG\n") == "NEG"

    def test_invalid_mode_raises(self):
        with pytest.raises(ValidationError) as exc_info:
            validate_mode("INVALID")
        assert "Invalid ionization mode" in str(exc_info.value)

    def test_none_raises(self):
        with pytest.raises(ValidationError) as exc_info:
            validate_mode(None)
        assert "cannot be None" in str(exc_info.value)

    def test_error_has_suggestion(self):
        with pytest.raises(ValidationError) as exc_info:
            validate_mode("INVALID")
        assert exc_info.value.suggestion is not None


class TestValidateFormula:
    """Tests for validate_formula function."""

    def test_valid_formula(self):
        assert validate_formula("C6H12O6") == "C6H12O6"
        assert validate_formula("H2O") == "H2O"
        assert validate_formula("C8H10N4O2") == "C8H10N4O2"

    def test_whitespace_normalization(self):
        assert validate_formula("C6 H12 O6") == "C6H12O6"
        assert validate_formula("  H2O  ") == "H2O"

    def test_none_raises(self):
        with pytest.raises(ValidationError):
            validate_formula(None)

    def test_empty_raises(self):
        with pytest.raises(ValidationError):
            validate_formula("")
        with pytest.raises(ValidationError):
            validate_formula("   ")

    def test_invalid_formula_raises(self):
        with pytest.raises(ValidationError) as exc_info:
            validate_formula("InvalidXYZ123!!!")
        assert "Invalid chemical formula" in str(exc_info.value)


class TestValidateFormulas:
    """Tests for validate_formulas function."""

    def test_all_valid(self):
        formulas = ["C6H12O6", "H2O", "C8H10N4O2"]
        result = validate_formulas(formulas)
        assert result == formulas

    def test_filters_invalid(self):
        formulas = ["C6H12O6", "InvalidXYZ", "H2O"]
        result = validate_formulas(formulas)
        assert "C6H12O6" in result
        assert "H2O" in result
        assert "InvalidXYZ" not in result

    def test_empty_list(self):
        assert validate_formulas([]) == []


class TestValidatePpmTolerance:
    """Tests for validate_ppm_tolerance function."""

    def test_valid_ppm(self):
        assert validate_ppm_tolerance(10.0) == 10.0
        assert validate_ppm_tolerance(5) == 5.0
        assert validate_ppm_tolerance("10") == 10.0

    def test_negative_raises(self):
        with pytest.raises(ValidationError):
            validate_ppm_tolerance(-5)

    def test_zero_raises(self):
        with pytest.raises(ValidationError):
            validate_ppm_tolerance(0)

    def test_unusually_high_raises(self):
        with pytest.raises(ValidationError) as exc_info:
            validate_ppm_tolerance(150)
        assert "unusually high" in str(exc_info.value)

    def test_invalid_type_raises(self):
        with pytest.raises(ValidationError):
            validate_ppm_tolerance("not a number")


class TestValidationError:
    """Tests for ValidationError class."""

    def test_message_only(self):
        error = ValidationError("Test message")
        assert "Test message" in str(error)
        assert error.suggestion is None

    def test_with_suggestion(self):
        error = ValidationError("Test message", suggestion="Try this")
        assert "Test message" in str(error)
        assert "Try this" in str(error)
        assert error.suggestion == "Try this"
