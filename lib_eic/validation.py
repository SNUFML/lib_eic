"""Input validation for LCMS Adduct Finder."""

from pathlib import Path
from typing import List, Optional

from .chemistry.mass import normalize_formula_str


class ValidationError(ValueError):
    """Exception raised for validation errors with helpful suggestions."""

    def __init__(self, message: str, suggestion: Optional[str] = None):
        self.suggestion = suggestion
        full_message = message
        if suggestion:
            full_message = f"{message}\nSuggestion: {suggestion}"
        super().__init__(full_message)


def validate_mode(mode: str) -> str:
    """Validate and normalize ionization mode.

    Args:
        mode: Ionization mode string.

    Returns:
        Normalized mode ("POS" or "NEG").

    Raises:
        ValidationError: If mode is invalid.
    """
    if mode is None:
        raise ValidationError(
            "Mode cannot be None",
            suggestion="Use 'POS' for positive mode or 'NEG' for negative mode.",
        )

    mode_normalized = str(mode).upper().strip()

    if mode_normalized in ("POS", "POSITIVE", "+"):
        return "POS"
    elif mode_normalized in ("NEG", "NEGATIVE", "-"):
        return "NEG"
    else:
        raise ValidationError(
            f"Invalid ionization mode: {mode!r}",
            suggestion="Use 'POS' for positive mode or 'NEG' for negative mode.",
        )


def validate_formula(formula_str: str) -> str:
    """Validate and normalize chemical formula.

    Args:
        formula_str: Chemical formula string.

    Returns:
        Normalized formula string.

    Raises:
        ValidationError: If formula is invalid.
    """
    if formula_str is None:
        raise ValidationError(
            "Formula cannot be None",
            suggestion="Provide a valid chemical formula like 'C6H12O6'.",
        )

    formula_normalized = normalize_formula_str(formula_str)

    if not formula_normalized:
        raise ValidationError(
            f"Empty formula after normalization: {formula_str!r}",
            suggestion="Provide a valid chemical formula like 'C6H12O6'.",
        )

    # Try to parse with molmass to validate
    try:
        from molmass import Formula

        f = Formula(formula_normalized)
        # Access mass property to trigger actual validation
        _ = f.isotope.mass
    except ImportError:
        # Can't validate without molmass, return as-is
        return formula_normalized
    except Exception as e:
        raise ValidationError(
            f"Invalid chemical formula: {formula_str!r}",
            suggestion=f"Check formula syntax. Error: {e}",
        )

    return formula_normalized


def validate_formulas(formulas: List[str]) -> List[str]:
    """Validate and normalize a list of formulas.

    Invalid formulas are logged and skipped.

    Args:
        formulas: List of formula strings.

    Returns:
        List of valid, normalized formulas.
    """
    import logging

    logger = logging.getLogger(__name__)

    valid_formulas = []
    for formula in formulas:
        try:
            valid_formulas.append(validate_formula(formula))
        except ValidationError as e:
            logger.warning("Skipping invalid formula %r: %s", formula, e)

    return valid_formulas


def validate_raw_file_path(file_path: str) -> Path:
    """Validate raw file path.

    Args:
        file_path: Path to raw file.

    Returns:
        Path object for the file.

    Raises:
        ValidationError: If file doesn't exist or has wrong extension.
    """
    path = Path(file_path)

    if not path.exists():
        raise ValidationError(
            f"Raw file not found: {file_path}",
            suggestion="Check that the file path is correct and the file exists.",
        )

    # Check extension
    suffix = path.suffix.lower()
    if suffix == ".mzml":
        raise ValidationError(
            f"mzML files are not supported: {file_path}",
            suggestion="Use Thermo .raw files instead.",
        )

    if suffix != ".raw":
        raise ValidationError(
            f"Unsupported file format: {suffix}",
            suggestion="Only Thermo .raw files are supported.",
        )

    return path


def validate_excel_file(file_path: str) -> Path:
    """Validate Excel file path.

    Args:
        file_path: Path to Excel file.

    Returns:
        Path object for the file.

    Raises:
        ValidationError: If file doesn't exist or has wrong extension.
    """
    path = Path(file_path)

    if not path.exists():
        raise ValidationError(
            f"Excel file not found: {file_path}",
            suggestion="Check that the file path is correct.",
        )

    suffix = path.suffix.lower()
    if suffix not in (".xlsx", ".xls"):
        raise ValidationError(
            f"Invalid Excel file format: {suffix}",
            suggestion="Use .xlsx or .xls Excel files.",
        )

    return path


def validate_ppm_tolerance(ppm: float) -> float:
    """Validate ppm tolerance value.

    Args:
        ppm: PPM tolerance value.

    Returns:
        Validated ppm value.

    Raises:
        ValidationError: If ppm is invalid.
    """
    try:
        ppm = float(ppm)
    except (TypeError, ValueError):
        raise ValidationError(
            f"PPM tolerance must be a number: {ppm!r}",
            suggestion="Use a numeric value like 10.0 or 5.0.",
        )

    if ppm <= 0:
        raise ValidationError(
            f"PPM tolerance must be positive: {ppm}",
            suggestion="Use a positive value like 10.0.",
        )

    if ppm > 100:
        raise ValidationError(
            f"PPM tolerance seems unusually high: {ppm}",
            suggestion="Typical values are 5-20 ppm for high-resolution MS.",
        )

    return ppm


def validate_config(config) -> None:
    """Validate configuration object.

    Args:
        config: Config instance.

    Raises:
        ValidationError: If config has invalid values.
    """
    # Validate ppm_tolerance
    validate_ppm_tolerance(config.ppm_tolerance)

    # Validate min_peak_intensity
    if config.min_peak_intensity < 0:
        raise ValidationError(
            f"min_peak_intensity cannot be negative: {config.min_peak_intensity}",
            suggestion="Use 0 or a positive value.",
        )

    # Validate ms2_match_mode
    if config.ms2_match_mode not in ("rt_linked", "global"):
        raise ValidationError(
            f"Invalid ms2_match_mode: {config.ms2_match_mode!r}",
            suggestion="Use 'rt_linked' or 'global'.",
        )

    # Validate area_method
    if config.area_method not in ("sum", "trapz"):
        raise ValidationError(
            f"Invalid area_method: {config.area_method!r}",
            suggestion="Use 'sum' or 'trapz'.",
        )

    # Validate raw_data_folder if specified
    if config.raw_data_folder:
        if not config.raw_data_folder.exists():
            raise ValidationError(
                f"Raw data folder not found: {config.raw_data_folder}",
                suggestion="Create the folder or update the path.",
            )
