"""Excel I/O operations for LCMS Adduct Finder."""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import pandas as pd

logger = logging.getLogger(__name__)

EXCEL_SHEETNAME_MAX_LEN = 31


def _sanitize_sheet_name(value: Any, *, max_len: int = 30) -> str:
    """Convert an arbitrary label into a safe Excel sheet name.

    Args:
        value: Label to convert (e.g., a formula, compound name, etc.).
        max_len: Maximum length to keep before uniqueness suffixing.

    Returns:
        A sanitized sheet name (never empty).

    Notes:
        Excel sheet names are limited to 31 characters. We keep the historical
        30-character truncation as the default to avoid changing existing outputs.
    """
    safe = "".join(c for c in str(value) if c.isalnum())[:max_len]
    return safe or "Target"


def _make_unique_sheet_name(base: str, used_names: Set[str]) -> str:
    """Ensure ``base`` is unique within ``used_names`` (and <=31 chars).

    Args:
        base: Base sheet name (should already be sanitized).
        used_names: Set of sheet names already present; mutated in-place.

    Returns:
        A unique sheet name within the Excel 31-character limit.
    """
    base = str(base) or "Target"

    if base not in used_names:
        used_names.add(base)
        return base

    counter = 1
    while True:
        suffix = f"_{counter}"
        trimmed = base[: EXCEL_SHEETNAME_MAX_LEN - len(suffix)]
        candidate = f"{trimmed}{suffix}"
        if candidate not in used_names:
            used_names.add(candidate)
            return candidate
        counter += 1


def read_input_excel(
    file_path: str,
    sheet_name: str = "Final",
    required_columns: Optional[Set[str]] = None,
) -> pd.DataFrame:
    """Read input Excel file with compound information.

    Args:
        file_path: Path to the Excel file.
        sheet_name: Name of the sheet to read.
        required_columns: Set of required column names.

    Returns:
        DataFrame with compound information.

    Raises:
        FileNotFoundError: If Excel file doesn't exist.
        ValueError: If required columns are missing.
    """
    if required_columns is None:
        required_columns = {"RawFile", "Mode", "Formula"}

    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Input Excel file not found: {path}")

    logger.info("Reading Excel file: %s (sheet: %s)", path, sheet_name)

    try:
        df = pd.read_excel(path, sheet_name=sheet_name)
    except Exception as e:
        raise ValueError(f"Failed to read Excel file: {e}") from e

    # Check required columns
    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns in sheet '{sheet_name}': {sorted(missing)}. "
            f"Available columns: {list(df.columns)}"
        )

    logger.info("Read %d rows from Excel file", len(df))
    return df


def read_input_excel_direct_mz(
    file_path: str,
    sheet_name: str = "Final",
    lc_mode: Optional[str] = None,
) -> pd.DataFrame:
    """Read input Excel file with direct m/z values.

    Expected columns (after stripping whitespace from headers):
        - "num" (optional; used for plot filename ordering)
        - "File name"
        - "mixture"
        - "Compound name"
        - "Polarity"
        - "m/z"

    Notes:
        The input Excel may contain merged cells in the first row; actual
        headers start from row 2, so we read with ``skiprows=1``.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"Input Excel file not found: {path}")

    logger.info("Reading Excel file (direct m/z): %s (sheet: %s)", path, sheet_name)

    try:
        df = pd.read_excel(path, sheet_name=sheet_name, skiprows=1)
    except Exception as e:
        raise ValueError(f"Failed to read Excel file: {e}") from e

    # Normalize column labels (common with Excel exports)
    df.columns = [str(c).strip() for c in df.columns]

    # Normalize optional numbering column
    for col in list(df.columns):
        if str(col).strip().lower() == "num" and col != "num":
            df = df.rename(columns={col: "num"})
            break

    required = {"File name", "mixture", "Compound name", "Polarity", "m/z"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns in sheet '{sheet_name}': {sorted(missing)}. "
            f"Available columns: {list(df.columns)}"
        )

    # Coerce m/z to numeric and filter invalid rows
    df["m/z"] = pd.to_numeric(df["m/z"], errors="coerce")
    df = df[df["m/z"].notna() & (df["m/z"] != 0)]

    # Keep only required columns in a stable order
    ordered_cols = ["File name", "mixture", "Compound name", "Polarity", "m/z"]
    if "num" in df.columns:
        ordered_cols = ["num"] + ordered_cols
    df = df[ordered_cols].copy()

    if lc_mode is not None:
        lc_mode_text = str(lc_mode).strip()
        if lc_mode_text:
            if "lc_mode" in df.columns:
                df["lc_mode"] = lc_mode_text
            else:
                df.insert(0, "lc_mode", lc_mode_text)

    logger.info("Read %d rows from Excel file (direct m/z)", len(df))
    return df


def read_all_lc_mode_sheets(
    file_path: str,
    lc_modes: Optional[List[str]] = None,
) -> Dict[str, pd.DataFrame]:
    """Read input Excel with separate LC mode sheets (e.g., RP and HILIC).

    Args:
        file_path: Path to the Excel file.
        lc_modes: List of sheet names to read. Defaults to ["RP", "HILIC"].

    Returns:
        Dict mapping lc_mode -> DataFrame with an added "lc_mode" column.
    """
    if lc_modes is None:
        lc_modes = ["RP", "HILIC"]

    result: Dict[str, pd.DataFrame] = {}
    for lc_mode in lc_modes:
        lc_mode_name = str(lc_mode).strip()
        if not lc_mode_name:
            continue

        try:
            df = read_input_excel_direct_mz(
                file_path,
                sheet_name=lc_mode_name,
                lc_mode=lc_mode_name,
            )
            result[lc_mode_name] = df
            logger.info("Read %d rows from %s sheet", len(df), lc_mode_name)
        except Exception as e:
            logger.warning("Could not read sheet '%s': %s", lc_mode_name, e)

    return result


def write_results_excel(
    results: List[Dict],
    output_path: str,
    include_pivot_tables: bool = False,
    *,
    status_rows: Optional[List[Dict]] = None,
) -> None:
    """Write analysis results to Excel file.

    Args:
        results: List of result dictionaries.
        output_path: Path to output Excel file.
        include_pivot_tables: Whether to include per-target pivot tables.
        status_rows: Optional list of per-target status rows (includes filtered and
            failed extractions) to be appended to the All_Features sheet.
    """
    if not results and not status_rows:
        logger.warning("No results to save")
        return

    logger.info("Saving results to: %s", output_path)

    df_results = pd.DataFrame(results) if results else pd.DataFrame()
    df_status = pd.DataFrame(status_rows) if status_rows else pd.DataFrame()

    df_all = df_status if not df_status.empty else df_results

    def _maybe_normalize_num_column(df: pd.DataFrame) -> None:
        if df.empty or "num" not in df.columns:
            return

        num_series = df["num"]
        num_numeric = pd.to_numeric(num_series, errors="coerce")
        convertible = pd.isna(num_series) | num_numeric.notna()
        if not bool(convertible.all()):
            return

        num_non_na = num_numeric.dropna()
        if not num_non_na.empty and bool((num_non_na % 1 == 0).all()):
            df["num"] = num_numeric.astype("Int64")
        else:
            df["num"] = num_numeric

    def _maybe_move_num_first(df: pd.DataFrame) -> pd.DataFrame:
        if df.empty or "num" not in df.columns:
            return df
        cols = ["num"] + [c for c in df.columns if c != "num"]
        return df[cols]

    _maybe_normalize_num_column(df_results)
    df_results = _maybe_move_num_first(df_results)
    _maybe_normalize_num_column(df_all)
    df_all = _maybe_move_num_first(df_all)

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        # Main results sheet
        df_all.to_excel(writer, sheet_name="All_Features", index=False)
        logger.debug("Wrote %d rows to All_Features sheet", len(df_all))

        if not include_pivot_tables or df_results.empty:
            return

        # Create per-target sheets with pivot tables
        if "Formula" in df_results.columns and "Adduct" in df_results.columns:
            target_col = "Formula"
            column_col = "Adduct"
        elif "Compound name" in df_results.columns and "Polarity" in df_results.columns:
            target_col = "Compound name"
            column_col = "Polarity"
        else:
            logger.warning(
                "Pivot tables skipped: required columns not found. "
                "Expected either (Formula, Adduct) or (Compound name, Polarity). "
                "Available columns: %s",
                list(df_results.columns),
            )
            return

        unique_targets = df_results[target_col].dropna().unique()
        used_sheet_names: Set[str] = set(writer.sheets.keys())

        for target in unique_targets:
            f_data = df_results[df_results[target_col] == target]

            # Create pivot tables
            pivot_area = f_data.pivot_table(index="RawFile", columns=column_col, values="Area")
            pivot_rt = f_data.pivot_table(index="RawFile", columns=column_col, values="RT_min")
            pivot_intensity = f_data.pivot_table(index="RawFile", columns=column_col, values="Intensity")

            base_name = _sanitize_sheet_name(target, max_len=30)
            safe_name = _make_unique_sheet_name(base_name, used_sheet_names)

            # Write Area Table
            pivot_area.to_excel(writer, sheet_name=safe_name, startrow=0)
            writer.sheets[safe_name].cell(row=1, column=1).value = "Area Table"

            # Write Retention Time Table
            current_row = len(pivot_area) + 3
            writer.sheets[safe_name].cell(row=current_row, column=1).value = "Retention Time (min)"
            pivot_rt.to_excel(writer, sheet_name=safe_name, startrow=current_row + 1)

            # Write Peak Height (Intensity) Table
            current_row = len(pivot_area) + len(pivot_rt) + 6
            writer.sheets[safe_name].cell(row=current_row, column=1).value = "Peak Height (Intensity)"
            pivot_intensity.to_excel(writer, sheet_name=safe_name, startrow=current_row + 1)

            logger.debug("Wrote pivot tables for %s: %s", target_col, target)

    logger.info("Results saved successfully")
