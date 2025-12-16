"""Excel I/O operations for LCMS Adduct Finder."""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd

logger = logging.getLogger(__name__)


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

    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Input Excel file not found: {file_path}")

    logger.info("Reading Excel file: %s (sheet: %s)", file_path, sheet_name)

    try:
        df = pd.read_excel(file_path, sheet_name=sheet_name)
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


def write_results_excel(
    results: List[Dict],
    output_path: str,
    include_pivot_tables: bool = True,
) -> None:
    """Write analysis results to Excel file.

    Args:
        results: List of result dictionaries.
        output_path: Path to output Excel file.
        include_pivot_tables: Whether to include per-formula pivot tables.
    """
    if not results:
        logger.warning("No results to save")
        return

    logger.info("Saving results to: %s", output_path)

    df_results = pd.DataFrame(results)

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        # Main results sheet
        df_results.to_excel(writer, sheet_name="All_Features", index=False)
        logger.debug("Wrote %d rows to All_Features sheet", len(df_results))

        if not include_pivot_tables:
            return

        # Create per-formula sheets with pivot tables
        unique_formulas = df_results["Formula"].unique()

        for formula in unique_formulas:
            f_data = df_results[df_results["Formula"] == formula]

            # Create pivot tables
            pivot_area = f_data.pivot_table(
                index="RawFile", columns="Adduct", values="Area"
            )
            pivot_rt = f_data.pivot_table(
                index="RawFile", columns="Adduct", values="RT_min"
            )
            pivot_intensity = f_data.pivot_table(
                index="RawFile", columns="Adduct", values="Intensity"
            )

            # Generate safe sheet name (max 31 chars for Excel)
            safe_name = "".join(c for c in formula if c.isalnum())[:30]

            # Write Area Table
            pivot_area.to_excel(writer, sheet_name=safe_name, startrow=0)
            writer.sheets[safe_name].cell(row=1, column=1).value = "Area Table"

            # Write Retention Time Table
            current_row = len(pivot_area) + 3
            writer.sheets[safe_name].cell(row=current_row, column=1).value = (
                "Retention Time (min)"
            )
            pivot_rt.to_excel(writer, sheet_name=safe_name, startrow=current_row + 1)

            # Write Peak Height (Intensity) Table
            current_row = len(pivot_area) + len(pivot_rt) + 6
            writer.sheets[safe_name].cell(row=current_row, column=1).value = (
                "Peak Height (Intensity)"
            )
            pivot_intensity.to_excel(
                writer, sheet_name=safe_name, startrow=current_row + 1
            )

            logger.debug("Wrote pivot tables for formula: %s", formula)

    logger.info("Results saved successfully")
