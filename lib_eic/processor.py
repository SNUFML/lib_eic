"""Main processing logic for LCMS Adduct Finder."""

import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from .config import Config
from .chemistry.mass import normalize_formula_str
from .io.raw_file import RawFileReader, sanitize_filename_component
from .io.excel import (
    read_all_lc_mode_sheets,
    read_input_excel,
    read_input_excel_direct_mz,
    write_results_excel,
)
from .io.plotting import save_eic_plot
from .analysis.eic import (
    Target,
    build_targets,
    calculate_area,
    dedupe_preserve_order,
)
from .analysis.fitting import fit_gaussian_and_score, score_to_quality_label
from .analysis.ms2 import build_ms2_index, match_ms2
from .validation import validate_mode

logger = logging.getLogger(__name__)


def _normalize_mixture_value(value: Any) -> str:
    """Normalize mixture values coming from Excel into a stable string."""
    if value is None or pd.isna(value):
        return ""
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, float) and float(value).is_integer():
        return str(int(value))

    text = str(value).strip()
    if not text:
        return ""

    # Handle Excel numeric values stored as strings (e.g., "121.0")
    try:
        f = float(text)
        if np.isfinite(f) and f.is_integer():
            return str(int(f))
    except Exception:
        pass

    return text


def find_matching_raw_files(partial_filename: str, raw_folder: Path) -> List[Path]:
    """Find all .raw entries that start with the given partial filename.

    Example: "Library_POS_Mix121" matches:
      - Library_POS_Mix121.raw
      - Library_POS_Mix121_2nd.raw
      - Library_POS_Mix121_3rd.raw
    """
    raw_folder = Path(raw_folder)
    partial = str(partial_filename or "").strip()
    if partial.lower().endswith(".raw"):
        partial = partial[:-4]

    if not partial:
        return []

    if not raw_folder.exists():
        logger.warning("Raw folder not found: %s", raw_folder)
        return []

    matches: List[Path] = []
    partial_lower = partial.lower()

    try:
        for entry in raw_folder.iterdir():
            name = entry.name
            if not name.lower().endswith(".raw"):
                continue
            stem = entry.stem
            if stem.startswith(partial) or stem.lower().startswith(partial_lower):
                matches.append(entry)
    except Exception as e:
        logger.warning("Failed listing raw folder %s: %s", raw_folder, e)
        return []

    return sorted(
        matches,
        key=lambda p: (p.stem.lower() != partial_lower, p.name.lower()),
    )


def extract_file_suffix(raw_file_path: Path, partial_filename: str) -> str:
    """Extract suffix from matched raw file.

    Example:
      raw_file_path = "Library_POS_Mix121_2nd.raw"
      partial_filename = "Library_POS_Mix121"
      returns: "_2nd"
    """
    stem = Path(raw_file_path).stem
    partial = str(partial_filename or "").strip()
    if partial.lower().endswith(".raw"):
        partial = partial[:-4]

    if not partial:
        return ""

    if stem.startswith(partial):
        return stem[len(partial):]

    # Case-insensitive fallback
    if stem.lower().startswith(partial.lower()):
        return stem[len(partial):]

    return ""


def process_raw_file(
    reader: RawFileReader,
    formulas: List[str],
    mode: str,
    config: Config,
) -> List[Dict[str, Any]]:
    """Process a single raw file and extract features.

    Args:
        reader: RawFileReader instance.
        formulas: List of chemical formulas to search.
        mode: Ionization mode ("POS" or "NEG").
        config: Configuration settings.

    Returns:
        List of result dictionaries.
    """
    filename = reader.filename
    logger.info("Processing: %s (%s)", filename, mode)

    results: List[Dict[str, Any]] = []

    # Build MS2 index if available
    ms2_index = None
    ms2_enabled = False
    if reader.has_scan_events_api():
        try:
            ms2_index = build_ms2_index(reader)
            ms2_enabled = True
            logger.debug("Built MS2 index with %d entries", len(ms2_index.get("entries", [])))
        except Exception as e:
            logger.warning("MS2 index build failed: %s", e)
            ms2_index = None

    # Build targets
    targets = build_targets(
        formulas,
        mode,
        enabled_adducts=config.enabled_adducts or None,
    )
    if not targets:
        logger.warning("No valid targets for file: %s", filename)
        return results

    logger.debug("Built %d targets", len(targets))

    # Try batch extraction first
    eic_dict: Optional[Dict] = None
    use_multi = False
    if reader.has_multi_chromatogram_api():
        try:
            target_mzs = [t.mz for t in targets]
            eic_results = reader.get_chromatograms_batch(
                target_mzs, config.ppm_tolerance
            )
            eic_dict = {t.key: result for t, result in zip(targets, eic_results)}
            use_multi = True
            logger.debug("Using batch chromatogram extraction")
        except Exception as e:
            logger.debug("Batch extraction failed, falling back to single: %s", e)
            eic_dict = None
            use_multi = False

    plots_saved = 0

    for target in targets:
        formula = target.formula
        adduct_name = target.adduct
        target_mz = target.mz

        try:
            # Get EIC
            if use_multi and eic_dict is not None:
                eic_data = eic_dict.get(target.key)
                if eic_data is None:
                    logger.warning("Missing EIC for target: %s %s", formula, adduct_name)
                    continue
                eic_rt, eic_int = eic_data
            else:
                eic_rt, eic_int = reader.get_chromatogram(target_mz, config.ppm_tolerance)
        except Exception as e:
            logger.warning(
                "EIC extraction error (%s, %s, %s): %s",
                filename, formula, adduct_name, e
            )
            continue

        max_intensity = float(np.max(eic_int)) if eic_int.size else 0.0

        # Skip peaks below threshold
        if max_intensity < config.min_peak_intensity:
            continue

        total_area = calculate_area(eic_rt, eic_int, method=config.area_method)

        best_rt_min = 0.0
        gauss_score = 0.0
        quality_label = "Noise"
        rt_apex_for_ms2 = None
        fit_params = None

        if max_intensity > 1000:
            apex_idx = np.argmax(eic_int)
            best_rt_min = float(eic_rt[apex_idx])
            rt_apex_for_ms2 = best_rt_min

            if config.enable_fitting:
                gauss_score, fit_params = fit_gaussian_and_score(
                    eic_rt, eic_int, fit_rt_window_min=config.fit_rt_window_min
                )
                quality_label = score_to_quality_label(gauss_score, fitted=True)
            else:
                gauss_score = 0.0
                quality_label = score_to_quality_label(0.0, fitted=False)

            # Save plot if enabled
            if config.enable_plotting and plots_saved < config.max_plots_per_file:
                save_eic_plot(
                    rt_arr=eic_rt,
                    int_arr=eic_int,
                    formula=formula,
                    adduct=adduct_name,
                    raw_filename=filename,
                    mz_val=target_mz,
                    output_folder=config.export_plot_folder,
                    fit_params=fit_params,
                    score=gauss_score,
                    dpi=config.plot_dpi,
                )
                plots_saved += 1

        # MS2 matching
        has_ms2 = None
        ms2_match = None
        if ms2_enabled:
            has_ms2, ms2_match = match_ms2(
                ms2_index,
                target_mz,
                config.ppm_tolerance,
                rt_apex_min=rt_apex_for_ms2,
                rt_window_min=config.ms2_rt_window_min,
                mode=config.ms2_match_mode,
            )

        # Build result row
        row_out: Dict[str, Any] = {
            "RawFile": filename,
            "Mode": mode,
            "Formula": formula,
            "Adduct": adduct_name,
            "mz_theoretical": target_mz,
            "RT_min": round(best_rt_min, 3),
            "Intensity": max_intensity,
            "Area": total_area,
            "GaussianScore": round(gauss_score, 3),
            "PeakQuality": quality_label,
            "HasMS2": bool(has_ms2) if ms2_enabled else None,
        }

        if ms2_enabled and config.store_ms2_match_details:
            row_out["MS2ScanNo"] = ms2_match.get("scan_no") if ms2_match else None
            row_out["MS2RT_min"] = ms2_match.get("rt_min") if ms2_match else None
            row_out["MS2Precursor_mz"] = ms2_match.get("precursor_mz") if ms2_match else None

        results.append(row_out)

    logger.info("Found %d features in %s", len(results), filename)
    return results


def process_raw_file_direct_mz(
    reader: RawFileReader,
    targets_df: pd.DataFrame,
    lc_mode: str,
    polarity: str,
    partial_filename: str,
    file_suffix: str,
    config: Config,
) -> List[Dict[str, Any]]:
    """Process a single raw file using direct m/z targets from input Excel."""
    filename = reader.filename
    logger.info("Processing: %s (%s)", filename, polarity)

    results: List[Dict[str, Any]] = []

    # Build MS2 index if available
    ms2_index = None
    ms2_enabled = False
    if reader.has_scan_events_api():
        try:
            ms2_index = build_ms2_index(reader)
            ms2_enabled = True
            logger.debug(
                "Built MS2 index with %d entries", len(ms2_index.get("entries", []))
            )
        except Exception as e:
            logger.warning("MS2 index build failed: %s", e)
            ms2_index = None

    # Ensure m/z is numeric and filtered
    targets_df = targets_df.copy()
    targets_df["m/z"] = pd.to_numeric(targets_df["m/z"], errors="coerce")
    targets_df = targets_df[targets_df["m/z"].notna() & (targets_df["m/z"] != 0)]
    targets_df = targets_df.reset_index(drop=True)

    if targets_df.empty:
        return results

    mz_list = targets_df["m/z"].astype(float).tolist()

    eic_results: Optional[List] = None
    use_multi = False
    if reader.has_multi_chromatogram_api():
        try:
            eic_results = reader.get_chromatograms_batch(mz_list, config.ppm_tolerance)
            use_multi = True
            logger.debug("Using batch chromatogram extraction")
        except Exception as e:
            logger.debug("Batch extraction failed, falling back to single: %s", e)
            eic_results = None
            use_multi = False

    plots_saved = 0

    for i, row in targets_df.iterrows():
        mz_val = float(row["m/z"])

        compound_raw = row.get("Compound name")
        compound_name = "" if pd.isna(compound_raw) else str(compound_raw).strip()
        if not compound_name:
            compound_name = "Unknown"

        mixture = _normalize_mixture_value(row.get("mixture"))

        try:
            if use_multi and eic_results is not None:
                eic_rt, eic_int = eic_results[i]
            else:
                eic_rt, eic_int = reader.get_chromatogram(mz_val, config.ppm_tolerance)
        except Exception as e:
            logger.warning(
                "EIC extraction error (%s, %s, %.4f): %s",
                filename,
                compound_name,
                mz_val,
                e,
            )
            continue

        max_intensity = float(np.max(eic_int)) if eic_int.size else 0.0
        if max_intensity < config.min_peak_intensity:
            continue

        total_area = calculate_area(eic_rt, eic_int, method=config.area_method)

        best_rt_min = 0.0
        gauss_score = 0.0
        quality_label = "Noise"
        rt_apex_for_ms2 = None
        fit_params = None

        if max_intensity > 1000:
            apex_idx = int(np.argmax(eic_int))
            best_rt_min = float(eic_rt[apex_idx]) if eic_rt.size else 0.0
            rt_apex_for_ms2 = best_rt_min

            if config.enable_fitting:
                gauss_score, fit_params = fit_gaussian_and_score(
                    eic_rt, eic_int, fit_rt_window_min=config.fit_rt_window_min
                )
                quality_label = score_to_quality_label(gauss_score, fitted=True)
            else:
                gauss_score = 0.0
                quality_label = score_to_quality_label(0.0, fitted=False)

            if config.enable_plotting and plots_saved < config.max_plots_per_file:
                from .io.plotting import save_eic_plot_direct_mz

                save_eic_plot_direct_mz(
                    rt_arr=eic_rt,
                    int_arr=eic_int,
                    compound_name=compound_name,
                    polarity=polarity,
                    raw_filename=filename,
                    mz_val=mz_val,
                    mixture=mixture,
                    output_folder=config.export_plot_folder,
                    lc_mode=lc_mode,
                    file_suffix=file_suffix,
                    partial_filename=partial_filename,
                    fit_params=fit_params,
                    score=gauss_score,
                    dpi=config.plot_dpi,
                )
                plots_saved += 1

        # MS2 matching
        has_ms2 = None
        ms2_match = None
        if ms2_enabled:
            has_ms2, ms2_match = match_ms2(
                ms2_index,
                mz_val,
                config.ppm_tolerance,
                rt_apex_min=rt_apex_for_ms2,
                rt_window_min=config.ms2_rt_window_min,
                mode=config.ms2_match_mode,
            )

        row_out: Dict[str, Any] = {
            "RawFile": filename,
            "File name": str(partial_filename),
            "lc_mode": str(lc_mode),
            "mixture": mixture,
            "Compound name": compound_name,
            "Polarity": polarity,
            "mz_target": mz_val,
            "RT_min": round(best_rt_min, 3),
            "Intensity": max_intensity,
            "Area": total_area,
            "GaussianScore": round(gauss_score, 3),
            "PeakQuality": quality_label,
            "HasMS2": bool(has_ms2) if ms2_enabled else None,
        }

        if ms2_enabled and config.store_ms2_match_details:
            row_out["MS2ScanNo"] = ms2_match.get("scan_no") if ms2_match else None
            row_out["MS2RT_min"] = ms2_match.get("rt_min") if ms2_match else None
            row_out["MS2Precursor_mz"] = (
                ms2_match.get("precursor_mz") if ms2_match else None
            )

        results.append(row_out)

    logger.info("Found %d features in %s", len(results), filename)
    return results


def process_all_formula_based(config: Config) -> None:
    """Process all files using formula + adduct based target generation.

    Args:
        config: Configuration settings.

    Raises:
        FileNotFoundError: If input file doesn't exist.
        ValueError: If input file format is invalid.
    """
    # Read input Excel file
    if not os.path.exists(config.input_excel):
        raise FileNotFoundError(
            f"Input Excel file not found: {config.input_excel}"
        )

    logger.info("Reading input file: %s", config.input_excel)
    meta_data = read_input_excel(
        config.input_excel,
        sheet_name=config.input_sheet,
        required_columns={"RawFile", "Mode", "Formula"},
    )

    all_results: List[Dict[str, Any]] = []

    # Group by RawFile and Mode
    grouped = meta_data.groupby(["RawFile", "Mode"])
    total_files = len(grouped)

    logger.info("Processing %d file groups", total_files)

    for idx, ((raw_filename, mode), group_df) in enumerate(grouped):
        raw_filename = str(raw_filename).strip()
        raw_filename = sanitize_filename_component(raw_filename)

        # Validate mode
        try:
            mode = validate_mode(mode)
        except ValueError as e:
            logger.error("Invalid mode for file %s: %s", raw_filename, e)
            continue

        # Check file extension
        raw_filename_lower = raw_filename.lower()
        if raw_filename_lower.endswith(".mzml"):
            logger.error(
                "[%d/%d] File: %s - mzML files are not supported. Use Thermo .raw files.",
                idx + 1, total_files, raw_filename
            )
            continue

        if not raw_filename_lower.endswith(".raw"):
            raw_filename = raw_filename + ".raw"

        # Build full path
        full_file_path = os.path.join(str(config.raw_data_folder), raw_filename)

        # Parse formulas
        formulas_raw = group_df["Formula"].dropna().astype(str).tolist()
        formulas_norm = [normalize_formula_str(f) for f in formulas_raw]
        formulas_norm = [f for f in formulas_norm if f]
        formulas = dedupe_preserve_order(formulas_norm)

        logger.info("[%d/%d] File: %s", idx + 1, total_files, raw_filename)

        if not os.path.exists(full_file_path):
            logger.error("File not found: %s", full_file_path)
            continue

        try:
            with RawFileReader(full_file_path) as reader:
                file_results = process_raw_file(reader, formulas, mode, config)
                all_results.extend(file_results)
        except Exception as e:
            logger.error("Failed to process file %s: %s", raw_filename, e)
            continue

    # Save results
    if all_results:
        logger.info("Saving %d results to: %s", len(all_results), config.output_excel)
        write_results_excel(all_results, config.output_excel, include_pivot_tables=True)
        logger.info("Processing complete: %s", config.output_excel)
    else:
        logger.warning("No results to save")


def process_all_direct_mz(config: Config) -> None:
    """Process all matching raw files using direct m/z targets from input Excel."""
    if not os.path.exists(config.input_excel):
        raise FileNotFoundError(f"Input Excel file not found: {config.input_excel}")

    lc_mode_data = read_all_lc_mode_sheets(
        config.input_excel,
        lc_modes=config.input_sheets,
    )

    # Backward-compatible fallback: attempt single configured sheet if no LC-mode
    # sheets were readable (e.g., legacy "Final" direct-m/z format).
    if not lc_mode_data:
        meta_data = read_input_excel_direct_mz(
            config.input_excel,
            sheet_name=config.input_sheet,
            lc_mode=config.input_sheet,
        )
        lc_mode_data = {str(config.input_sheet): meta_data}

    all_results: List[Dict[str, Any]] = []

    total_groups = sum(len(df.groupby(["File name", "Polarity"])) for df in lc_mode_data.values())

    logger.info("Processing %d file groups (direct m/z)", total_groups)

    group_counter = 0

    for lc_mode, meta_data in lc_mode_data.items():
        logger.info("Processing LC mode: %s (%d rows)", lc_mode, len(meta_data))

        grouped = meta_data.groupby(["File name", "Polarity"])

        for (partial_filename, polarity), group_df in grouped:
            group_counter += 1
            partial_filename = str(partial_filename).strip()
            if not partial_filename:
                logger.warning("Skipping empty 'File name' at group %d", group_counter)
                continue

            try:
                polarity_norm = validate_mode(polarity)
            except ValueError as e:
                logger.error(
                    "Invalid polarity for file %s (%s): %s", partial_filename, lc_mode, e
                )
                continue

            matching_files = find_matching_raw_files(
                partial_filename, config.raw_data_folder
            )
            if not matching_files:
                logger.warning(
                    "No matching raw files for: %s (%s)", partial_filename, lc_mode
                )
                continue

            for raw_file_path in matching_files:
                file_suffix = extract_file_suffix(raw_file_path, partial_filename)

                logger.info(
                    "[%d/%d] File: %s | %s (matched: %s)",
                    group_counter,
                    total_groups,
                    lc_mode,
                    partial_filename,
                    raw_file_path.name,
                )

                try:
                    with RawFileReader(str(raw_file_path)) as reader:
                        file_results = process_raw_file_direct_mz(
                            reader=reader,
                            targets_df=group_df,
                            lc_mode=lc_mode,
                            polarity=polarity_norm,
                            partial_filename=partial_filename,
                            file_suffix=file_suffix,
                            config=config,
                        )
                        all_results.extend(file_results)
                except Exception as e:
                    logger.error("Failed to process file %s: %s", raw_file_path, e)
                    continue

    if all_results:
        logger.info("Saving %d results to: %s", len(all_results), config.output_excel)
        write_results_excel(all_results, config.output_excel, include_pivot_tables=True)
        logger.info("Processing complete: %s", config.output_excel)
    else:
        logger.warning("No results to save")


def process_all(config: Config) -> None:
    """Process all files according to configuration.

    Automatically detects the input Excel format:
      - Direct m/z format: columns include "File name", "Compound name", "Polarity", "m/z"
      - Formula-based format: columns include "RawFile", "Mode", "Formula"
    """
    try:
        process_all_direct_mz(config)
        return
    except ValueError as e:
        # Only fall back if the direct-m/z reader couldn't find required columns.
        msg = str(e)
        if "Missing required columns" not in msg:
            raise

    process_all_formula_based(config)
