"""Main processing logic for LCMS Adduct Finder."""

import logging
import os
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from .config import Config
from .chemistry.mass import normalize_formula_str
from .io.raw_file import RawFileReader, sanitize_filename_component
from .io.excel import read_input_excel, write_results_excel
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


def process_all(config: Config) -> None:
    """Process all files according to configuration.

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
