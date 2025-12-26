"""Main processing logic for LCMS Adduct Finder."""

from bisect import bisect_left, bisect_right
from concurrent.futures import as_completed
import logging
import os
import re
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from .config import Config
from .chemistry.mass import normalize_formula_str
from .io.raw_file import RawFileReader
from .io.excel import (
    read_all_lc_mode_sheets,
    read_input_excel,
    read_input_excel_direct_mz,
    write_results_excel,
)
from .io.plotting import create_eic_plotter, save_eic_plot
from .analysis.eic import (
    build_targets,
    calculate_area,
    dedupe_preserve_order,
)
from .analysis.fitting import fit_gaussian_and_score, score_to_quality_label
from .analysis.ms2 import build_ms2_index, match_ms2
from .parallel import create_process_pool, resolve_max_workers, should_use_process_pool
from .progress import progress_bar, should_show_progress
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


def _normalize_num_value(value: Any) -> str:
    """Normalize numbering values (e.g., 'num' column) into a stable string."""
    if value is None or pd.isna(value):
        return ""
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, float) and float(value).is_integer():
        return str(int(value))

    text = str(value).strip()
    if not text:
        return ""

    try:
        f = float(text)
        if np.isfinite(f) and f.is_integer():
            return str(int(f))
    except Exception:
        pass

    return text


def _collect_raw_entries_recursive(raw_folder: Path) -> List[Path]:
    """Collect ``.raw`` file/directory entries under ``raw_folder`` (recursive).

    Notes:
        Thermo ``.raw`` can be either a file (common on Linux exports) or a
        directory (common on Windows). If a ``.raw`` directory is encountered,
        we treat it as a leaf and **do not** recurse into it.
    """
    raw_folder = Path(raw_folder)
    if not raw_folder.exists():
        return []

    matches: List[Path] = []
    try:
        for root, dirs, files in os.walk(raw_folder):
            # Treat "*.raw" directories as leaf nodes and do not recurse into them.
            raw_dirs = [d for d in dirs if str(d).lower().endswith(".raw")]
            for d in raw_dirs:
                matches.append(Path(root) / d)
            dirs[:] = [d for d in dirs if d not in raw_dirs]

            for name in files:
                if str(name).lower().endswith(".raw"):
                    matches.append(Path(root) / name)
    except Exception as e:
        logger.warning("Failed walking raw folder %s: %s", raw_folder, e)
        return []

    return sorted(matches, key=lambda p: str(p).lower())


def _build_raw_entry_name_index(entries: List[Path]) -> Dict[str, List[Path]]:
    """Build a case-insensitive index of ``Path.name`` -> matching entries."""
    index: Dict[str, List[Path]] = {}
    for entry in entries:
        index.setdefault(entry.name.lower(), []).append(entry)
    return index


def _build_raw_entry_stem_index(entries: List[Path]) -> Tuple[List[str], List[Path], Dict[str, List[Path]]]:
    """Build a prefix-searchable index for matching ``Path.stem`` values.

    Returns:
        Tuple of (sorted_stems, sorted_paths, exact_map).
    """
    items: List[Tuple[str, Path]] = []
    exact: Dict[str, List[Path]] = {}

    for entry in entries:
        if not entry.name.lower().endswith(".raw"):
            continue
        stem_lower = entry.stem.lower()
        items.append((stem_lower, entry))
        exact.setdefault(stem_lower, []).append(entry)

    items.sort(key=lambda t: t[0])
    stems = [stem for stem, _ in items]
    paths = [path for _, path in items]
    return stems, paths, exact


def _resolve_lc_mode_raw_folder(raw_root: Path, lc_mode: str) -> Path:
    """Resolve an LC-mode-specific raw folder if present.

    If ``raw_root/<lc_mode>`` exists (case-insensitive), returns it; otherwise
    returns ``raw_root``.
    """
    raw_root = Path(raw_root)
    lc_mode_text = str(lc_mode or "").strip()
    if not lc_mode_text:
        return raw_root

    direct = raw_root / lc_mode_text
    if direct.exists():
        return direct

    try:
        lc_mode_lower = lc_mode_text.lower()
        for entry in raw_root.iterdir():
            if entry.is_dir() and entry.name.lower() == lc_mode_lower:
                return entry
    except Exception:
        return raw_root

    return raw_root


def _infer_run_label(raw_file_path: Path, search_root: Path) -> str:
    """Infer run label (e.g., '1st', '2nd') from folder structure."""
    try:
        rel = raw_file_path.relative_to(search_root)
    except Exception:
        return ""

    if len(rel.parts) < 2:
        return ""

    candidate = str(rel.parts[0]).strip()
    if not candidate:
        return ""

    if re.fullmatch(r"\d+(st|nd|rd|th)", candidate.lower()):
        return candidate

    return ""


def _build_raw_file_id(raw_file_path: Path, raw_root: Path) -> str:
    """Build a stable RawFile identifier for output tables."""
    try:
        return raw_file_path.relative_to(raw_root).as_posix()
    except Exception:
        return raw_file_path.name


def find_matching_raw_files(
    partial_filename: str,
    raw_folder: Path,
    *,
    raw_entries: Optional[List[Path]] = None,
    raw_entry_stem_index: Optional[Tuple[List[str], List[Path], Dict[str, List[Path]]]] = None,
) -> List[Path]:
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

    entries = raw_entries
    if entries is None and raw_entry_stem_index is None:
        if not raw_folder.exists():
            logger.warning("Raw folder not found: %s", raw_folder)
            return []

        try:
            entries = [p for p in raw_folder.iterdir()]
        except Exception as e:
            logger.warning("Failed listing raw folder %s: %s", raw_folder, e)
            return []

    matches: List[Path] = []
    partial_lower = partial.lower()

    if raw_entry_stem_index is not None:
        stems, paths, exact = raw_entry_stem_index
        matches.extend(exact.get(partial_lower, []))

        prefix = f"{partial_lower}_"
        lo = bisect_left(stems, prefix)
        hi = bisect_right(stems, f"{prefix}\uffff")
        matches.extend(paths[lo:hi])

        return sorted(
            matches,
            key=lambda p: (p.stem.lower() != partial_lower, p.name.lower(), str(p).lower()),
        )

    if entries is None:
        return []

    for entry in entries:
        name = entry.name
        if not name.lower().endswith(".raw"):
            continue
        stem = entry.stem
        stem_lower = stem.lower()
        if stem_lower == partial_lower:
            matches.append(entry)
            continue
        if stem_lower.startswith(partial_lower):
            # Only treat as a "repeat" match if the suffix is underscore-delimited
            # (e.g., prevents "Mixture_1" matching "Mixture_10").
            next_ch = stem_lower[len(partial_lower) : len(partial_lower) + 1]
            if next_ch == "_":
                matches.append(entry)

    return sorted(
        matches,
        key=lambda p: (p.stem.lower() != partial_lower, p.name.lower(), str(p).lower()),
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
        return stem[len(partial) :]

    # Case-insensitive fallback
    if stem.lower().startswith(partial.lower()):
        return stem[len(partial) :]

    return ""


def _status_rows_for_formula_file_failure(
    *,
    raw_file_id: str,
    mode: str,
    formulas: List[str],
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for formula in formulas:
        rows.append(
            {
                "RawFile": raw_file_id,
                "Mode": mode,
                "Formula": formula,
                "Adduct": None,
                "mz_theoretical": None,
                "RT_min": None,
                "Intensity": None,
                "Area": None,
                "GaussianScore": None,
                "PeakQuality": None,
                "HasMS2": None,
                "EICGenerated": False,
                "FilteredOut": False,
            }
        )
    return rows


def _status_rows_for_direct_mz_file_failure(
    *,
    raw_file_id: str,
    targets_records: List[Dict[str, Any]],
    lc_mode: str,
    polarity: str,
    partial_filename: str,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for record in targets_records:
        num_value_out: Optional[Union[int, str]] = None
        if "num" in record:
            num_text = _normalize_num_value(record.get("num"))
            if num_text and num_text.isdigit():
                num_value_out = int(num_text)
            elif num_text:
                num_value_out = num_text

        mz_raw = record.get("m/z")
        mz_target: Optional[float]
        if mz_raw is None:
            mz_target = None
        else:
            mz_val = pd.to_numeric(mz_raw, errors="coerce")
            mz_target = float(mz_val) if pd.notna(mz_val) else None

        compound_raw = record.get("Compound name")
        compound_name = "" if pd.isna(compound_raw) else str(compound_raw).strip()
        if not compound_name:
            compound_name = "Unknown"

        rows.append(
            {
                "RawFile": raw_file_id,
                "num": num_value_out,
                "File name": str(partial_filename),
                "lc_mode": str(lc_mode),
                "mixture": _normalize_mixture_value(record.get("mixture")),
                "Compound name": compound_name,
                "Polarity": polarity,
                "mz_target": mz_target,
                "RT_min": None,
                "Intensity": None,
                "Area": None,
                "GaussianScore": None,
                "PeakQuality": None,
                "HasMS2": None,
                "EICGenerated": False,
                "FilteredOut": False,
            }
        )
    return rows


def _process_single_file_formula_worker(
    *,
    raw_file_path: str,
    formulas: List[str],
    mode: str,
    config_dict: Dict[str, Any],
    raw_file_id: str,
) -> Dict[str, Any]:
    config = Config.from_dict(config_dict)
    status_rows: List[Dict[str, Any]] = []
    with RawFileReader(raw_file_path) as reader:
        results = process_raw_file(
            reader,
            formulas,
            mode,
            config,
            raw_file_id=raw_file_id,
            status_rows=status_rows,
        )
    return {"results": results, "status_rows": status_rows, "raw_file_id": raw_file_id}


def _process_single_file_direct_mz_worker(
    *,
    raw_file_path: str,
    targets_records: List[Dict[str, Any]],
    config_dict: Dict[str, Any],
    lc_mode: str,
    polarity: str,
    partial_filename: str,
    file_suffix: str,
    raw_file_id: str,
) -> Dict[str, Any]:
    config = Config.from_dict(config_dict)
    targets_df = pd.DataFrame.from_records(targets_records)
    status_rows: List[Dict[str, Any]] = []
    with RawFileReader(str(raw_file_path)) as reader:
        results = process_raw_file_direct_mz(
            reader=reader,
            targets_df=targets_df,
            lc_mode=lc_mode,
            polarity=polarity,
            partial_filename=partial_filename,
            file_suffix=file_suffix,
            config=config,
            raw_file_id=raw_file_id,
            status_rows=status_rows,
        )
    return {"results": results, "status_rows": status_rows, "raw_file_id": raw_file_id}


def _process_common_batch(
    *,
    config: Config,
    work_items: List[Dict[str, Any]],
    worker_fn: Callable[..., Dict[str, Any]],
    failure_status_fn: Callable[[Dict[str, Any]], List[Dict[str, Any]]],
    progress_desc: str = "Processing raw files",
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Run a batch of work items, optionally using multiprocessing.

    This helper centralizes the pool/sequential orchestration and ensures output
    order matches the input ``work_items`` order (important for deterministic
    reports and tests).

    Args:
        config: Processing configuration (controls parallel settings and logging).
        work_items: List of item dictionaries passed to ``worker_fn`` via kwargs.
        worker_fn: Callable that accepts ``**item`` and returns a dict with
            optional "results" and "status_rows" keys.
        failure_status_fn: Callable that builds per-target status rows for an
            item when ``worker_fn`` raises.
        progress_desc: Label for the progress bar.

    Returns:
        Tuple of (results, status_rows).
    """
    all_results: List[Dict[str, Any]] = []
    all_status_rows: List[Dict[str, Any]] = []

    if not work_items:
        return all_results, all_status_rows

    use_pool = should_use_process_pool(
        parallel_mode=config.parallel_mode,
        num_workers=config.num_workers,
        num_items=len(work_items),
    )

    max_workers = None
    if use_pool:
        max_workers = resolve_max_workers(config.num_workers, num_items=len(work_items))
        logger.info(
            "Processing %d raw files with %d workers (spawn)",
            len(work_items),
            max_workers or 1,
        )
    else:
        logger.info(
            "Processing %d raw files sequentially (parallel_mode=%s, num_workers=%d)",
            len(work_items),
            config.parallel_mode,
            config.num_workers,
        )

    if use_pool:
        try:
            results_by_item: List[List[Dict[str, Any]]] = [[] for _ in work_items]
            status_by_item: List[List[Dict[str, Any]]] = [[] for _ in work_items]
            with create_process_pool(
                max_workers=max_workers,
                log_level=config.log_level,
                log_handlers=list(logging.getLogger().handlers),
            ) as executor:
                futures = {executor.submit(worker_fn, **item): idx for idx, item in enumerate(work_items)}
                succeeded = 0
                failed = 0
                show_progress = should_show_progress(bool(config.show_progress))
                with progress_bar(
                    total=len(work_items),
                    desc=progress_desc,
                    enabled=show_progress,
                    unit="file",
                ) as pbar:
                    for future in as_completed(futures):
                        idx = futures[future]
                        item = work_items[idx]
                        pbar.set_postfix_str(str(item.get("raw_file_id", "")))
                        try:
                            result = future.result()
                        except Exception as e:
                            failed += 1
                            logger.error(
                                "Failed to process file %s: %s",
                                item.get("raw_file_id", ""),
                                e,
                            )
                            status_by_item[idx] = failure_status_fn(item)
                        else:
                            succeeded += 1
                            results_by_item[idx] = result.get("results", [])
                            status_by_item[idx] = result.get("status_rows", [])
                        finally:
                            pbar.update(1)

                logger.info(
                    "Processed %d/%d raw files (%d failed)",
                    succeeded,
                    len(work_items),
                    failed,
                )

            for idx in range(len(work_items)):
                all_results.extend(results_by_item[idx])
                all_status_rows.extend(status_by_item[idx])

            return all_results, all_status_rows
        except Exception as e:
            logger.warning("Parallel processing failed, falling back to sequential: %s", e)
            use_pool = False

    if not use_pool:
        succeeded = 0
        failed = 0
        show_progress = should_show_progress(bool(config.show_progress))
        with progress_bar(
            total=len(work_items),
            desc=progress_desc,
            enabled=show_progress,
            unit="file",
        ) as pbar:
            for item in work_items:
                pbar.set_postfix_str(str(item.get("raw_file_id", "")))
                try:
                    result = worker_fn(**item)
                except Exception as e:
                    failed += 1
                    logger.error("Failed to process file %s: %s", item.get("raw_file_id", ""), e)
                    all_status_rows.extend(failure_status_fn(item))
                else:
                    succeeded += 1
                    all_results.extend(result.get("results", []))
                    all_status_rows.extend(result.get("status_rows", []))
                finally:
                    pbar.update(1)

        logger.info(
            "Processed %d/%d raw files (%d failed)",
            succeeded,
            len(work_items),
            failed,
        )

    return all_results, all_status_rows


def _write_results_if_any(*, config: Config, results: List[Dict[str, Any]], status_rows: List[Dict[str, Any]]) -> None:
    """Write results/status rows to the configured output Excel.

    Args:
        config: Processing configuration (includes output path and flags).
        results: List of extracted feature rows.
        status_rows: Per-target status rows, including filtered and failed
            extractions.
    """
    if results or status_rows:
        logger.info(
            "Saving %d results (%d targets) to: %s",
            len(results),
            len(status_rows),
            config.output_excel,
        )
        write_results_excel(
            results,
            config.output_excel,
            include_pivot_tables=bool(config.include_pivot_tables),
            status_rows=status_rows,
        )
        logger.info("Processing complete: %s", config.output_excel)
    else:
        logger.warning("No results to save")


class _LazyMS2Context:
    """Lazily build the MS2 index only when it is first needed."""

    def __init__(self, reader: RawFileReader, config: Config) -> None:
        self._reader = reader
        self._config = config
        self._index: Optional[Dict[str, Any]] = None
        self._built = False
        self._enabled = bool(config.enable_ms2 and reader.has_scan_events_api())

    @property
    def enabled(self) -> bool:
        return self._enabled

    def get_index(self) -> Optional[Dict[str, Any]]:
        if not self._enabled:
            return None

        if self._built:
            return self._index

        self._built = True
        try:
            self._index = build_ms2_index(self._reader)
        except Exception as exc:
            logger.warning("MS2 index build failed: %s", exc)
            self._enabled = False
            self._index = None
        else:
            logger.debug("Built MS2 index with %d entries", len(self._index.get("entries", [])))
        return self._index


def _build_ms2_context(reader: RawFileReader, config: Config) -> Tuple[bool, Optional[Dict[str, Any]]]:
    """Build an MS2 index when enabled and supported by the reader.

    Args:
        reader: Raw file reader.
        config: Processing configuration.

    Returns:
        Tuple of (ms2_enabled, ms2_index). If MS2 is disabled or unavailable,
        returns (False, None).
    """
    if not (config.enable_ms2 and reader.has_scan_events_api()):
        return False, None

    try:
        ms2_index = build_ms2_index(reader)
    except Exception as e:
        logger.warning("MS2 index build failed: %s", e)
        return False, None

    logger.debug("Built MS2 index with %d entries", len(ms2_index.get("entries", [])))
    return True, ms2_index


def _require_batch_chromatogram_api(reader: RawFileReader) -> None:
    """Validate that batch chromatogram extraction is available.

    Args:
        reader: Raw file reader.

    Raises:
        RuntimeError: If the batch chromatogram API is unavailable.
    """
    if not reader.has_multi_chromatogram_api():
        raise RuntimeError("Batch chromatogram extraction is required (multi-chromatogram API unavailable).")


def _get_chromatograms_batch(
    reader: RawFileReader,
    mz_list: List[float],
    ppm_tolerance: float,
    *,
    batch_size: int,
) -> List[Any]:
    """Extract chromatograms for multiple m/z values.

    Args:
        reader: Raw file reader.
        mz_list: Target m/z list.
        ppm_tolerance: m/z tolerance in ppm.

    Returns:
        List of (rt_arr, int_arr) tuples aligned to ``mz_list``.

    Raises:
        RuntimeError: If extraction fails.
    """
    try:
        try:
            return reader.get_chromatograms_batch(mz_list, ppm_tolerance, batch_size=batch_size)
        except TypeError as exc:
            # Backward-compatible fallback for test doubles / legacy readers.
            if "batch_size" in str(exc) and "unexpected keyword argument" in str(exc):
                return reader.get_chromatograms_batch(mz_list, ppm_tolerance)
            raise
    except Exception as e:
        raise RuntimeError(f"Batch chromatogram extraction failed: {e}") from e


def _summarize_eic_peak(
    rt_arr: np.ndarray,
    int_arr: np.ndarray,
    config: Config,
    *,
    max_intensity: float,
    filtered_out: bool,
    apex_idx: Optional[int] = None,
) -> Tuple[
    float,
    float,
    float,
    str,
    Optional[Tuple[float, float, float]],
    Optional[float],
]:
    """Summarize a single EIC into peak metrics.

    Args:
        rt_arr: Retention time array (minutes).
        int_arr: Intensity array.
        config: Processing configuration.
        max_intensity: Precomputed maximum intensity for ``int_arr``.
        filtered_out: Whether this peak is below the configured intensity threshold.
        apex_idx: Optional precomputed apex index into ``rt_arr``/``int_arr``.

    Returns:
        Tuple of (total_area, best_rt_min, gauss_score, quality_label, fit_params,
        rt_apex_for_ms2).
    """
    total_area = calculate_area(rt_arr, int_arr, method=config.area_method)

    best_rt_min = 0.0
    gauss_score = 0.0
    quality_label = "Noise"
    rt_apex_for_ms2 = None
    fit_params = None

    if not filtered_out and max_intensity > 1000:
        idx = int(apex_idx) if apex_idx is not None else int(np.argmax(int_arr)) if int_arr.size else 0
        if rt_arr.size:
            idx = max(0, min(idx, int(rt_arr.size - 1)))
            best_rt_min = float(rt_arr[idx])
        else:
            best_rt_min = 0.0
        rt_apex_for_ms2 = best_rt_min

        if config.enable_fitting:
            gauss_score, fit_params = fit_gaussian_and_score(
                rt_arr, int_arr, fit_rt_window_min=config.fit_rt_window_min
            )
            quality_label = score_to_quality_label(gauss_score, fitted=True)
        else:
            gauss_score = 0.0
            quality_label = score_to_quality_label(0.0, fitted=False)

    return (
        total_area,
        best_rt_min,
        gauss_score,
        quality_label,
        fit_params,
        rt_apex_for_ms2,
    )


def _maybe_match_ms2(
    *,
    ms2_context: _LazyMS2Context,
    mz_val: float,
    config: Config,
    filtered_out: bool,
    rt_apex_for_ms2: Optional[float],
) -> Tuple[Optional[bool], Optional[Dict[str, Any]]]:
    """Match MS2 events to an extracted EIC peak when enabled.

    Args:
        ms2_context: Lazy MS2 context (index is built on-demand).
        mz_val: Target m/z value.
        config: Processing configuration.
        filtered_out: Whether the EIC peak was filtered out.
        rt_apex_for_ms2: Optional apex retention time used for RT-window matching.

    Returns:
        Tuple of (has_ms2, ms2_match). If matching is skipped, returns (None, None).
    """
    if filtered_out or not ms2_context.enabled:
        return None, None

    ms2_index = ms2_context.get_index()
    if ms2_index is None:
        return None, None

    has_ms2, ms2_match = match_ms2(
        ms2_index,
        mz_val,
        config.ppm_tolerance,
        rt_apex_min=rt_apex_for_ms2,
        rt_window_min=config.ms2_rt_window_min,
        mode=config.ms2_match_mode,
    )
    return has_ms2, ms2_match


def _maybe_add_ms2_match_details(
    row_out: Dict[str, Any],
    *,
    filtered_out: bool,
    ms2_enabled: bool,
    ms2_match: Optional[Dict[str, Any]],
    config: Config,
) -> None:
    """Attach MS2 match fields to the output row when configured.

    Args:
        row_out: Output row to mutate in-place.
        filtered_out: Whether the peak was filtered out.
        ms2_enabled: Whether MS2 matching is enabled and supported.
        ms2_match: Optional MS2 match details.
        config: Processing configuration.
    """
    if filtered_out or not ms2_enabled or not config.store_ms2_match_details:
        return

    row_out["MS2ScanNo"] = ms2_match.get("scan_no") if ms2_match else None
    row_out["MS2RT_min"] = ms2_match.get("rt_min") if ms2_match else None
    row_out["MS2Precursor_mz"] = ms2_match.get("precursor_mz") if ms2_match else None


def process_raw_file(
    reader: RawFileReader,
    formulas: List[str],
    mode: str,
    config: Config,
    *,
    raw_file_id: Optional[str] = None,
    status_rows: Optional[List[Dict[str, Any]]] = None,
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
    filename = str(raw_file_id).strip() if raw_file_id else reader.filename
    logger.info("Processing: %s (%s)", filename, mode)

    results: List[Dict[str, Any]] = []

    ms2_context = _LazyMS2Context(reader, config)

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

    _require_batch_chromatogram_api(reader)

    target_mzs = [t.mz for t in targets]
    eic_results = _get_chromatograms_batch(
        reader,
        target_mzs,
        config.ppm_tolerance,
        batch_size=config.chromatogram_batch_size,
    )

    eic_dict = {t.key: result for t, result in zip(targets, eic_results)}
    logger.debug("Using batch chromatogram extraction")

    plotter: Any = None
    plotter_attempted = False

    for target in targets:
        formula = target.formula
        adduct_name = target.adduct
        target_mz = target.mz

        eic_data = eic_dict.get(target.key)
        if eic_data is None:
            logger.warning("Missing EIC for target: %s %s", formula, adduct_name)
            if status_rows is not None:
                status_rows.append(
                    {
                        "RawFile": filename,
                        "Mode": mode,
                        "Formula": formula,
                        "Adduct": adduct_name,
                        "mz_theoretical": target_mz,
                        "RT_min": None,
                        "Intensity": None,
                        "Area": None,
                        "GaussianScore": None,
                        "PeakQuality": None,
                        "HasMS2": None,
                        "EICGenerated": False,
                        "FilteredOut": False,
                    }
                )
            continue

        eic_rt, eic_int = eic_data

        max_intensity = float(np.max(eic_int)) if eic_int.size else 0.0

        # Skip peaks below threshold
        filtered_out = max_intensity < config.min_peak_intensity
        if filtered_out and status_rows is None:
            continue

        apex_idx: Optional[int] = None
        if not filtered_out and max_intensity > 1000 and eic_int.size:
            apex_idx = int(np.argmax(eic_int))

        (
            total_area,
            best_rt_min,
            gauss_score,
            quality_label,
            fit_params,
            rt_apex_for_ms2,
        ) = _summarize_eic_peak(
            eic_rt,
            eic_int,
            config,
            max_intensity=max_intensity,
            filtered_out=filtered_out,
            apex_idx=apex_idx,
        )

        if not filtered_out and max_intensity > 1000:
            # Save plot if enabled
            if config.enable_plotting:
                if not plotter_attempted:
                    plotter = create_eic_plotter()
                    plotter_attempted = True
                if plotter is not None:
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
                        apex_rt=best_rt_min,
                        max_intensity=max_intensity,
                        apex_idx=apex_idx,
                        plotter=plotter,
                    )

        # MS2 matching
        has_ms2, ms2_match = _maybe_match_ms2(
            ms2_context=ms2_context,
            mz_val=target_mz,
            config=config,
            filtered_out=filtered_out,
            rt_apex_for_ms2=rt_apex_for_ms2,
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
            "HasMS2": bool(has_ms2) if ms2_context.enabled else None,
            "EICGenerated": True,
            "FilteredOut": filtered_out,
        }

        _maybe_add_ms2_match_details(
            row_out,
            filtered_out=filtered_out,
            ms2_enabled=ms2_context.enabled,
            ms2_match=ms2_match,
            config=config,
        )

        if status_rows is not None:
            status_rows.append(row_out)
        if not filtered_out:
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
    *,
    raw_file_id: Optional[str] = None,
    status_rows: Optional[List[Dict[str, Any]]] = None,
) -> List[Dict[str, Any]]:
    """Process a single raw file using direct m/z targets from input Excel."""
    filename = str(raw_file_id).strip() if raw_file_id else reader.filename
    logger.info("Processing: %s (%s)", filename, polarity)

    results: List[Dict[str, Any]] = []

    ms2_context = _LazyMS2Context(reader, config)

    # Ensure m/z is numeric and filtered
    targets_df = targets_df.copy()
    targets_df["m/z"] = pd.to_numeric(targets_df["m/z"], errors="coerce")
    targets_df = targets_df.reset_index(drop=True)

    if targets_df.empty:
        return results

    num_width = 0
    if "num" in targets_df.columns:
        num_candidates = []
        for v in targets_df["num"].tolist():
            text = _normalize_num_value(v)
            if text and text.isdigit():
                num_candidates.append(text)
        if num_candidates:
            num_width = max(len(t) for t in num_candidates)

    valid_mz_mask = targets_df["m/z"].notna() & (targets_df["m/z"] != 0)
    valid_mz_indices = targets_df.index[valid_mz_mask].tolist()
    valid_index_to_pos = {idx: pos for pos, idx in enumerate(valid_mz_indices)}
    mz_list = targets_df.loc[valid_mz_mask, "m/z"].astype(float).tolist()

    eic_results_valid: List[Any] = []
    if mz_list:
        _require_batch_chromatogram_api(reader)
        eic_results_valid = _get_chromatograms_batch(
            reader,
            mz_list,
            config.ppm_tolerance,
            batch_size=config.chromatogram_batch_size,
        )

        logger.debug("Using batch chromatogram extraction")

    plotter = None
    plotter_attempted = False

    mz_values = targets_df["m/z"].to_numpy()
    num_values = targets_df["num"].to_numpy() if "num" in targets_df.columns else None
    compound_values = (
        targets_df["Compound name"].to_numpy()
        if "Compound name" in targets_df.columns
        else np.full(len(targets_df), np.nan, dtype=object)
    )
    mixture_values = (
        targets_df["mixture"].to_numpy()
        if "mixture" in targets_df.columns
        else np.full(len(targets_df), np.nan, dtype=object)
    )

    for i in range(len(targets_df)):
        mz_raw = mz_values[i]
        mz_val = float(mz_raw) if pd.notna(mz_raw) else None

        num_prefix = ""
        num_value_out: Optional[Union[int, str]] = None
        num_text = ""
        if num_values is not None:
            num_text = _normalize_num_value(num_values[i])
            if num_text and num_text.isdigit():
                num_value_out = int(num_text)
            elif num_text:
                num_value_out = num_text

        if num_width and num_values is not None:
            if num_text and num_text.isdigit():
                num_prefix = num_text.zfill(num_width)
            else:
                num_prefix = num_text

        compound_raw = compound_values[i]
        compound_name = "" if pd.isna(compound_raw) else str(compound_raw).strip()
        if not compound_name:
            compound_name = "Unknown"

        mixture = _normalize_mixture_value(mixture_values[i])

        row_out: Dict[str, Any] = {
            "RawFile": filename,
            "num": num_value_out,
            "File name": str(partial_filename),
            "lc_mode": str(lc_mode),
            "mixture": mixture,
            "Compound name": compound_name,
            "Polarity": polarity,
            "mz_target": mz_val,
            "RT_min": None,
            "Intensity": None,
            "Area": None,
            "GaussianScore": None,
            "PeakQuality": None,
            "HasMS2": None,
            "EICGenerated": False,
            "FilteredOut": False,
        }

        valid_pos = valid_index_to_pos.get(i)
        if mz_val is None or valid_pos is None:
            if status_rows is not None:
                status_rows.append(row_out)
            continue

        try:
            eic_rt, eic_int = eic_results_valid[valid_pos]
        except Exception as e:
            logger.warning(
                "EIC extraction error (%s, %s, %.4f): %s",
                filename,
                compound_name,
                float(mz_val) if mz_val is not None else float("nan"),
                e,
            )
            if status_rows is not None:
                status_rows.append(row_out)
            continue

        max_intensity = float(np.max(eic_int)) if eic_int.size else 0.0
        filtered_out = max_intensity < config.min_peak_intensity
        if filtered_out and status_rows is None:
            continue

        apex_idx: Optional[int] = None
        if not filtered_out and max_intensity > 1000 and eic_int.size:
            apex_idx = int(np.argmax(eic_int))

        (
            total_area,
            best_rt_min,
            gauss_score,
            quality_label,
            fit_params,
            rt_apex_for_ms2,
        ) = _summarize_eic_peak(
            eic_rt,
            eic_int,
            config,
            max_intensity=max_intensity,
            filtered_out=filtered_out,
            apex_idx=apex_idx,
        )

        if not filtered_out and max_intensity > 1000:
            if config.enable_plotting:
                from .io.plotting import save_eic_plot_direct_mz

                if not plotter_attempted:
                    plotter = create_eic_plotter()
                    plotter_attempted = True
                if plotter is not None:
                    save_eic_plot_direct_mz(
                        rt_arr=eic_rt,
                        int_arr=eic_int,
                        compound_name=compound_name,
                        polarity=polarity,
                        raw_filename=filename,
                        mz_val=mz_val,
                        mixture=mixture,
                        output_folder=config.export_plot_folder,
                        num_prefix=num_prefix,
                        lc_mode=lc_mode,
                        file_suffix=file_suffix,
                        partial_filename=partial_filename,
                        fit_params=fit_params,
                        score=gauss_score,
                        dpi=config.plot_dpi,
                        apex_rt=best_rt_min,
                        max_intensity=max_intensity,
                        apex_idx=apex_idx,
                        plotter=plotter,
                    )

        # MS2 matching
        has_ms2, ms2_match = _maybe_match_ms2(
            ms2_context=ms2_context,
            mz_val=float(mz_val) if mz_val is not None else float("nan"),
            config=config,
            filtered_out=filtered_out,
            rt_apex_for_ms2=rt_apex_for_ms2,
        )

        row_out.update(
            {
                "RT_min": round(best_rt_min, 3),
                "Intensity": max_intensity,
                "Area": total_area,
                "GaussianScore": round(gauss_score, 3),
                "PeakQuality": quality_label,
                "HasMS2": bool(has_ms2) if ms2_context.enabled else None,
                "EICGenerated": True,
                "FilteredOut": filtered_out,
            }
        )

        _maybe_add_ms2_match_details(
            row_out,
            filtered_out=filtered_out,
            ms2_enabled=ms2_context.enabled,
            ms2_match=ms2_match,
            config=config,
        )

        if status_rows is not None:
            status_rows.append(row_out)
        if not filtered_out:
            results.append(row_out)

    logger.info("Found %d features in %s", len(results), filename)
    return results


def _build_formula_work_items(
    *,
    grouped: Any,
    total_files: int,
    raw_entries: List[Path],
    config: Config,
    config_dict: Dict[str, Any],
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Build work items for formula-based processing.

    Args:
        grouped: Input rows grouped by (RawFile, Mode).
        total_files: Total number of file groups (for progress logs).
        raw_entries: Pre-indexed ``.raw`` entries under the raw data folder.
        config: Processing configuration.
        config_dict: Serialized config passed to worker processes.

    Returns:
        Tuple of (work_items, status_rows). The returned status rows contain
        pre-processing failures (e.g., missing raw files).
    """
    work_items: List[Dict[str, Any]] = []
    status_rows: List[Dict[str, Any]] = []
    raw_name_index = _build_raw_entry_name_index(raw_entries)

    for idx, ((raw_filename, mode), group_df) in enumerate(grouped):
        raw_filename = str(raw_filename).strip()
        if not raw_filename:
            logger.warning("Skipping empty RawFile at group %d", idx + 1)
            continue

        raw_filename = raw_filename.replace("\\", "/").strip()
        raw_filename_path = Path(raw_filename)
        if raw_filename_path.is_absolute() or any(part == ".." for part in raw_filename_path.parts):
            logger.error("Invalid RawFile path (must be relative): %s", raw_filename)
            continue

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
                idx + 1,
                total_files,
                raw_filename,
            )
            continue

        if not raw_filename_lower.endswith(".raw"):
            raw_filename = raw_filename + ".raw"

        # Build full path (supports either relative paths, or fallback search
        # within a nested raw folder structure).
        full_path = Path(config.raw_data_folder) / raw_filename
        full_file_path = str(full_path)
        raw_file_id = _build_raw_file_id(full_path, config.raw_data_folder)

        # Parse formulas
        formulas_raw = group_df["Formula"].dropna().astype(str).tolist()
        formulas_norm = [normalize_formula_str(f) for f in formulas_raw]
        formulas_norm = [f for f in formulas_norm if f]
        formulas = dedupe_preserve_order(formulas_norm)

        logger.info("[%d/%d] File: %s", idx + 1, total_files, raw_filename)

        if not os.path.exists(full_file_path):
            target_name = Path(raw_filename).name
            candidates = raw_name_index.get(target_name.lower(), [])
            if len(candidates) == 1:
                full_file_path = str(candidates[0])
                raw_file_id = _build_raw_file_id(candidates[0], config.raw_data_folder)
                logger.info("Resolved raw file under nested folders: %s", full_file_path)
            elif not candidates:
                logger.error("File not found: %s", full_file_path)
                status_rows.extend(
                    _status_rows_for_formula_file_failure(
                        raw_file_id=raw_file_id,
                        mode=mode,
                        formulas=formulas,
                    )
                )
                continue
            else:
                logger.error(
                    "Multiple raw files match %s under %s; specify a relative path in the input Excel. Matches: %s",
                    target_name,
                    config.raw_data_folder,
                    ", ".join(str(p) for p in candidates[:10]),
                )
                status_rows.extend(
                    _status_rows_for_formula_file_failure(
                        raw_file_id=raw_file_id,
                        mode=mode,
                        formulas=formulas,
                    )
                )
                continue

        work_items.append(
            {
                "raw_file_path": str(full_file_path),
                "formulas": formulas,
                "mode": mode,
                "config_dict": config_dict,
                "raw_file_id": raw_file_id,
            }
        )

    return work_items, status_rows


def _build_direct_mz_work_items(
    *,
    lc_mode_data: Dict[str, pd.DataFrame],
    total_groups: int,
    config: Config,
    config_dict: Dict[str, Any],
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Build work items for direct m/z processing.

    Args:
        lc_mode_data: Mapping of LC mode name to the corresponding target table.
        total_groups: Total number of (File name, Polarity) groups across all LC modes.
        config: Processing configuration.
        config_dict: Serialized config passed to worker processes.

    Returns:
        Tuple of (work_items, status_rows). The returned status rows contain
        pre-processing failures (e.g., missing raw file matches).
    """
    work_items: List[Dict[str, Any]] = []
    status_rows: List[Dict[str, Any]] = []
    group_counter = 0

    for lc_mode, meta_data in lc_mode_data.items():
        logger.info("Processing LC mode: %s (%d rows)", lc_mode, len(meta_data))

        raw_search_root = _resolve_lc_mode_raw_folder(config.raw_data_folder, lc_mode)
        raw_entries = _collect_raw_entries_recursive(raw_search_root)
        raw_entry_stem_index = _build_raw_entry_stem_index(raw_entries)
        logger.info("Indexed %d raw entries under: %s", len(raw_entries), raw_search_root)

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
                    "Invalid polarity for file %s (%s): %s",
                    partial_filename,
                    lc_mode,
                    e,
                )
                continue

            matching_files = find_matching_raw_files(
                partial_filename,
                raw_search_root,
                raw_entry_stem_index=raw_entry_stem_index,
            )
            if not matching_files:
                logger.warning("No matching raw files for: %s (%s)", partial_filename, lc_mode)
                for _, row in group_df.iterrows():
                    num_value_out: Optional[Union[int, str]] = None
                    if "num" in group_df.columns:
                        num_text = _normalize_num_value(row.get("num"))
                        if num_text and num_text.isdigit():
                            num_value_out = int(num_text)
                        elif num_text:
                            num_value_out = num_text

                    mz_raw = row.get("m/z")
                    mz_target: Optional[float]
                    if mz_raw is None:
                        mz_target = None
                    else:
                        mz_val = pd.to_numeric(mz_raw, errors="coerce")
                        mz_target = float(mz_val) if pd.notna(mz_val) else None

                    compound_raw = row.get("Compound name")
                    compound_name = "" if pd.isna(compound_raw) else str(compound_raw).strip()
                    if not compound_name:
                        compound_name = "Unknown"

                    status_rows.append(
                        {
                            "RawFile": None,
                            "num": num_value_out,
                            "File name": str(partial_filename),
                            "lc_mode": str(lc_mode),
                            "mixture": _normalize_mixture_value(row.get("mixture")),
                            "Compound name": compound_name,
                            "Polarity": polarity_norm,
                            "mz_target": mz_target,
                            "RT_min": None,
                            "Intensity": None,
                            "Area": None,
                            "GaussianScore": None,
                            "PeakQuality": None,
                            "HasMS2": None,
                            "EICGenerated": False,
                            "FilteredOut": False,
                        }
                    )
                continue

            targets_records = group_df.to_dict(orient="records")

            for raw_file_path in matching_files:
                run_label = _infer_run_label(raw_file_path, raw_search_root)
                file_suffix = extract_file_suffix(raw_file_path, partial_filename)
                if run_label:
                    run_suffix = f"_{run_label}"
                    if not file_suffix:
                        file_suffix = run_suffix
                    elif run_suffix.lower() not in file_suffix.lower():
                        file_suffix = f"{file_suffix}{run_suffix}"

                raw_file_id = _build_raw_file_id(raw_file_path, config.raw_data_folder)

                logger.info(
                    "[%d/%d] File: %s | %s (matched: %s)",
                    group_counter,
                    total_groups,
                    lc_mode,
                    partial_filename,
                    raw_file_id,
                )

                work_items.append(
                    {
                        "raw_file_path": str(raw_file_path),
                        "targets_records": targets_records,
                        "config_dict": config_dict,
                        "lc_mode": str(lc_mode),
                        "polarity": polarity_norm,
                        "partial_filename": partial_filename,
                        "file_suffix": file_suffix,
                        "raw_file_id": raw_file_id,
                    }
                )

    return work_items, status_rows


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
        raise FileNotFoundError(f"Input Excel file not found: {config.input_excel}")

    logger.info("Reading input file: %s", config.input_excel)
    meta_data = read_input_excel(
        config.input_excel,
        sheet_name=config.input_sheet,
        required_columns={"RawFile", "Mode", "Formula"},
    )

    all_results: List[Dict[str, Any]] = []
    all_status_rows: List[Dict[str, Any]] = []

    # Group by RawFile and Mode
    grouped = meta_data.groupby(["RawFile", "Mode"])
    total_files = len(grouped)

    logger.info("Processing %d file groups", total_files)

    raw_entries = _collect_raw_entries_recursive(config.raw_data_folder)

    config_dict = config.to_dict()
    work_items, pre_status_rows = _build_formula_work_items(
        grouped=grouped,
        total_files=total_files,
        raw_entries=raw_entries,
        config=config,
        config_dict=config_dict,
    )
    all_status_rows.extend(pre_status_rows)

    def _failure_status(item: Dict[str, Any]) -> List[Dict[str, Any]]:
        return _status_rows_for_formula_file_failure(
            raw_file_id=item["raw_file_id"],
            mode=item["mode"],
            formulas=item["formulas"],
        )

    batch_results, batch_status = _process_common_batch(
        config=config,
        work_items=work_items,
        worker_fn=_process_single_file_formula_worker,
        failure_status_fn=_failure_status,
    )
    all_results.extend(batch_results)
    all_status_rows.extend(batch_status)

    _write_results_if_any(config=config, results=all_results, status_rows=all_status_rows)


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
    all_status_rows: List[Dict[str, Any]] = []

    total_groups = sum(len(df.groupby(["File name", "Polarity"])) for df in lc_mode_data.values())

    logger.info("Processing %d file groups (direct m/z)", total_groups)
    config_dict = config.to_dict()
    work_items, pre_status_rows = _build_direct_mz_work_items(
        lc_mode_data=lc_mode_data,
        total_groups=total_groups,
        config=config,
        config_dict=config_dict,
    )
    all_status_rows.extend(pre_status_rows)

    def _failure_status(item: Dict[str, Any]) -> List[Dict[str, Any]]:
        return _status_rows_for_direct_mz_file_failure(
            raw_file_id=item["raw_file_id"],
            targets_records=item["targets_records"],
            lc_mode=item["lc_mode"],
            polarity=item["polarity"],
            partial_filename=item["partial_filename"],
        )

    batch_results, batch_status = _process_common_batch(
        config=config,
        work_items=work_items,
        worker_fn=_process_single_file_direct_mz_worker,
        failure_status_fn=_failure_status,
    )
    all_results.extend(batch_results)
    all_status_rows.extend(batch_status)

    _write_results_if_any(config=config, results=all_results, status_rows=all_status_rows)


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
