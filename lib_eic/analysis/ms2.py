"""MS2 precursor matching for LC-MS analysis."""

import bisect
import logging
import math
from typing import Any, Dict, List, Optional, Tuple

from ..chemistry.mass import ppm_to_da
from ..io.raw_file import RawFileReader

logger = logging.getLogger(__name__)


# Type alias for MS2 index structure
MS2Index = Dict[str, Any]


def build_ms2_index(reader: RawFileReader) -> MS2Index:
    """Build MS2 precursor index from raw file.

    Args:
        reader: RawFileReader instance.

    Returns:
        Dictionary with:
            - "entries": List of {scan_no, rt_min, precursor_mz}
            - "mz_list": Sorted list of precursor m/z values

    Raises:
        RuntimeError: If MS2 index cannot be built.
    """
    try:
        from fisher_py.data.filter_enums import MsOrderType
    except ImportError:
        raise ImportError("fisher_py is required for MS2 indexing")

    if not reader.has_scan_events_api():
        raise RuntimeError("Scan events API is not available")

    first_scan, last_scan = reader.get_scan_range()

    entries: List[Dict[str, Any]] = []
    scan_events = reader.get_scan_events(first_scan, last_scan)

    for i, scan_event in enumerate(scan_events):
        if getattr(scan_event, "ms_order", None) != MsOrderType.Ms2:
            continue

        scan_no = first_scan + i
        rt_min = reader.get_retention_time_from_scan(scan_no)

        try:
            mass_count = int(getattr(scan_event, "mass_count", 0))
        except (TypeError, ValueError):
            mass_count = 0

        for j in range(max(0, mass_count)):
            try:
                reaction = scan_event.get_reaction(j)
                pmz = getattr(reaction, "precursor_mass", None)
                if pmz is None:
                    continue
                pmz = float(pmz)
                if not math.isfinite(pmz):
                    continue
            except Exception:
                continue

            entries.append(
                {
                    "scan_no": int(scan_no),
                    "rt_min": float(rt_min) if rt_min is not None else None,
                    "precursor_mz": float(pmz),
                }
            )

    # Sort by m/z for binary search
    entries.sort(key=lambda e: e["precursor_mz"])
    mz_list = [e["precursor_mz"] for e in entries]

    logger.debug("Built MS2 index with %d entries", len(entries))
    return {"entries": entries, "mz_list": mz_list}


def find_candidates_by_mz(
    ms2_index: Optional[MS2Index],
    target_mz: float,
    tolerance_da: float,
) -> List[Dict[str, Any]]:
    """Find MS2 scan candidates within mass tolerance.

    Args:
        ms2_index: MS2 index from build_ms2_index().
        target_mz: Target m/z value.
        tolerance_da: Mass tolerance in Daltons.

    Returns:
        List of matching MS2 entries.
    """
    if ms2_index is None:
        return []

    entries = ms2_index.get("entries", [])
    mz_list = ms2_index.get("mz_list", [])

    if not entries or not mz_list:
        return []

    # Binary search for candidates
    lo = bisect.bisect_left(mz_list, target_mz - tolerance_da)
    hi = bisect.bisect_right(mz_list, target_mz + tolerance_da)

    return entries[lo:hi]


def match_ms2(
    ms2_index: Optional[MS2Index],
    target_mz: float,
    ppm_tolerance: float,
    rt_apex_min: Optional[float] = None,
    rt_window_min: float = 0.30,
    mode: str = "rt_linked",
) -> Tuple[bool, Optional[Dict[str, Any]]]:
    """Match MS2 scans to a target m/z.

    Args:
        ms2_index: MS2 index from build_ms2_index().
        target_mz: Target m/z value.
        ppm_tolerance: Mass tolerance in ppm.
        rt_apex_min: EIC apex retention time (minutes).
        rt_window_min: RT window for rt_linked mode.
        mode: Matching mode ("rt_linked" or "global").

    Returns:
        Tuple of (has_match, best_match) where:
            - has_match: True if a match was found.
            - best_match: Best matching MS2 entry or None.

    Raises:
        ValueError: If mode is invalid.
    """
    if target_mz is None or not math.isfinite(float(target_mz)):
        return False, None

    tolerance_da = ppm_to_da(float(target_mz), float(ppm_tolerance))
    candidates = find_candidates_by_mz(ms2_index, float(target_mz), tolerance_da)

    if not candidates:
        return False, None

    mode = (mode or "rt_linked").lower().strip()

    if mode == "rt_linked":
        if rt_apex_min is None or not math.isfinite(float(rt_apex_min)):
            return False, None

        rt_apex_min = float(rt_apex_min)

        # Filter to scans within RT window
        in_window = []
        for entry in candidates:
            rt = entry.get("rt_min")
            if rt is None:
                continue
            try:
                rt = float(rt)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(rt):
                continue
            if abs(rt - rt_apex_min) <= float(rt_window_min):
                in_window.append(entry)

        if not in_window:
            return False, None

        # Find best match (closest by RT, then by m/z)
        best = min(
            in_window,
            key=lambda e: (
                abs(float(e["rt_min"]) - rt_apex_min),
                abs(float(e["precursor_mz"]) - float(target_mz)),
            ),
        )
        return True, best

    elif mode == "global":
        # Find closest by m/z only
        best = min(
            candidates,
            key=lambda e: abs(float(e["precursor_mz"]) - float(target_mz)),
        )
        return True, best

    else:
        raise ValueError(f"Unknown MS2 match mode: {mode!r}. Expected 'rt_linked' or 'global'.")


def extract_ms2_spectrum(
    reader: RawFileReader,
    scan_no: int,
) -> Tuple[Any, Any]:
    """Extract MS2 spectrum data from a scan.

    Args:
        reader: RawFileReader instance.
        scan_no: MS2 scan number.

    Returns:
        Tuple of (mz_array, intensity_array).
    """
    mz, intensity, _charge, _event = reader.get_scan(scan_no)
    return mz, intensity
