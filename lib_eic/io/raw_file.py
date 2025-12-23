"""Raw file reader abstraction for Thermo .raw files."""

import logging
import os
import re
from typing import Any, List, Optional, Tuple, cast

import numpy as np

logger = logging.getLogger(__name__)

RT_UNIT_SECONDS_THRESHOLD = 200.0


def sanitize_filename_component(value: str) -> str:
    """Sanitize a string for use as a filename component.

    Args:
        value: String to sanitize.

    Returns:
        Sanitized string safe for use in filenames.
    """
    value = os.path.basename(str(value))
    value = value.replace("/", "_").replace("\\", "_")
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("._")
    return value or "unknown"


class RawFileReader:
    """Abstraction layer for reading Thermo .raw files.

    This class wraps fisher_py's low-level RawFileAccess to provide a cleaner interface
    and avoid expensive eager indexing done by fisher_py's higher-level RawFile wrapper.
    """

    _raw_access: Any

    def __init__(self, file_path: str):
        """Open a raw file.

        Args:
            file_path: Path to the Thermo .raw file.

        Raises:
            FileNotFoundError: If file doesn't exist.
            RuntimeError: If file cannot be opened.
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Raw file not found: {file_path}")

        try:
            from fisher_py.data import Device
            from fisher_py.raw_file_reader import RawFileReaderAdapter

            raw_access = RawFileReaderAdapter.file_factory(file_path)
            raw_access.select_instrument(Device.MS, 1)
            self._raw_access: Any = raw_access
        except ImportError as e:
            raise ImportError(
                "fisher-py is required for reading .raw files. " "Install with: pip install fisher-py pythonnet"
            ) from e
        except Exception as e:
            raise RuntimeError(f"Failed to open raw file: {e}") from e

        self._file_path = file_path
        self._rt_unit_inferred: Optional[str] = None

    @property
    def file_path(self) -> str:
        """Get the file path."""
        return self._file_path

    @property
    def filename(self) -> str:
        """Get the sanitized filename."""
        return sanitize_filename_component(os.path.basename(self._file_path))

    def close(self) -> None:
        """Close the raw file."""
        dispose_fn = getattr(self._raw_access, "dispose", None)
        if callable(dispose_fn):
            dispose_fn()

    def __enter__(self) -> "RawFileReader":
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """Context manager exit."""
        self.close()

    def _infer_rt_unit(self, rt_arr: np.ndarray) -> str:
        """Infer whether RT values are in seconds or minutes.

        Args:
            rt_arr: Array of retention time values.

        Returns:
            "seconds" or "minutes".
        """
        rt_arr = np.asarray(rt_arr, dtype=float)
        if rt_arr.size == 0:
            return "minutes"

        rt_max = float(np.nanmax(rt_arr))
        # Heuristic: values above typical LC run lengths (in minutes) likely mean
        # the vendor API is returning RT in seconds.
        if rt_max > RT_UNIT_SECONDS_THRESHOLD:
            return "seconds"

        if rt_arr.size > 1:
            rt_sorted = np.sort(rt_arr)
            diffs = np.diff(rt_sorted)
            diffs = diffs[~np.isnan(diffs)]
            if diffs.size > 0 and float(np.nanmedian(diffs)) > 0.1:
                return "seconds"

        return "minutes"

    def _normalize_rt_to_minutes(self, rt_arr: np.ndarray) -> np.ndarray:
        """Normalize RT array to minutes.

        Args:
            rt_arr: Array of retention time values.

        Returns:
            Array with RT values in minutes.
        """
        rt_arr = np.asarray(rt_arr, dtype=float)
        if rt_arr.size == 0:
            return rt_arr

        # Infer unit on first call
        if self._rt_unit_inferred is None:
            self._rt_unit_inferred = self._infer_rt_unit(rt_arr)
            rt_max = float(np.nanmax(rt_arr))
            if self._rt_unit_inferred == "seconds":
                logger.debug(
                    "RT unit inferred as seconds (max=%.1f); converting to minutes",
                    rt_max,
                )
            else:
                logger.debug(
                    "RT unit inferred as minutes (max=%.1f); using as-is",
                    rt_max,
                )

        if self._rt_unit_inferred == "seconds":
            return rt_arr / 60.0
        return rt_arr

    def get_chromatogram(
        self,
        target_mz: float,
        ppm_tolerance: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Extract an EIC (Extracted Ion Chromatogram) for a target m/z.

        Args:
            target_mz: Target m/z value.
            ppm_tolerance: Mass tolerance in ppm.

        Returns:
            Tuple of (rt_minutes, intensity) arrays.
        """
        from fisher_py.data.business import ChromatogramTraceSettings, TraceType
        from fisher_py.data.business.mass_options import MassOptions
        from fisher_py.data.business.range import Range
        from fisher_py.data.tolerance_units import ToleranceUnits

        mass_opts = MassOptions()
        mass_opts.tolerance_units = ToleranceUnits.ppm
        mass_opts.tolerance = float(ppm_tolerance)

        settings = ChromatogramTraceSettings()
        settings.trace = TraceType.MassRange
        settings.filter = "ms"
        settings.mass_ranges = [Range.create(float(target_mz), float(target_mz))]

        chrom_data = self._raw_access.get_chromatogram_data([settings], -1, -1, mass_opts)
        rt_arr = chrom_data.positions_array[0]
        int_arr = chrom_data.intensities_array[0]

        rt_min = self._normalize_rt_to_minutes(np.asarray(rt_arr, dtype=float))
        intensity = np.asarray(int_arr, dtype=float)

        # Handle size mismatch
        n = min(rt_min.size, intensity.size)
        rt_min = rt_min[:n]
        intensity = intensity[:n]

        # Sort by RT if needed
        if rt_min.size > 1 and np.any(np.diff(rt_min) < 0):
            order = np.argsort(rt_min)
            rt_min = rt_min[order]
            intensity = intensity[order]

        return rt_min, intensity

    def has_multi_chromatogram_api(self) -> bool:
        """Check if the multi-chromatogram API is available.

        Returns:
            True if batch chromatogram extraction is supported.
        """
        get_chrom_data = getattr(self._raw_access, "get_chromatogram_data", None)
        return callable(get_chrom_data)

    def get_chromatograms_batch(
        self,
        target_mzs: List[float],
        ppm_tolerance: float,
        batch_size: int = 256,
    ) -> List[Tuple[np.ndarray, np.ndarray]]:
        """Extract multiple EICs in batch mode.

        Args:
            target_mzs: List of target m/z values.
            ppm_tolerance: Mass tolerance in ppm.
            batch_size: Number of chromatograms per batch.

        Returns:
            List of (rt_minutes, intensity) tuples.

        Raises:
            RuntimeError: If batch API is not available.
        """
        from fisher_py.data.business import ChromatogramTraceSettings, TraceType
        from fisher_py.data.business.mass_options import MassOptions
        from fisher_py.data.business.range import Range
        from fisher_py.data.tolerance_units import ToleranceUnits

        get_chrom_data = cast(Any, getattr(self._raw_access, "get_chromatogram_data", None))
        if not callable(get_chrom_data):
            raise RuntimeError("get_chromatogram_data is unavailable on this raw file access object.")

        mass_opts = MassOptions()
        mass_opts.tolerance_units = ToleranceUnits.ppm
        mass_opts.tolerance = float(ppm_tolerance)

        results: List[Tuple[np.ndarray, np.ndarray]] = []

        # Process in batches
        for i in range(0, len(target_mzs), batch_size):
            batch_mzs = target_mzs[i : i + batch_size]

            # Create settings for each m/z
            settings_list = []
            for mz in batch_mzs:
                settings = ChromatogramTraceSettings()
                settings.trace = TraceType.MassRange
                settings.filter = "ms"
                settings.mass_ranges = [Range.create(float(mz), float(mz))]
                settings_list.append(settings)

            # Extract chromatograms
            chrom_data: Any = get_chrom_data(settings_list, -1, -1, mass_opts)
            positions = chrom_data.positions_array
            intensities = chrom_data.intensities_array

            if len(positions) != len(batch_mzs) or len(intensities) != len(batch_mzs):
                raise RuntimeError(f"ChromatogramData size mismatch: expected {len(batch_mzs)} traces")

            shared_rt_min: Optional[np.ndarray] = None
            shared_rt_min_unsorted: Optional[np.ndarray] = None
            shared_sort_order: Optional[np.ndarray] = None
            can_share_rt = False

            if len(positions) > 0:
                try:
                    first_rt_raw = np.asarray(positions[0], dtype=float)
                except Exception:
                    first_rt_raw = None

                if first_rt_raw is not None:
                    shared_rt_min_unsorted = self._normalize_rt_to_minutes(first_rt_raw)
                    shared_rt_min = shared_rt_min_unsorted

                    if shared_rt_min.size > 1 and np.any(np.diff(shared_rt_min) < 0):
                        shared_sort_order = np.argsort(shared_rt_min)
                        shared_rt_min = shared_rt_min[shared_sort_order]

                    shared_len = int(shared_rt_min.size)
                    if shared_len == 0:
                        can_share_rt = True
                    else:
                        # Quick sampling check to ensure all traces share the same RT axis
                        sample_indices: List[int] = [0, shared_len - 1]
                        if shared_len >= 4:
                            sample_indices.extend([shared_len // 2, shared_len // 4, (3 * shared_len) // 4])
                        sample_indices = sorted(set(i for i in sample_indices if 0 <= i < shared_len))

                        scale = 1.0
                        if self._rt_unit_inferred == "seconds":
                            scale = 1.0 / 60.0

                        can_share_rt = True
                        for trace_idx, (pos_arr, int_arr) in enumerate(zip(positions, intensities)):
                            if len(pos_arr) != shared_len or len(int_arr) != shared_len:
                                can_share_rt = False
                                break
                            if trace_idx == 0 or not sample_indices:
                                continue
                            try:
                                for sample_idx in sample_indices:
                                    rt_val = float(pos_arr[sample_idx]) * scale
                                    if not np.isclose(
                                        rt_val,
                                        float(shared_rt_min_unsorted[sample_idx]),
                                        rtol=0.0,
                                        atol=1e-6,
                                    ):
                                        can_share_rt = False
                                        break
                            except Exception:
                                can_share_rt = False
                            if not can_share_rt:
                                break

            if can_share_rt and shared_rt_min is not None:
                for int_arr in intensities:
                    intensity = np.asarray(int_arr, dtype=float)
                    if shared_sort_order is not None:
                        intensity = intensity[shared_sort_order]
                    results.append((shared_rt_min, intensity))
            else:
                # Fall back to per-trace processing (handles mismatched RT/intensity sizes safely)
                for rt_arr, int_arr in zip(positions, intensities):
                    rt_min = self._normalize_rt_to_minutes(np.asarray(rt_arr, dtype=float))
                    intensity = np.asarray(int_arr, dtype=float)

                    # Handle size mismatch
                    n = min(rt_min.size, intensity.size)
                    rt_min = rt_min[:n]
                    intensity = intensity[:n]

                    # Sort by RT if needed
                    if rt_min.size > 1 and np.any(np.diff(rt_min) < 0):
                        order = np.argsort(rt_min)
                        rt_min = rt_min[order]
                        intensity = intensity[order]

                    results.append((rt_min, intensity))

        return results

    def has_scan_events_api(self) -> bool:
        """Check if the scan events API is available.

        Returns:
            True if MS2 scan event access is supported.
        """
        get_scan_events = cast(Any, getattr(self._raw_access, "get_scan_events", None))
        return callable(get_scan_events)

    def get_scan_range(self) -> Tuple[int, int]:
        """Get the first and last scan numbers.

        Returns:
            Tuple of (first_scan, last_scan).

        Raises:
            RuntimeError: If scan range cannot be determined.
        """
        try:
            run_header = getattr(self._raw_access, "run_header", None)
            if run_header is not None:
                first = getattr(run_header, "first_spectrum", None)
                last = getattr(run_header, "last_spectrum", None)
                if first is not None and last is not None:
                    return int(first), int(last)
        except Exception:
            pass

        raise RuntimeError("Unable to determine scan range")

    def get_retention_time_from_scan(self, scan_no: int) -> Optional[float]:
        """Get retention time for a scan number.

        Args:
            scan_no: Scan number.

        Returns:
            Retention time in minutes, or None if unavailable.
        """
        import math

        rt_fn = cast(Any, getattr(self._raw_access, "retention_time_from_scan_number", None))
        if callable(rt_fn):
            try:
                rt_raw: Any = rt_fn(int(scan_no))
                rt = float(rt_raw)
                if not math.isfinite(rt):
                    return None
                return rt if rt <= 200 else rt / 60.0
            except Exception:
                pass

        return None

    def get_scan_events(self, first_scan: int, last_scan: int) -> List[Any]:
        """Get scan events for a range of scans.

        Args:
            first_scan: First scan number.
            last_scan: Last scan number.

        Returns:
            List of scan events.

        Raises:
            RuntimeError: If scan events API is unavailable.
        """
        get_scan_events = cast(Any, getattr(self._raw_access, "get_scan_events", None))
        if not callable(get_scan_events):
            raise RuntimeError("get_scan_events is unavailable")

        return cast(List[Any], get_scan_events(first_scan, last_scan))

    def get_scan(self, scan_no: int) -> Tuple[np.ndarray, np.ndarray, Any, Any]:
        """Get spectrum data for a scan number.

        Args:
            scan_no: Scan number.

        Returns:
            Tuple of (mz_array, intensity_array, charge, scan_event).
        """
        from fisher_py.data.filter_enums import MassAnalyzerType

        scan_no = int(scan_no)
        scan_event = self._raw_access.get_scan_event_for_scan_number(scan_no)

        if getattr(scan_event, "mass_analyzer", None) == MassAnalyzerType.MassAnalyzerFTMS:
            spectrum = self._raw_access.get_centroid_stream(scan_no, False)
            mz_array = np.array(spectrum.masses)
            intensity_array = np.array(spectrum.intensities)
            charges = np.array(spectrum.charges)
        else:
            stats = self._raw_access.get_scan_stats_for_scan_number(scan_no)
            spectrum = self._raw_access.get_segmented_scan_from_scan_number(scan_no, stats)
            mz_array = np.array(spectrum.positions)
            intensity_array = np.array(spectrum.intensities)
            charges = np.zeros(mz_array.shape)

        scan_event_str = self._raw_access.get_scan_event_string_for_scan_number(scan_no)
        return mz_array, intensity_array, charges, scan_event_str
