import os
import sys
import bisect
import math
import re
import warnings
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

# ==========================================
# [안전 장치] 라이브러리 설치 여부 확인
# ==========================================
try:
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from fisher_py import RawFile
    from fisher_py.data.business import TraceType, ChromatogramTraceSettings
    from fisher_py.data.business.mass_options import MassOptions
    from fisher_py.data.business.range import Range
    from molmass import Formula
    from scipy.optimize import curve_fit, OptimizeWarning
    from fisher_py.data.tolerance_units import ToleranceUnits
except ImportError as e:
    print("\n[!] 필수 라이브러리가 설치되지 않았거나, 경로가 잘못되었습니다.")
    print(f"    에러 상세: {e}")
    print("\n[해결 방법]")
    print("터미널(CMD 또는 PowerShell)에서 아래 명령어를 실행하여 설치해주세요:")
    print(
        ">> pip install pandas openpyxl molmass scipy matplotlib numpy fisher-py pythonnet"
    )
    print("\n(이미 설치했는데도 이 오류가 뜬다면, 파이썬 실행 경로를 확인해주세요.)")
    sys.exit(1)

# ==========================================
# 1. 환경 설정 및 상수 정의
# ==========================================

# [사용자 설정 필요] Thermo .raw 파일들이 들어있는 폴더 경로
# "./"는 현재 이 코드가 있는 폴더를 의미합니다.
# 데이터를 다른 곳에 두셨다면 아래 경로를 수정하세요. (예: r"C:\Data\MyProject")
RAW_DATA_FOLDER = r"./raw"

INPUT_EXCEL_FILE = "file_list.xlsx"  # 입력 엑셀 파일명
INPUT_SHEET_NAME = "Final"  # 읽어올 시트 이름
EXPORT_PLOT_FOLDER = "EIC_Plots_Export"  # 그래프 저장 폴더명
OUTPUT_EXCEL_FILE = "Final_Result_With_Plots.xlsx"  # 결과 엑셀 파일명

# 질량 상수
MASS_H = 1.0073
MASS_Na = 22.9892
MASS_ACN = 41.0265
MASS_FA = 46.0055
MASS_HCOO = 44.9983
MASS_NH4 = 18.0338
MASS_H2O = 18.0106
MASS_E = 0.00054858  # 전자 질량

# Adduct 정의
ADDUCT_DEFINITIONS = {
    # [n=2] Dimer Adducts
    # "[2M-2H+Na]-": {"multiplier": 2, "delta": -2 * MASS_H + MASS_Na, "net_charge": -1},
    # "[2M-H]-": {"multiplier": 2, "delta": -MASS_H, "net_charge": -1},
    # "[2M+H]+": {"multiplier": 2, "delta": +MASS_H, "net_charge": +1},
    # [n=1] Monomer Adducts
    # "[M-2H2O+H]+": {"multiplier": 1, "delta": -(2 * MASS_H2O) + MASS_H, "net_charge": +1},
    # "[M-3H2O+H]+": {"multiplier": 1, "delta": -(3 * MASS_H2O) + MASS_H, "net_charge": +1},
    "[M-H]-": {"multiplier": 1, "delta": -MASS_H, "net_charge": -1},
    # "[M-H2O-H]-": {"multiplier": 1, "delta": -(MASS_H2O + MASS_H), "net_charge": -1},
    # "[M-H2O+H]+": {"multiplier": 1, "delta": -MASS_H2O + MASS_H, "net_charge": +1},
    # "[M]-": {"multiplier": 1, "delta": 0, "net_charge": -1},
    # "[M]+": {"multiplier": 1, "delta": 0, "net_charge": +1},
    # "[M+ACN+H]+": {"multiplier": 1, "delta": +MASS_ACN + MASS_H, "net_charge": +1},
    # "[M+FA-H]-": {"multiplier": 1, "delta": +(MASS_FA - MASS_H), "net_charge": -1},
    "[M+H]+": {"multiplier": 1, "delta": +MASS_H, "net_charge": +1},
    # "[M+Na]+": {"multiplier": 1, "delta": +MASS_Na, "net_charge": +1},
    # "[M+NH4]+": {"multiplier": 1, "delta": +MASS_NH4, "net_charge": +1}
}

PPM_TOLERANCE = 10.0

# Gaussian fitting / plotting controls
ENABLE_FITTING = True
ENABLE_PLOTTING = True
MAX_PLOTS_PER_FILE = 25
FIT_RT_WINDOW_MIN = 0.30
PLOT_DPI = 120

# MS2 precursor matching
# - "rt_linked": only counts MS2 within ±MS2_RT_WINDOW_MIN around the EIC apex RT
# - "global": counts any MS2 match anywhere in the file (RT ignored)
MS2_MATCH_MODE = "rt_linked"  # "rt_linked" or "global"
MS2_RT_WINDOW_MIN = 0.30
STORE_MS2_MATCH_DETAILS = True  # include MS2ScanNo/MS2RT_min/MS2Precursor_mz in output

# Extension-ready (disabled by default)
EXPORT_MS2_MGF = False

# Peak area calculation method (must be consistent across runs).
# - "sum": legacy-compatible (matches previous sum(intensity) behavior)
# - "trapz": time-aware (trapz(intensity, rt_minutes))
AREA_METHOD = "sum"

# Minimum peak intensity threshold - peaks below this value are filtered out
MIN_PEAK_INTENSITY = 100000

# ---------
# Caches
# ---------
_FORMULA_EXACT_MASS_CACHE: Dict[str, Optional[float]] = {}


# ==========================================
# 2. 계산 및 분석 함수
# ==========================================

_FORMULA_WS_RE = re.compile(r"\s+")


def normalize_formula_str(value: str) -> str:
    value = "" if value is None else str(value)
    value = _FORMULA_WS_RE.sub("", value).strip()
    return value


def dedupe_preserve_order(items: Iterable[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for item in items:
        if item in seen:
            continue
        seen.add(item)
        out.append(item)
    return out


def get_exact_mass(formula_str: str) -> Optional[float]:
    formula_norm = normalize_formula_str(formula_str)
    if not formula_norm:
        return None

    cached = _FORMULA_EXACT_MASS_CACHE.get(formula_norm, None)
    if cached is not None or formula_norm in _FORMULA_EXACT_MASS_CACHE:
        return cached

    try:
        mass = float(Formula(formula_norm).isotope.mass)
    except Exception:
        mass = None

    _FORMULA_EXACT_MASS_CACHE[formula_norm] = mass
    return mass


def calculate_target_mz_from_mass(exact_mass: float, adduct_info) -> Optional[float]:
    try:
        exact_mass = float(exact_mass)
        if not math.isfinite(exact_mass):
            return None
        mult = adduct_info["multiplier"]
        delta = adduct_info["delta"]
        charge = adduct_info["net_charge"]
        ion_mass = (exact_mass * mult) + delta - (charge * MASS_E)
        return float(ion_mass / abs(charge))
    except Exception:
        return None


def calculate_target_mz(formula_str, adduct_info):
    exact_mass = get_exact_mass(formula_str)
    if exact_mass is None:
        return None
    return calculate_target_mz_from_mass(exact_mass, adduct_info)


def gaussian_func(x, a, x0, sigma):
    return a * np.exp(-((x - x0) ** 2) / (2 * sigma**2))


def fit_gaussian_and_score(
    rt_array, int_array, fit_rt_window_min: float = FIT_RT_WINDOW_MIN
):
    rt_array = np.asarray(rt_array, dtype=float)
    int_array = np.asarray(int_array, dtype=float)

    if rt_array.size < 5 or int_array.size < 5 or float(np.nanmax(int_array)) == 0.0:
        return 0.0, None
    try:
        max_idx = int(np.nanargmax(int_array))
        a_guess = float(int_array[max_idx])
        x0_guess = float(rt_array[max_idx])
        sigma_guess = 10.0 / 60.0

        finite_mask = np.isfinite(rt_array) & np.isfinite(int_array)
        if float(fit_rt_window_min) > 0:
            finite_mask &= (rt_array >= (x0_guess - float(fit_rt_window_min))) & (
                rt_array <= (x0_guess + float(fit_rt_window_min))
            )

        mask = finite_mask & (int_array > (a_guess * 0.1))
        if int(np.sum(mask)) < 5:
            mask = finite_mask

        x_data, y_data = rt_array[mask], int_array[mask]
        if x_data.size < 5 or float(np.nanmax(y_data)) == 0.0:
            return 0.0, None

        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            popt, _ = curve_fit(
                gaussian_func,
                x_data,
                y_data,
                p0=[a_guess, x0_guess, sigma_guess],
                maxfev=2000,
            )

        residuals = y_data - gaussian_func(x_data, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0.0
        return max(0.0, r_squared), popt
    except (ValueError, RuntimeError, OptimizeWarning):
        return 0.0, None


def save_eic_plot_png(
    rt_arr, int_arr, fit_params, score, formula, adduct, raw_filename, mz_val
):
    try:
        os.makedirs(EXPORT_PLOT_FOLDER, exist_ok=True)

        # Convert to relative abundance (0-100%)
        max_int = np.max(int_arr) if int_arr.size > 0 else 1.0
        if max_int == 0:
            max_int = 1.0
        rel_abundance = (int_arr / max_int) * 100.0

        # Calculate apex RT and peak height
        apex_idx = np.argmax(int_arr) if int_arr.size > 0 else 0
        apex_rt = rt_arr[apex_idx] if rt_arr.size > 0 else 0.0
        peak_height = max_int

        # Create figure with extra space at the bottom for annotations
        fig, ax = plt.subplots(figsize=(9, 4))

        ax.plot(rt_arr, rel_abundance, "b-", label="Raw EIC", linewidth=0.8, alpha=0.7)

        # Mark the apex point
        ax.plot(apex_rt, 100.0, "ro", markersize=6, label=f"Apex")

        # Add vertical line at apex RT
        ax.axvline(x=apex_rt, color="r", linestyle="--", linewidth=0.5, alpha=0.5)

        # Gaussian curve fitting disabled
        # if fit_params is not None:
        #     fitted_curve = gaussian_func(rt_arr, *fit_params)
        #     rel_fitted = (fitted_curve / max_int) * 100.0
        #     ax.plot(rt_arr, rel_fitted, 'r--', label=f'Fit (R2={score:.2f})', linewidth=0.8)

        # Title with m/z value
        ax.set_title(
            f"{formula} {adduct}  |  m/z = {mz_val:.4f}\nFile: {raw_filename}",
            fontsize=10,
        )
        ax.set_xlabel("Retention Time (min)")
        ax.set_ylabel("Relative Abundance (%)")

        ax.legend(loc="upper right", fontsize="small")
        ax.grid(True, linestyle=":", alpha=0.6)

        # Adjust subplot to make room at the bottom for annotations
        plt.subplots_adjust(bottom=0.22)

        # Add RT and peak height outside the chromatogram (below x-axis)
        annotation_text = (
            f"Apex RT: {apex_rt:.3f} min    |    Peak Height: {peak_height:.2e}"
        )
        fig.text(
            0.5,
            0.06,
            annotation_text,
            ha="center",
            va="top",
            fontsize=10,
            bbox=dict(
                boxstyle="round,pad=0.4",
                facecolor="lightyellow",
                edgecolor="gray",
                alpha=0.9,
            ),
        )

        safe_form = "".join(c for c in formula if c.isalnum())
        safe_add = (
            adduct.replace("[", "").replace("]", "").replace("+", "p").replace("-", "m")
        )
        safe_raw = sanitize_filename_component(raw_filename)
        save_name = f"{safe_form}_{safe_add}_{safe_raw}.png"

        save_path = os.path.join(EXPORT_PLOT_FOLDER, save_name)
        fig.savefig(save_path, dpi=PLOT_DPI)
    except Exception as e:
        print(f"Plot error: {e}")
    finally:
        plt.close("all")


_RT_UNIT_INFERRED = None


def sanitize_filename_component(value: str) -> str:
    value = os.path.basename(str(value))
    value = value.replace("/", "_").replace("\\", "_")
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("._")
    return value or "unknown"


def _infer_rt_unit(rt_arr: np.ndarray) -> str:
    rt_arr = np.asarray(rt_arr, dtype=float)
    if rt_arr.size == 0:
        return "minutes"

    rt_max = float(np.nanmax(rt_arr))
    if rt_max > 200:
        return "seconds"

    if rt_arr.size > 1:
        rt_sorted = np.sort(rt_arr)
        diffs = np.diff(rt_sorted)
        diffs = diffs[~np.isnan(diffs)]
        if diffs.size > 0 and float(np.nanmedian(diffs)) > 0.1:
            return "seconds"

    return "minutes"


def normalize_rt_to_minutes(rt_arr: np.ndarray) -> np.ndarray:
    global _RT_UNIT_INFERRED
    rt_arr = np.asarray(rt_arr, dtype=float)
    if rt_arr.size == 0:
        return rt_arr

    if _RT_UNIT_INFERRED is None:
        _RT_UNIT_INFERRED = _infer_rt_unit(rt_arr)
        rt_max = float(np.nanmax(rt_arr))
        if _RT_UNIT_INFERRED == "seconds":
            print(f"[RT] Inferred seconds (max={rt_max:.1f}); converting to minutes.")
        else:
            print(f"[RT] Inferred minutes (max={rt_max:.1f}); using minutes as-is.")

    if _RT_UNIT_INFERRED == "seconds":
        return rt_arr / 60.0
    return rt_arr


def calculate_area(rt_min: np.ndarray, intensity: np.ndarray) -> float:
    if AREA_METHOD == "trapz":
        if not intensity.size:
            return 0.0
        # Use np.trapezoid (NumPy 2.0+) if available, otherwise fall back to np.trapz
        trapz_func = getattr(np, "trapezoid", np.trapz)
        return float(trapz_func(intensity, rt_min))
    return float(np.sum(intensity)) if intensity.size else 0.0


def get_eic(raw_file_obj, target_mz, ppm_tol) -> Tuple[np.ndarray, np.ndarray]:
    rt_arr, int_arr = raw_file_obj.get_chromatogram(
        target_mz, ppm_tol, TraceType.MassRange
    )
    rt_min = normalize_rt_to_minutes(np.asarray(rt_arr, dtype=float))
    intensity = np.asarray(int_arr, dtype=float)

    n = min(rt_min.size, intensity.size)
    rt_min = rt_min[:n]
    intensity = intensity[:n]

    if rt_min.size > 1 and np.any(np.diff(rt_min) < 0):
        order = np.argsort(rt_min)
        rt_min = rt_min[order]
        intensity = intensity[order]

    return rt_min, intensity


# ==========================================
# 2.1 EIC extraction (multi-chromatogram)
# ==========================================


@dataclass(frozen=True)
class Target:
    formula: str
    adduct: str
    mz: float

    @property
    def key(self) -> Tuple[str, str]:
        return self.formula, self.adduct


def build_targets_for_file(formulas: Iterable[str], mode: str) -> List[Target]:
    targets: List[Target] = []
    for formula in formulas:
        exact_mass = get_exact_mass(formula)
        if exact_mass is None:
            continue

        for adduct_name, adduct_info in ADDUCT_DEFINITIONS.items():
            if (mode == "POS" and adduct_info["net_charge"] < 0) or (
                mode == "NEG" and adduct_info["net_charge"] > 0
            ):
                continue

            target_mz = calculate_target_mz_from_mass(exact_mass, adduct_info)
            if target_mz is None:
                continue
            targets.append(
                Target(formula=formula, adduct=adduct_name, mz=float(target_mz))
            )

    return targets


def _make_mass_range_settings(mz: float) -> ChromatogramTraceSettings:
    settings = ChromatogramTraceSettings()
    settings.trace = TraceType.MassRange
    settings.filter = "ms"
    settings.mass_ranges = [Range.create(float(mz), float(mz))]
    return settings


def _has_multi_chrom_api(raw_file_obj) -> bool:
    raw_access = getattr(raw_file_obj, "_raw_file_access", None)
    get_chrom_data = (
        getattr(raw_access, "get_chromatogram_data", None)
        if raw_access is not None
        else None
    )
    return callable(get_chrom_data)


def _chunked(items: List[Target], size: int) -> Iterable[List[Target]]:
    size = max(1, int(size))
    for i in range(0, len(items), size):
        yield items[i : i + size]


def extract_eics_multi(
    raw_file_obj,
    targets: List[Target],
    ppm_tol: float,
    batch_size: int = 256,
) -> Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]]:
    raw_access = getattr(raw_file_obj, "_raw_file_access", None)
    get_chrom_data = (
        getattr(raw_access, "get_chromatogram_data", None)
        if raw_access is not None
        else None
    )
    if not callable(get_chrom_data):
        raise RuntimeError(
            "get_chromatogram_data is unavailable on this RawFile object."
        )

    mass_opts = MassOptions()
    mass_opts.tolerance_units = ToleranceUnits.ppm
    mass_opts.tolerance = float(ppm_tol)

    out: Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]] = {}
    for batch in _chunked(targets, batch_size):
        settings_list = [_make_mass_range_settings(t.mz) for t in batch]
        chrom_data = get_chrom_data(settings_list, -1, -1, mass_opts)
        pos = chrom_data.positions_array
        ints = chrom_data.intensities_array

        if len(pos) != len(batch) or len(ints) != len(batch):
            raise RuntimeError(
                f"ChromatogramData size mismatch: expected {len(batch)} traces, got "
                f"positions={len(pos)}, intensities={len(ints)}."
            )

        for t, rt_arr, int_arr in zip(batch, pos, ints):
            rt_min = normalize_rt_to_minutes(np.asarray(rt_arr, dtype=float))
            intensity = np.asarray(int_arr, dtype=float)

            n = min(rt_min.size, intensity.size)
            rt_min = rt_min[:n]
            intensity = intensity[:n]

            if rt_min.size > 1 and np.any(np.diff(rt_min) < 0):
                order = np.argsort(rt_min)
                rt_min = rt_min[order]
                intensity = intensity[order]

            out[t.key] = (rt_min, intensity)

    return out


# ==========================================
# 2.2 MS2 (structured metadata only)
# ==========================================


def ppm_to_da(mz: float, ppm: float) -> float:
    return float(mz) * (float(ppm) / 1e6)


def _get_scan_range(raw_file_obj):
    first = getattr(raw_file_obj, "first_scan", None)
    last = getattr(raw_file_obj, "last_scan", None)
    if first is None or last is None:
        first_fn = getattr(raw_file_obj, "first_scan_number", None)
        last_fn = getattr(raw_file_obj, "last_scan_number", None)
        if callable(first_fn) and callable(last_fn):
            first, last = first_fn(), last_fn()
    if first is None or last is None:
        raise RuntimeError("Unable to determine scan range (first/last scan).")
    return int(first), int(last)


def _safe_get_retention_time_min(raw_file_obj, scan_no: int):
    fn = getattr(raw_file_obj, "get_retention_time_from_scan_number", None)
    if callable(fn):
        try:
            rt = fn(int(scan_no))
            rt = float(rt)
            if not math.isfinite(rt):
                return None
            return rt if rt <= 200 else rt / 60.0
        except Exception:
            return None

    raw_access = getattr(raw_file_obj, "_raw_file_access", None)
    rt_fn = (
        getattr(raw_access, "retention_time_from_scan_number", None)
        if raw_access is not None
        else None
    )
    if callable(rt_fn):
        try:
            rt = float(rt_fn(int(scan_no)))
            if not math.isfinite(rt):
                return None
            return rt if rt <= 200 else rt / 60.0
        except Exception:
            return None
    return None


def _has_scan_events_api(raw_file_obj) -> bool:
    raw_access = getattr(raw_file_obj, "_raw_file_access", None)
    get_scan_events = (
        getattr(raw_access, "get_scan_events", None) if raw_access is not None else None
    )
    return callable(get_scan_events)


def build_ms2_index_structured(raw_file_obj):
    """
    Build MS2 precursor index using structured scan-event metadata only.

    Returns:
        dict: {"entries": list[dict], "mz_list": list[float]}
    """
    from fisher_py.data.filter_enums import MsOrderType

    first_scan, last_scan = _get_scan_range(raw_file_obj)

    raw_access = getattr(raw_file_obj, "_raw_file_access", None)
    get_scan_events = (
        getattr(raw_access, "get_scan_events", None) if raw_access is not None else None
    )
    if not callable(get_scan_events):
        raise RuntimeError(
            "Structured scan-event access is unavailable on this RawFile object."
        )

    entries = []
    scan_events = get_scan_events(first_scan, last_scan)
    for i, scan_event in enumerate(scan_events):
        if getattr(scan_event, "ms_order", None) != MsOrderType.Ms2:
            continue

        scan_no = first_scan + i
        rt_min = _safe_get_retention_time_min(raw_file_obj, scan_no)

        try:
            mass_count = int(getattr(scan_event, "mass_count", 0))
        except Exception:
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

    entries.sort(key=lambda e: e["precursor_mz"])
    mz_list = [e["precursor_mz"] for e in entries]
    return {"entries": entries, "mz_list": mz_list}


def find_candidates_by_mz(ms2_index, target_mz: float, tol_da: float):
    if ms2_index is None:
        return []

    if isinstance(ms2_index, dict):
        entries = ms2_index.get("entries", [])
        mz_list = ms2_index.get("mz_list", [])
    else:
        entries = list(ms2_index)
        mz_list = [e.get("precursor_mz") for e in entries]

    if not entries or not mz_list:
        return []

    lo = bisect.bisect_left(mz_list, target_mz - tol_da)
    hi = bisect.bisect_right(mz_list, target_mz + tol_da)
    return entries[lo:hi]


def match_ms2(
    ms2_index,
    target_mz: float,
    ppm: float,
    rt_apex_min=None,
    rt_window_min: float = 0.30,
    mode: str = "rt_linked",
):
    if target_mz is None or not math.isfinite(float(target_mz)):
        return False, None

    tol_da = ppm_to_da(float(target_mz), float(ppm))
    candidates = find_candidates_by_mz(ms2_index, float(target_mz), tol_da)
    if not candidates:
        return False, None

    mode = (mode or "rt_linked").lower().strip()
    if mode == "rt_linked":
        if rt_apex_min is None or not math.isfinite(float(rt_apex_min)):
            return False, None

        rt_apex_min = float(rt_apex_min)
        in_window = []
        for e in candidates:
            rt = e.get("rt_min", None)
            if rt is None:
                continue
            try:
                rt = float(rt)
            except Exception:
                continue
            if not math.isfinite(rt):
                continue
            if abs(rt - rt_apex_min) <= float(rt_window_min):
                in_window.append(e)

        if not in_window:
            return False, None

        best = min(
            in_window,
            key=lambda e: (
                abs(float(e["rt_min"]) - rt_apex_min),
                abs(float(e["precursor_mz"]) - float(target_mz)),
            ),
        )
        return True, best

    if mode == "global":
        best = min(
            candidates, key=lambda e: abs(float(e["precursor_mz"]) - float(target_mz))
        )
        return True, best

    raise ValueError(
        f"Unknown MS2_MATCH_MODE: {mode!r} (expected 'rt_linked' or 'global')"
    )


def extract_ms2_spectrum(raw_file_obj, scan_no: int):
    """
    Extension hook: return (mz_array, intensity_array) for an MS2 scan number.
    """
    mz, intensity, _charge, _event = raw_file_obj.get_scan_from_scan_number(
        int(scan_no)
    )
    return mz, intensity


def export_ms2_to_mgf(*_args, **_kwargs):
    """
    Extension hook (optional): export matched MS2 to MGF.
    """
    raise NotImplementedError("MS2 MGF export is not implemented yet.")


def process_raw_file(
    raw_file_obj, raw_filename, formulas: List[str], output_data, mode
):
    filename = sanitize_filename_component(raw_filename)
    print(f"Processing: {filename} ({mode})")

    ms2_index = None
    ms2_enabled = False
    if _has_scan_events_api(raw_file_obj):
        try:
            ms2_index = build_ms2_index_structured(raw_file_obj)
            ms2_enabled = True
        except Exception as e:
            print(f"  - MS2 structured index build failed ({filename}): {e}")
            ms2_index = None

    targets = build_targets_for_file(formulas, mode)
    if not targets:
        return

    eic_dict: Optional[Dict[Tuple[str, str], Tuple[np.ndarray, np.ndarray]]] = None
    use_multi = False
    if _has_multi_chrom_api(raw_file_obj):
        try:
            eic_dict = extract_eics_multi(raw_file_obj, targets, PPM_TOLERANCE)
            use_multi = True
        except Exception:
            eic_dict = None
            use_multi = False

    plots_saved = 0

    for target in targets:
        formula = target.formula
        adduct_name = target.adduct
        target_mz = target.mz

        try:
            if use_multi:
                if eic_dict is None:
                    raise RuntimeError("Multi-chrom EIC dict is unavailable.")
                eic_rt, eic_int = eic_dict.get(target.key, (None, None))
                if eic_rt is None or eic_int is None:
                    raise RuntimeError("Missing EIC for target in multi-chrom results.")
            else:
                eic_rt, eic_int = get_eic(raw_file_obj, target_mz, PPM_TOLERANCE)
        except Exception as e:
            print(
                f"  - EIC extraction error ({filename}, {formula}, {adduct_name}): {e}"
            )
            continue

        max_intensity = float(np.max(eic_int)) if eic_int.size else 0.0

        # Skip peaks with max intensity below threshold
        if max_intensity < MIN_PEAK_INTENSITY:
            continue

        total_area = calculate_area(eic_rt, eic_int)

        best_rt_min = 0.0
        gauss_score = 0.0
        quality_label = "Noise"
        rt_apex_for_ms2 = None

        if max_intensity > 1000:
            best_rt_min = float(eic_rt[np.argmax(eic_int)])
            rt_apex_for_ms2 = best_rt_min
            fit_params = None
            if ENABLE_FITTING:
                gauss_score, fit_params = fit_gaussian_and_score(
                    eic_rt, eic_int, fit_rt_window_min=FIT_RT_WINDOW_MIN
                )

                if gauss_score > 0.8:
                    quality_label = "Excellent"
                elif gauss_score > 0.5:
                    quality_label = "Good"
                else:
                    quality_label = "Poor Shape"
            else:
                gauss_score = 0.0
                quality_label = "Not Fitted"

            if ENABLE_PLOTTING and plots_saved < int(MAX_PLOTS_PER_FILE):
                save_eic_plot_png(
                    eic_rt,
                    eic_int,
                    fit_params,
                    gauss_score,
                    formula,
                    adduct_name,
                    filename,
                    target_mz,
                )
                plots_saved += 1

        has_ms2 = None
        ms2_match = None
        if ms2_enabled:
            has_ms2, ms2_match = match_ms2(
                ms2_index,
                target_mz,
                PPM_TOLERANCE,
                rt_apex_min=rt_apex_for_ms2,
                rt_window_min=MS2_RT_WINDOW_MIN,
                mode=MS2_MATCH_MODE,
            )

        row_out = {
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

        if ms2_enabled and STORE_MS2_MATCH_DETAILS:
            row_out["MS2ScanNo"] = ms2_match.get("scan_no") if ms2_match else None
            row_out["MS2RT_min"] = ms2_match.get("rt_min") if ms2_match else None
            row_out["MS2Precursor_mz"] = (
                ms2_match.get("precursor_mz") if ms2_match else None
            )

        output_data.append(row_out)


# ==========================================
# 3. 메인 실행 (폴더 경로 적용)
# ==========================================


def main():

    if not os.path.exists(INPUT_EXCEL_FILE):
        print(
            f"오류: 엑셀 파일 '{INPUT_EXCEL_FILE}'이 스크립트와 같은 폴더에 없습니다."
        )
        return

    print("엑셀 파일을 읽는 중...")
    try:
        meta_data = pd.read_excel(INPUT_EXCEL_FILE, sheet_name=INPUT_SHEET_NAME)
    except Exception as e:
        print(f"오류: {e}")
        return

    required_cols = {"RawFile", "Mode", "Formula"}
    missing = sorted(required_cols.difference(set(meta_data.columns)))
    if missing:
        print("\n[오류] 입력 엑셀 형식이 올바르지 않습니다. (Long format만 지원)")
        print(f"  - 시트: {INPUT_SHEET_NAME!r}")
        print(f"  - 누락된 필수 컬럼: {missing}")
        print(f"  - 현재 컬럼: {list(meta_data.columns)}")
        print("  - 필요 컬럼: RawFile, Mode, Formula")
        return

    all_results = []

    # Group by RawFile and Mode to process each file once with all its formulas
    grouped = meta_data.groupby(["RawFile", "Mode"])
    total_files = len(grouped)

    for idx, ((raw_filename, mode), group_df) in enumerate(grouped):
        raw_filename = str(raw_filename).strip()
        raw_filename = sanitize_filename_component(raw_filename)
        raw_filename_lower = raw_filename.lower()
        if raw_filename_lower.endswith(".mzml"):
            print(f"\n[{idx + 1}/{total_files}] File: {raw_filename}")
            print(
                "  - (Error) .mzML 입력은 지원하지 않습니다. Thermo .raw 파일명을 사용해주세요."
            )
            continue
        if not raw_filename_lower.endswith(".raw"):
            raw_filename = raw_filename + ".raw"

        # [경로 결합] 지정된 폴더 + 파일명
        full_file_path = os.path.join(RAW_DATA_FOLDER, raw_filename)

        # 화학식 파싱 + 정규화/중복 제거(순서 유지)
        formulas_raw = group_df["Formula"].dropna().astype(str).tolist()
        formulas_norm = [normalize_formula_str(f) for f in formulas_raw]
        formulas_norm = [f for f in formulas_norm if f]
        formulas = dedupe_preserve_order(formulas_norm)

        print(f"\n[{idx + 1}/{total_files}] File: {raw_filename}")

        if os.path.exists(full_file_path):
            try:
                raw = RawFile(full_file_path)
            except Exception as e:
                print(f"  - (Error) RAW 파일을 열 수 없습니다: {e}")
                continue

            try:
                process_raw_file(raw, raw_filename, formulas, all_results, mode)
            finally:
                close_fn = getattr(raw, "close", None)
                if callable(close_fn):
                    close_fn()
        else:
            print(f"  - (Error) 파일을 찾을 수 없습니다.")
            print(f"  - 경로 확인: {full_file_path}")

    if all_results:
        print("\n데이터 저장 중...")
        df_results = pd.DataFrame(all_results)

        with pd.ExcelWriter(OUTPUT_EXCEL_FILE, engine="openpyxl") as writer:
            # 1. 전체 Raw Data 시트
            df_results.to_excel(writer, sheet_name="All_Features", index=False)

            unique_formulas = df_results["Formula"].unique()
            for formula in unique_formulas:
                f_data = df_results[df_results["Formula"] == formula]

                # 피벗 테이블 생성
                pivot_area = f_data.pivot_table(
                    index="RawFile", columns="Adduct", values="Area"
                )
                pivot_rt = f_data.pivot_table(
                    index="RawFile", columns="Adduct", values="RT_min"
                )
                pivot_intensity = f_data.pivot_table(
                    index="RawFile", columns="Adduct", values="Intensity"
                )

                # 시트명 생성
                safe_name = "".join(c for c in formula if c.isalnum())[:30]

                # -------------------------------------------------------
                # [1] Area Table 저장
                # -------------------------------------------------------
                pivot_area.to_excel(writer, sheet_name=safe_name, startrow=0)
                writer.sheets[safe_name].cell(row=1, column=1).value = "Area Table"

                # -------------------------------------------------------
                # [2] Retention Time Table 저장
                # -------------------------------------------------------

                # 다음 표가 시작될 위치 계산: (Area 표 길이) + (여백 3줄)
                current_row = len(pivot_area) + 3

                writer.sheets[safe_name].cell(row=current_row, column=1).value = (
                    "Retention Time (min)"
                )
                pivot_rt.to_excel(
                    writer, sheet_name=safe_name, startrow=current_row + 1
                )

                # -------------------------------------------------------
                # [3] Peak Height (Intensity) Table 저장
                # -------------------------------------------------------

                # 다음 표가 시작될 위치 계산: (Area 표 길이) + (RT 표 길이) + (여백 6줄)
                current_row = len(pivot_area) + len(pivot_rt) + 6

                writer.sheets[safe_name].cell(row=current_row, column=1).value = (
                    "Peak Height (Intensity)"
                )
                pivot_intensity.to_excel(
                    writer, sheet_name=safe_name, startrow=current_row + 1
                )

        print(f"\n[완료] 결과 엑셀: {OUTPUT_EXCEL_FILE}")
    else:
        print("\n[알림] 저장할 데이터가 없습니다.")


if __name__ == "__main__":
    main()
