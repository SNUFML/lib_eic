"""EIC plot generation for LCMS Adduct Finder."""

import logging
import os
from typing import Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

_DEFAULT_PNG_COMPRESS_LEVEL = 1  # faster (larger files) than matplotlib default


def _create_figure(figsize: Tuple[float, float] = (9, 4)):
    try:
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        from matplotlib.figure import Figure
    except ImportError:
        logger.error("matplotlib is required for plotting")
        return None, None

    fig = Figure(figsize=figsize)
    FigureCanvas(fig)
    ax = fig.add_subplot(111)
    return fig, ax


def _save_png(fig, save_path: str, *, dpi: int) -> None:
    try:
        fig.savefig(
            save_path,
            dpi=dpi,
            pil_kwargs={"compress_level": _DEFAULT_PNG_COMPRESS_LEVEL},
        )
    except TypeError:
        fig.savefig(save_path, dpi=dpi)


class _EICPlotter:
    def __init__(self) -> None:
        fig, ax = _create_figure(figsize=(9, 4))
        if fig is None or ax is None:
            raise RuntimeError("matplotlib is required for plotting")

        self.fig = fig
        self.ax = ax

        (self._line,) = ax.plot(
            [],
            [],
            "b-",
            label="Raw EIC",
            linewidth=0.8,
            alpha=0.7,
        )
        (self._apex_marker,) = ax.plot(
            [],
            [],
            "ro",
            markersize=6,
            label="Apex",
        )
        self._vline = ax.axvline(
            x=0.0,
            color="r",
            linestyle="--",
            linewidth=0.5,
            alpha=0.5,
        )
        self._title = ax.set_title("", fontsize=10)

        ax.set_xlabel("Retention Time (min)")
        ax.set_ylabel("Relative Abundance (%)")
        ax.set_ylim(0.0, 105.0)

        ax.legend(loc="upper right", fontsize="small")
        ax.grid(True, linestyle=":", alpha=0.6)

        fig.subplots_adjust(bottom=0.22)
        self._annotation = fig.text(
            0.5,
            0.06,
            "",
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

    def update(
        self,
        rt_arr: np.ndarray,
        rel_abundance: np.ndarray,
        *,
        apex_rt: float,
        title: str,
        annotation_text: str,
    ) -> None:
        rt_arr = np.asarray(rt_arr, dtype=float)
        rel_abundance = np.asarray(rel_abundance, dtype=float)
        apex_rt = float(apex_rt)

        self._line.set_data(rt_arr, rel_abundance)
        self._apex_marker.set_data([apex_rt], [100.0])
        self._vline.set_xdata([apex_rt, apex_rt])
        self._title.set_text(title)
        self._annotation.set_text(annotation_text)

        if rt_arr.size:
            x_min = float(rt_arr[0])
            x_max = float(rt_arr[-1])
            if x_min > x_max:
                x_min, x_max = x_max, x_min
            self.ax.set_xlim(x_min, x_max)
        else:
            self.ax.set_xlim(0.0, 1.0)


_EIC_PLOTTER: Optional[_EICPlotter] = None


def _get_eic_plotter() -> Optional[_EICPlotter]:
    global _EIC_PLOTTER
    if _EIC_PLOTTER is None:
        try:
            _EIC_PLOTTER = _EICPlotter()
        except Exception:
            _EIC_PLOTTER = None
    return _EIC_PLOTTER


def save_eic_plot(
    rt_arr: np.ndarray,
    int_arr: np.ndarray,
    formula: str,
    adduct: str,
    raw_filename: str,
    mz_val: float,
    output_folder: str,
    fit_params: Optional[Tuple[float, float, float]] = None,
    score: float = 0.0,
    dpi: int = 120,
) -> Optional[str]:
    """Save an EIC plot as a PNG file.

    Args:
        rt_arr: Retention time array (minutes).
        int_arr: Intensity array.
        formula: Chemical formula.
        adduct: Adduct name.
        raw_filename: Name of the raw file.
        mz_val: Target m/z value.
        output_folder: Folder to save plots.
        fit_params: Optional Gaussian fit parameters (a, x0, sigma).
        score: Gaussian fit R-squared score.
        dpi: Plot resolution.

    Returns:
        Path to saved plot, or None if saving failed.
    """
    plotter = _get_eic_plotter()
    if plotter is None:
        return None

    try:
        os.makedirs(output_folder, exist_ok=True)

        # Convert to relative abundance (0-100%)
        max_int = np.max(int_arr) if int_arr.size > 0 else 1.0
        if max_int == 0:
            max_int = 1.0
        rel_abundance = (int_arr / max_int) * 100.0

        # Calculate apex RT and peak height
        apex_idx = np.argmax(int_arr) if int_arr.size > 0 else 0
        apex_rt = rt_arr[apex_idx] if rt_arr.size > 0 else 0.0
        peak_height = max_int

        title_text = f"{formula} {adduct}  |  m/z = {mz_val:.4f}\nFile: {raw_filename}"
        annotation_text = (
            f"Apex RT: {apex_rt:.3f} min    |    Peak Height: {peak_height:.2e}"
        )
        plotter.update(
            rt_arr=rt_arr,
            rel_abundance=rel_abundance,
            apex_rt=float(apex_rt),
            title=title_text,
            annotation_text=annotation_text,
        )

        # Generate safe filename
        safe_form = "".join(c for c in formula if c.isalnum())
        safe_add = (
            adduct.replace("[", "").replace("]", "").replace("+", "p").replace("-", "m")
        )
        safe_raw = _sanitize_component(raw_filename)
        save_name = f"{safe_form}_{safe_add}_{safe_raw}.png"

        save_path = os.path.join(output_folder, save_name)
        _save_png(plotter.fig, save_path, dpi=dpi)

        logger.debug("Saved plot: %s", save_path)
        return save_path

    except Exception as e:
        logger.error("Failed to save plot: %s", e)
        return None


def _sanitize_component(value: str) -> str:
    """Sanitize a string for use in filenames."""
    import re

    value = os.path.basename(str(value))
    value = value.replace("/", "_").replace("\\", "_")
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("._")
    return value or "unknown"


def build_direct_mz_plot_filename(
    compound_name: str,
    polarity: str,
    mixture: str,
    *,
    file_suffix: str = "",
    num_prefix: str = "",
) -> str:
    """Build a filesystem-safe direct-m/z plot filename (including .png)."""
    compound_name_disp = str(compound_name).strip() or "Unknown"
    safe_compound = _sanitize_component(compound_name_disp)

    polarity_safe = _sanitize_component(polarity)
    polarity_disp = str(polarity_safe).strip() or "UNK"

    safe_mixture = _sanitize_component(mixture)

    if file_suffix:
        import re

        suffix_text = str(file_suffix)
        suffix_text = suffix_text.replace("/", "_").replace("\\", "_")
        suffix_text = re.sub(r"[^A-Za-z0-9._-]+", "_", suffix_text)
        suffix_text = re.sub(r"_+", "_", suffix_text).strip(".")
        safe_suffix = suffix_text
    else:
        safe_suffix = ""

    prefix_text = str(num_prefix).strip()
    prefix_text = "" if not prefix_text or prefix_text.lower() == "nan" else prefix_text
    if prefix_text:
        safe_prefix = _sanitize_component(prefix_text)
        return f"{safe_prefix}_{safe_compound}_{polarity_disp}_{safe_mixture}{safe_suffix}.png"

    return f"{safe_compound}_{polarity_disp}_{safe_mixture}{safe_suffix}.png"


def save_eic_plot_direct_mz(
    rt_arr: np.ndarray,
    int_arr: np.ndarray,
    compound_name: str,
    polarity: str,
    raw_filename: str,
    mz_val: float,
    mixture: str,
    output_folder: str,
    lc_mode: str = "",
    file_suffix: str = "",
    partial_filename: str = "",
    fit_params: Optional[Tuple[float, float, float]] = None,
    score: float = 0.0,
    dpi: int = 120,
    num_prefix: str = "",
) -> Optional[str]:
    """Save an EIC plot as a PNG file (direct m/z input format).

    Directory structure:
        {output_folder}/{lc_mode}/{polarity}/{partial_filename}/

    Filename:
        {num_prefix}_{compound_name}_{polarity}_{mixture}{file_suffix}.png
    """
    plotter = _get_eic_plotter()
    if plotter is None:
        return None

    del fit_params, score  # reserved for future plot overlays

    try:
        lc_mode_safe = _sanitize_component(lc_mode) if str(lc_mode).strip() else ""
        polarity = _sanitize_component(polarity)
        partial_filename_safe = _sanitize_component(partial_filename)
        plot_dir = os.path.join(
            output_folder,
            lc_mode_safe,
            polarity,
            partial_filename_safe,
        )
        os.makedirs(plot_dir, exist_ok=True)

        # Convert to relative abundance (0-100%)
        max_int = np.max(int_arr) if int_arr.size > 0 else 1.0
        if max_int == 0:
            max_int = 1.0
        rel_abundance = (int_arr / max_int) * 100.0

        # Calculate apex RT and peak height
        apex_idx = np.argmax(int_arr) if int_arr.size > 0 else 0
        apex_rt = rt_arr[apex_idx] if rt_arr.size > 0 else 0.0
        peak_height = max_int

        compound_name_disp = str(compound_name).strip() or "Unknown"
        lc_mode_disp = str(lc_mode).strip() or "UNK"
        polarity_disp = str(polarity).strip() or "UNK"

        title_text = (
            f"{compound_name_disp} | {lc_mode_disp} | {polarity_disp} | m/z = {mz_val:.4f}\n"
            f"File: {raw_filename}"
        )
        annotation_text = f"Apex RT: {float(apex_rt):.3f} min    |    Peak Height: {float(peak_height):.2e}"
        plotter.update(
            rt_arr=rt_arr,
            rel_abundance=rel_abundance,
            apex_rt=float(apex_rt),
            title=title_text,
            annotation_text=annotation_text,
        )

        save_name = build_direct_mz_plot_filename(
            compound_name=compound_name_disp,
            polarity=polarity_disp,
            mixture=mixture,
            file_suffix=file_suffix,
            num_prefix=num_prefix,
        )
        save_path = os.path.join(plot_dir, save_name)
        _save_png(plotter.fig, save_path, dpi=dpi)

        logger.debug("Saved plot: %s", save_path)
        return save_path

    except Exception as e:
        logger.error("Failed to save plot: %s", e)
        return None
