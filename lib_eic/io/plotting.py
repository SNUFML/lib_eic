"""EIC plot generation for LCMS Adduct Finder."""

import logging
import os
from typing import Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


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
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.error("matplotlib is required for plotting")
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

        # Create figure
        fig, ax = plt.subplots(figsize=(9, 4))

        # Plot EIC
        ax.plot(
            rt_arr,
            rel_abundance,
            "b-",
            label="Raw EIC",
            linewidth=0.8,
            alpha=0.7,
        )

        # Mark apex point
        ax.plot(apex_rt, 100.0, "ro", markersize=6, label="Apex")

        # Add vertical line at apex RT
        ax.axvline(x=apex_rt, color="r", linestyle="--", linewidth=0.5, alpha=0.5)

        # Title with m/z value
        ax.set_title(
            f"{formula} {adduct}  |  m/z = {mz_val:.4f}\nFile: {raw_filename}",
            fontsize=10,
        )
        ax.set_xlabel("Retention Time (min)")
        ax.set_ylabel("Relative Abundance (%)")

        ax.legend(loc="upper right", fontsize="small")
        ax.grid(True, linestyle=":", alpha=0.6)

        # Adjust subplot for annotation space
        plt.subplots_adjust(bottom=0.22)

        # Add RT and peak height annotation
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

        # Generate safe filename
        safe_form = "".join(c for c in formula if c.isalnum())
        safe_add = (
            adduct.replace("[", "")
            .replace("]", "")
            .replace("+", "p")
            .replace("-", "m")
        )
        safe_raw = _sanitize_component(raw_filename)
        save_name = f"{safe_form}_{safe_add}_{safe_raw}.png"

        save_path = os.path.join(output_folder, save_name)
        fig.savefig(save_path, dpi=dpi)

        logger.debug("Saved plot: %s", save_path)
        return save_path

    except Exception as e:
        logger.error("Failed to save plot: %s", e)
        return None

    finally:
        plt.close("all")


def _sanitize_component(value: str) -> str:
    """Sanitize a string for use in filenames."""
    import re
    value = os.path.basename(str(value))
    value = value.replace("/", "_").replace("\\", "_")
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("._")
    return value or "unknown"
