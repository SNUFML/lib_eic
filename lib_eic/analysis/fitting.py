"""Gaussian peak fitting for LC-MS chromatograms."""

import logging
import warnings
from typing import Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Default fitting parameters
DEFAULT_FIT_RT_WINDOW_MIN = 0.30
DEFAULT_SIGMA_GUESS = 10.0 / 60.0  # ~10 seconds in minutes
MIN_DATA_POINTS = 5


def gaussian_func(x: np.ndarray, a: float, x0: float, sigma: float) -> np.ndarray:
    """Gaussian function for curve fitting.

    Args:
        x: X values (retention times).
        a: Amplitude (peak height).
        x0: Center position (apex RT).
        sigma: Standard deviation (peak width).

    Returns:
        Y values (intensities).
    """
    return a * np.exp(-((x - x0) ** 2) / (2 * sigma**2))


def fit_gaussian_and_score(
    rt_array: np.ndarray,
    int_array: np.ndarray,
    fit_rt_window_min: float = DEFAULT_FIT_RT_WINDOW_MIN,
) -> Tuple[float, Optional[Tuple[float, float, float]]]:
    """Fit a Gaussian curve to EIC data and calculate R-squared score.

    Args:
        rt_array: Retention time array (minutes).
        int_array: Intensity array.
        fit_rt_window_min: RT window around apex for fitting.

    Returns:
        Tuple of (r_squared, fit_params) where:
            - r_squared: Goodness of fit (0-1), 0 if fitting failed.
            - fit_params: (a, x0, sigma) or None if fitting failed.
    """
    try:
        from scipy.optimize import curve_fit, OptimizeWarning
    except ImportError:
        logger.error("scipy is required for Gaussian fitting")
        return 0.0, None

    rt_array = np.asarray(rt_array, dtype=float)
    int_array = np.asarray(int_array, dtype=float)

    # Check minimum data requirements
    if rt_array.size < MIN_DATA_POINTS or int_array.size < MIN_DATA_POINTS:
        return 0.0, None

    max_intensity = float(np.nanmax(int_array))
    if max_intensity == 0.0:
        return 0.0, None

    try:
        # Find apex
        max_idx = int(np.nanargmax(int_array))
        a_guess = float(int_array[max_idx])
        x0_guess = float(rt_array[max_idx])
        sigma_guess = DEFAULT_SIGMA_GUESS

        # Create mask for valid data points
        finite_mask = np.isfinite(rt_array) & np.isfinite(int_array)

        # Apply RT window filter if specified
        if float(fit_rt_window_min) > 0:
            finite_mask &= (rt_array >= (x0_guess - float(fit_rt_window_min))) & (
                rt_array <= (x0_guess + float(fit_rt_window_min))
            )

        # Filter to points above 10% of max intensity for better fitting
        intensity_mask = finite_mask & (int_array > (a_guess * 0.1))
        if int(np.sum(intensity_mask)) < MIN_DATA_POINTS:
            intensity_mask = finite_mask

        x_data = rt_array[intensity_mask]
        y_data = int_array[intensity_mask]

        if x_data.size < MIN_DATA_POINTS or float(np.nanmax(y_data)) == 0.0:
            return 0.0, None

        # Perform curve fitting
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            popt, _ = curve_fit(
                gaussian_func,
                x_data,
                y_data,
                p0=[a_guess, x0_guess, sigma_guess],
                maxfev=800,
            )

        # Calculate R-squared
        fitted_values = gaussian_func(x_data, *popt)
        residuals = y_data - fitted_values
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0.0

        fit_params = (float(popt[0]), float(popt[1]), float(popt[2]))
        return max(0.0, r_squared), fit_params

    except (ValueError, RuntimeError, OptimizeWarning) as e:
        logger.debug("Gaussian fitting failed: %s", e)
        return 0.0, None


def score_to_quality_label(score: float, fitted: bool = True) -> str:
    """Convert R-squared score to quality label.

    Args:
        score: R-squared score (0-1).
        fitted: Whether fitting was performed.

    Returns:
        Quality label string.
    """
    if not fitted:
        return "Not Fitted"

    if score > 0.8:
        return "Excellent"
    elif score > 0.5:
        return "Good"
    else:
        return "Poor Shape"
