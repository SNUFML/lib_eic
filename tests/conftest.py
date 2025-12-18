"""Pytest configuration and fixtures for LCMS Adduct Finder tests."""

import pytest
import numpy as np


@pytest.fixture
def sample_eic_data():
    """Generate sample EIC data for testing."""
    rng = np.random.default_rng(0)
    # Create a Gaussian-like peak
    rt = np.linspace(0, 10, 100)
    center = 5.0
    sigma = 0.5
    intensity = 1e6 * np.exp(-((rt - center) ** 2) / (2 * sigma**2))
    # Add some noise
    intensity += rng.normal(0, 1000, intensity.shape)
    intensity = np.maximum(intensity, 0)  # No negative intensities
    return rt, intensity


@pytest.fixture
def noisy_eic_data():
    """Generate noisy EIC data without clear peak."""
    rng = np.random.default_rng(1)
    rt = np.linspace(0, 10, 100)
    intensity = rng.normal(1000, 500, rt.shape)
    intensity = np.maximum(intensity, 0)
    return rt, intensity


@pytest.fixture
def sample_config():
    """Create a sample Config instance for testing."""
    from lib_eic.config import Config

    return Config(
        ppm_tolerance=10.0,
        min_peak_intensity=1000,
        enable_fitting=True,
        enable_plotting=False,
    )
