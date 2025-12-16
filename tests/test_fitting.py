"""Tests for analysis/fitting.py module."""

import pytest
import numpy as np
from lib_eic.analysis.fitting import (
    gaussian_func,
    fit_gaussian_and_score,
    score_to_quality_label,
)


class TestGaussianFunc:
    """Tests for gaussian_func function."""

    def test_peak_at_center(self):
        """Test Gaussian function peaks at center."""
        x = np.array([0, 1, 2, 3, 4])
        result = gaussian_func(x, a=100, x0=2, sigma=1)
        assert result[2] == pytest.approx(100, rel=1e-6)

    def test_symmetric(self):
        """Test Gaussian function is symmetric."""
        x = np.array([0, 1, 2, 3, 4])
        result = gaussian_func(x, a=100, x0=2, sigma=1)
        assert result[1] == pytest.approx(result[3], rel=1e-6)
        assert result[0] == pytest.approx(result[4], rel=1e-6)

    def test_sigma_affects_width(self):
        """Test sigma parameter affects peak width."""
        x = np.linspace(0, 10, 101)  # 101 points so center is exactly at index 50
        narrow = gaussian_func(x, a=100, x0=5, sigma=0.5)
        wide = gaussian_func(x, a=100, x0=5, sigma=2)
        # Both should have same height at exact center (x=5)
        center_idx = 50
        assert narrow[center_idx] == pytest.approx(100, rel=1e-6)
        assert wide[center_idx] == pytest.approx(100, rel=1e-6)
        # But narrow peak should fall off faster at offset positions
        offset_idx = 40  # x=4, which is 1 unit away from center
        assert narrow[offset_idx] < wide[offset_idx]


class TestFitGaussianAndScore:
    """Tests for fit_gaussian_and_score function."""

    def test_perfect_gaussian_fit(self):
        """Test fitting a perfect Gaussian peak."""
        # Generate perfect Gaussian data
        rt = np.linspace(0, 10, 100)
        intensity = gaussian_func(rt, a=1e6, x0=5.0, sigma=0.3)

        score, params = fit_gaussian_and_score(rt, intensity)

        assert score > 0.99, "Perfect Gaussian should have R² > 0.99"
        assert params is not None
        a, x0, sigma = params
        assert abs(x0 - 5.0) < 0.1, "Center should be near 5.0"

    def test_noisy_gaussian_fit(self, sample_eic_data):
        """Test fitting a noisy Gaussian peak."""
        rt, intensity = sample_eic_data

        score, params = fit_gaussian_and_score(rt, intensity)

        assert score > 0.5, "Noisy Gaussian should still have reasonable fit"
        assert params is not None

    def test_pure_noise_low_score(self, noisy_eic_data):
        """Test that pure noise has low score."""
        rt, intensity = noisy_eic_data

        score, params = fit_gaussian_and_score(rt, intensity)

        # Pure noise should have low R² (may vary)
        assert score < 0.8 or params is None

    def test_insufficient_data(self):
        """Test handling of insufficient data."""
        rt = np.array([1, 2, 3])
        intensity = np.array([100, 200, 100])

        score, params = fit_gaussian_and_score(rt, intensity)

        assert score == 0.0
        assert params is None

    def test_empty_data(self):
        """Test handling of empty data."""
        rt = np.array([])
        intensity = np.array([])

        score, params = fit_gaussian_and_score(rt, intensity)

        assert score == 0.0
        assert params is None

    def test_zero_intensity(self):
        """Test handling of all-zero intensity."""
        rt = np.linspace(0, 10, 100)
        intensity = np.zeros_like(rt)

        score, params = fit_gaussian_and_score(rt, intensity)

        assert score == 0.0
        assert params is None

    def test_rt_window_filtering(self):
        """Test RT window filtering around apex."""
        # Create two peaks
        rt = np.linspace(0, 20, 200)
        peak1 = gaussian_func(rt, a=1e6, x0=5.0, sigma=0.3)
        peak2 = gaussian_func(rt, a=0.5e6, x0=15.0, sigma=0.3)
        intensity = peak1 + peak2

        # Fit with window should focus on main peak
        score, params = fit_gaussian_and_score(rt, intensity, fit_rt_window_min=2.0)

        assert params is not None
        a, x0, sigma = params
        # Should fit the larger peak at x=5
        assert abs(x0 - 5.0) < 1.0


class TestScoreToQualityLabel:
    """Tests for score_to_quality_label function."""

    def test_excellent_label(self):
        assert score_to_quality_label(0.85) == "Excellent"
        assert score_to_quality_label(0.95) == "Excellent"
        assert score_to_quality_label(0.81) == "Excellent"

    def test_good_label(self):
        assert score_to_quality_label(0.6) == "Good"
        assert score_to_quality_label(0.75) == "Good"
        assert score_to_quality_label(0.51) == "Good"

    def test_poor_shape_label(self):
        assert score_to_quality_label(0.5) == "Poor Shape"
        assert score_to_quality_label(0.3) == "Poor Shape"
        assert score_to_quality_label(0.0) == "Poor Shape"

    def test_not_fitted_label(self):
        assert score_to_quality_label(0.5, fitted=False) == "Not Fitted"
        assert score_to_quality_label(0.9, fitted=False) == "Not Fitted"

    def test_boundary_values(self):
        # 0.8 is NOT > 0.8, but IS > 0.5, so it's "Good"
        assert score_to_quality_label(0.8) == "Good"
        assert score_to_quality_label(0.800001) == "Excellent"
        # 0.5 is NOT > 0.5, so it's "Poor Shape"
        assert score_to_quality_label(0.5) == "Poor Shape"
        assert score_to_quality_label(0.500001) == "Good"
