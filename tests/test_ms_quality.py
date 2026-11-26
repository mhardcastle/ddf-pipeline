"""Unit tests for the pure-numpy metric functions in scripts/ms_quality.py.

These do not need python-casacore: the metric arithmetic is tested on synthetic
arrays, so it runs on a laptop and in CI. The casacore IO path is exercised
separately (casa-backed tests are added where the dependency exists).
"""
import numpy as np
import pytest

from scripts import ms_quality as mq


def test_residual_rms_known_ratio():
    data = np.ones((2, 1, 2), dtype=complex)
    model = 0.9 * np.ones((2, 1, 2), dtype=complex)
    out = mq.residual_rms(data, model)
    assert out["rms_abs"] == pytest.approx([0.1, 0.1])
    assert out["rms_over_data"] == pytest.approx([0.1, 0.1])


def test_residual_rms_ignores_flagged():
    data = np.ones((2, 1, 1), dtype=complex)
    model = np.zeros_like(data)
    data[0] = 10.0  # flagged outlier must not enter the RMS
    flag = np.zeros((2, 1, 1), dtype=bool)
    flag[0] = True
    out = mq.residual_rms(data, model, flag)
    assert out["rms_abs"][0] == pytest.approx(1.0)


def test_normalised_residual_scales_with_sqrt_weight():
    data = np.array([[[1 + 1j]], [[1 - 1j]]])
    model = np.zeros_like(data)
    base = mq.normalised_residual(data, model, np.ones((2, 1, 1)))[0]
    assert base == pytest.approx(1.0)
    scaled = mq.normalised_residual(data, model, 4.0 * np.ones((2, 1, 1)))[0]
    assert scaled == pytest.approx(2.0)


def test_flag_fraction_per_corr_and_worst_channel():
    flag = np.array([[[True, False]], [[True, True]]])
    out = mq.flag_fraction(flag)
    assert out["total"] == pytest.approx(0.75)
    assert out["per_corr"] == pytest.approx([1.0, 0.5])
    assert out["worst_channel"] == pytest.approx(0.75)


def test_weight_stats_nonpositive_and_nonfinite():
    out = mq.weight_stats(np.array([0.0, 1.0, 2.0, -1.0, np.nan]))
    assert out["frac_nonpositive"] == pytest.approx(0.5)
    assert out["frac_nonfinite"] == pytest.approx(0.2)
    assert out["mean"] == pytest.approx(0.5)
    assert out["min"] == pytest.approx(-1.0)
    assert out["max"] == pytest.approx(2.0)


def test_solution_stats_unit_gains():
    G = np.ones((2, 1, 2, 1), dtype=complex)
    out = mq.solution_stats(G)
    assert out["amp_median"] == pytest.approx(1.0)
    assert out["amp_rms"] == pytest.approx(0.0)
    assert out["coverage"] == pytest.approx(1.0)
    assert out["frac_nonfinite"] == pytest.approx(0.0)
    assert out["phase_rms_deg"] == pytest.approx(0.0)


def test_solution_stats_counts_nan_gains():
    G = np.ones((2, 1, 2, 1), dtype=complex)
    G[0, 0, 0, 0] = np.nan
    out = mq.solution_stats(G)
    assert out["frac_nonfinite"] == pytest.approx(0.25)
    assert out["coverage"] == pytest.approx(0.75)


def test_solution_stats_circular_phase_spread():
    G = np.ones((2, 1, 1, 1), dtype=complex)
    G[1] = np.exp(1j * np.pi / 2)
    out = mq.solution_stats(G)
    # circular std of phases {0, pi/2}
    expected = np.degrees(np.sqrt(-2.0 * np.log(np.sqrt(2.0) / 2.0)))
    assert out["phase_rms_deg"] == pytest.approx(expected, abs=0.1)


def test_solution_stats_jones_matrix_amplitude():
    eye = np.eye(2, dtype=complex)
    G = np.ones((2, 1, 1, 1, 2, 2), dtype=complex) * eye
    out = mq.solution_stats(G)
    assert out["amp_median"] == pytest.approx(1.0)


def test_model_data_correlation():
    data = np.array([[[1.0]], [[1.0]]], dtype=complex)
    same = np.array([[[1.0]], [[1.0]]], dtype=complex)
    anti = np.array([[[1.0]], [[-1.0]]], dtype=complex)
    assert mq.model_data_correlation(data, same)[0] == pytest.approx(1.0)
    assert mq.model_data_correlation(data, anti)[0] == pytest.approx(0.0)


def test_main_returns_2_on_missing_path():
    assert mq.main(["ms_quality.py", "/nonexistent/for/sure"]) == 2
