"""Sensitivity sweeps: stress-range scaling and initial-crack-size trend."""

from __future__ import annotations

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_GEOMETRY,
    CANONICAL_PARIS_LAW,
    CANONICAL_TARGET_CRACK_LENGTH,
    initial_crack_length_sensitivity,
    stress_range_sensitivity,
)

SCALES = (0.5, 0.75, 1.0, 1.25, 1.5)
INITIAL_CRACK_LENGTHS = (0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3)


def _stress_sweep():
    return stress_range_sensitivity(
        SCALES,
        1.0e-3,
        CANONICAL_TARGET_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    )


def test_stress_range_sweep_preserves_stress_ratio_and_scales_the_range():
    for point, scale in zip(_stress_sweep(), SCALES):
        assert point.parameter_value == scale
        assert point.stress_range == pytest.approx(scale * 100.0e6)


def test_life_scales_as_stress_range_to_the_minus_m():
    """N proportional to delta_sigma ** -m for the constant-Y Paris model."""
    m = CANONICAL_PARIS_LAW.m
    points = _stress_sweep()
    reference = next(p for p in points if p.parameter_value == 1.0)
    for point in points:
        expected = reference.numerical_cycles * point.parameter_value**-m
        assert point.numerical_cycles == pytest.approx(expected, rel=1e-9)


def test_stress_range_sweep_is_monotonically_decreasing():
    lives = [p.numerical_cycles for p in _stress_sweep()]
    assert all(later < earlier for earlier, later in zip(lives, lives[1:]))


def test_sweep_numerical_and_analytical_agree():
    for point in _stress_sweep():
        assert point.numerical_cycles == pytest.approx(
            point.analytical_cycles, rel=1e-8
        )


def test_initial_crack_sweep_life_decreases_with_larger_initial_crack():
    points = initial_crack_length_sensitivity(
        INITIAL_CRACK_LENGTHS,
        CANONICAL_TARGET_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    )
    assert [p.parameter_value for p in points] == list(INITIAL_CRACK_LENGTHS)
    lives = [p.numerical_cycles for p in points]
    assert all(later < earlier for earlier, later in zip(lives, lives[1:]))
    for point in points:
        assert point.numerical_cycles == pytest.approx(
            point.analytical_cycles, rel=1e-8
        )


def test_initial_crack_sweep_matches_the_closed_form_a_dependence():
    """For m = 3, N is proportional to (a0**-0.5 - af**-0.5)."""
    points = initial_crack_length_sensitivity(
        INITIAL_CRACK_LENGTHS,
        CANONICAL_TARGET_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    )
    af = CANONICAL_TARGET_CRACK_LENGTH
    reference = points[0]
    ref_shape = reference.parameter_value**-0.5 - af**-0.5
    for point in points:
        shape = point.parameter_value**-0.5 - af**-0.5
        assert point.numerical_cycles == pytest.approx(
            reference.numerical_cycles * shape / ref_shape, rel=1e-8
        )


def test_sweep_rejects_non_positive_scale():
    with pytest.raises(ValueError, match="strictly positive"):
        stress_range_sensitivity(
            (1.0, 0.0),
            1.0e-3,
            CANONICAL_TARGET_CRACK_LENGTH,
            CANONICAL_CYCLE,
            CANONICAL_GEOMETRY,
            CANONICAL_PARIS_LAW,
        )
