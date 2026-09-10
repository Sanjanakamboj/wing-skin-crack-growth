"""B-E. Stress-cycle validation and the stress-range / mean / R identities."""

from __future__ import annotations

import dataclasses
import math

import pytest

from crackgrowth import StressCycle


# B. validation
def test_min_above_max_rejected():
    with pytest.raises(ValueError, match="sigma_min must not exceed sigma_max"):
        StressCycle(sigma_max=10.0e6, sigma_min=20.0e6)


def test_equal_max_and_min_is_allowed():
    cycle = StressCycle(sigma_max=50.0e6, sigma_min=50.0e6)
    assert cycle.stress_range == 0.0


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_non_finite_stresses_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        StressCycle(sigma_max=bad, sigma_min=0.0)
    with pytest.raises(ValueError, match="finite"):
        StressCycle(sigma_max=0.0, sigma_min=bad)


def test_cycle_is_immutable():
    cycle = StressCycle(sigma_max=1.0e6, sigma_min=0.0)
    with pytest.raises(dataclasses.FrozenInstanceError):
        cycle.sigma_max = 2.0e6


# C. delta-sigma identity
def test_stress_range_identity():
    cycle = StressCycle(sigma_max=120.0e6, sigma_min=20.0e6)
    assert cycle.stress_range == pytest.approx(100.0e6)


def test_compressive_minimum_is_not_clamped():
    """No crack closure in Milestone 1: compression widens the algebraic range."""
    cycle = StressCycle(sigma_max=100.0e6, sigma_min=-50.0e6)
    assert cycle.sigma_min == -50.0e6
    assert cycle.stress_range == pytest.approx(150.0e6)


# D. mean / alternating identities
def test_mean_and_alternating_stress():
    cycle = StressCycle(sigma_max=120.0e6, sigma_min=20.0e6)
    assert cycle.mean_stress == pytest.approx(70.0e6)
    assert cycle.alternating_stress == pytest.approx(50.0e6)
    assert cycle.alternating_stress == pytest.approx(cycle.stress_range / 2.0)
    assert cycle.mean_stress + cycle.alternating_stress == pytest.approx(
        cycle.sigma_max
    )
    assert cycle.mean_stress - cycle.alternating_stress == pytest.approx(
        cycle.sigma_min
    )


def test_mean_stress_negative_for_compression_dominated_cycle():
    cycle = StressCycle(sigma_max=10.0e6, sigma_min=-90.0e6)
    assert cycle.mean_stress == pytest.approx(-40.0e6)


# E. R-ratio hand calculations
def test_stress_ratio_hand_calculation():
    # R = 20 / 120 = 1/6 = 0.1666...
    cycle = StressCycle(sigma_max=120.0e6, sigma_min=20.0e6)
    assert cycle.stress_ratio == pytest.approx(1.0 / 6.0)


def test_stress_ratio_of_fully_reversed_cycle_is_minus_one():
    cycle = StressCycle(sigma_max=80.0e6, sigma_min=-80.0e6)
    assert cycle.stress_ratio == pytest.approx(-1.0)


def test_stress_ratio_zero_to_tension_is_zero():
    assert StressCycle(sigma_max=100.0e6, sigma_min=0.0).stress_ratio == 0.0


def test_stress_ratio_undefined_when_sigma_max_is_zero():
    cycle = StressCycle(sigma_max=0.0, sigma_min=-100.0e6)
    with pytest.raises(ZeroDivisionError, match="undefined"):
        cycle.stress_ratio
