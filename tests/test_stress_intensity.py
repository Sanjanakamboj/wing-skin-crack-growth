"""F-L. Mode-I stress intensity and stress-intensity range."""

from __future__ import annotations

import math

import pytest

from crackgrowth import (
    StressCycle,
    ThroughCrackGeometry,
    delta_stress_intensity,
    stress_intensity,
)


# F. K hand calculation
def test_stress_intensity_hand_calculation(unit_geometry):
    # K = 1.0 * 100e6 * sqrt(pi * 0.001)
    #   = 1e8 * sqrt(3.14159265e-3) = 1e8 * 0.056049912 = 5.6049912e6 Pa*sqrt(m)
    k = stress_intensity(100.0e6, 1.0e-3, unit_geometry)
    assert k == pytest.approx(5.604991216397929e6, rel=1e-12)


def test_stress_intensity_with_non_unit_geometry_factor():
    geometry = ThroughCrackGeometry(geometry_factor=1.12)
    # 1.12 * 200e6 * sqrt(pi * 0.004) = 1.12 * 2e8 * 0.112099824 = 2.5110361e7
    k = stress_intensity(200.0e6, 4.0e-3, geometry)
    assert k == pytest.approx(1.12 * 2.0e8 * math.sqrt(math.pi * 4.0e-3), rel=1e-12)
    assert k == pytest.approx(2.5110360e7, rel=1e-6)


def test_stress_intensity_is_signed_for_compression(unit_geometry):
    tensile = stress_intensity(100.0e6, 2.0e-3, unit_geometry)
    compressive = stress_intensity(-100.0e6, 2.0e-3, unit_geometry)
    assert tensile > 0.0
    assert compressive < 0.0
    assert compressive == pytest.approx(-tensile)


def test_zero_stress_gives_zero_stress_intensity(unit_geometry):
    assert stress_intensity(0.0, 1.0e-3, unit_geometry) == 0.0


# G. K proportional to sigma
def test_stress_intensity_scales_linearly_with_stress(unit_geometry):
    base = stress_intensity(50.0e6, 3.0e-3, unit_geometry)
    tripled = stress_intensity(150.0e6, 3.0e-3, unit_geometry)
    assert tripled == pytest.approx(3.0 * base)


# H. K proportional to sqrt(a)
def test_stress_intensity_scales_with_sqrt_of_crack_length(unit_geometry):
    base = stress_intensity(100.0e6, 1.0e-3, unit_geometry)
    quadrupled_length = stress_intensity(100.0e6, 4.0e-3, unit_geometry)
    assert quadrupled_length == pytest.approx(2.0 * base)


@pytest.mark.parametrize("bad", [0.0, -1.0e-3])
def test_non_positive_crack_length_rejected(unit_geometry, bad):
    with pytest.raises(ValueError, match="strictly positive"):
        stress_intensity(100.0e6, bad, unit_geometry)


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_non_finite_inputs_rejected(unit_geometry, bad):
    with pytest.raises(ValueError, match="finite"):
        stress_intensity(bad, 1.0e-3, unit_geometry)
    with pytest.raises(ValueError, match="finite"):
        stress_intensity(100.0e6, bad, unit_geometry)


# I. delta-K hand calculation
def test_delta_k_hand_calculation(unit_geometry, cycle_100mpa_range):
    # delta_K = 1.0 * 100e6 * sqrt(pi * 0.001) = 5.604991216e6 Pa*sqrt(m)
    dk = delta_stress_intensity(1.0e-3, cycle_100mpa_range, unit_geometry)
    assert dk == pytest.approx(5.604991216397929e6, rel=1e-12)
    assert dk / 1.0e6 == pytest.approx(5.604991216, rel=1e-9)  # MPa*sqrt(m)


def test_delta_k_at_target_crack_length(unit_geometry, cycle_100mpa_range):
    # delta_K = 1e8 * sqrt(pi * 0.01) = 1e8 * 0.1772453851 = 1.772453851e7
    dk = delta_stress_intensity(10.0e-3, cycle_100mpa_range, unit_geometry)
    assert dk == pytest.approx(1.7724538509055160e7, rel=1e-12)


def test_delta_k_equals_k_difference_over_the_cycle(unit_geometry):
    """delta_K is the algebraic difference of the signed endpoint intensities."""
    cycle = StressCycle(sigma_max=100.0e6, sigma_min=-40.0e6)
    a = 2.5e-3
    k_max = stress_intensity(cycle.sigma_max, a, unit_geometry)
    k_min = stress_intensity(cycle.sigma_min, a, unit_geometry)
    assert delta_stress_intensity(a, cycle, unit_geometry) == pytest.approx(
        k_max - k_min
    )


def test_compressive_minimum_still_contributes_to_delta_k(unit_geometry):
    """Explicit: Milestone 1 models no closure, so compression inflates delta_K."""
    tension_only = StressCycle(sigma_max=100.0e6, sigma_min=0.0)
    with_compression = StressCycle(sigma_max=100.0e6, sigma_min=-100.0e6)
    a = 1.0e-3
    assert delta_stress_intensity(a, with_compression, unit_geometry) == pytest.approx(
        2.0 * delta_stress_intensity(a, tension_only, unit_geometry)
    )


# J. delta-K proportional to delta-sigma
def test_delta_k_scales_linearly_with_stress_range(unit_geometry):
    small = StressCycle(sigma_max=50.0e6, sigma_min=0.0)
    large = StressCycle(sigma_max=200.0e6, sigma_min=0.0)
    a = 2.0e-3
    assert delta_stress_intensity(a, large, unit_geometry) == pytest.approx(
        4.0 * delta_stress_intensity(a, small, unit_geometry)
    )


# K. delta-K proportional to Y
def test_delta_k_scales_linearly_with_geometry_factor(cycle_100mpa_range):
    a = 2.0e-3
    base = delta_stress_intensity(a, cycle_100mpa_range, ThroughCrackGeometry(1.0))
    scaled = delta_stress_intensity(a, cycle_100mpa_range, ThroughCrackGeometry(1.5))
    assert scaled == pytest.approx(1.5 * base)


# L. delta-K proportional to sqrt(a)
def test_delta_k_scales_with_sqrt_of_crack_length(
    unit_geometry, cycle_100mpa_range
):
    base = delta_stress_intensity(1.0e-3, cycle_100mpa_range, unit_geometry)
    nine_times = delta_stress_intensity(9.0e-3, cycle_100mpa_range, unit_geometry)
    assert nine_times == pytest.approx(3.0 * base)


def test_zero_stress_range_gives_zero_delta_k(unit_geometry):
    cycle = StressCycle(sigma_max=80.0e6, sigma_min=80.0e6)
    assert delta_stress_intensity(1.0e-3, cycle, unit_geometry) == 0.0


def test_delta_k_rejects_wrong_types(unit_geometry, cycle_100mpa_range):
    with pytest.raises(TypeError):
        delta_stress_intensity(1.0e-3, "not a cycle", unit_geometry)
    with pytest.raises(TypeError):
        delta_stress_intensity(1.0e-3, cycle_100mpa_range, 1.0)
