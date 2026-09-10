"""A-T, AK-AM. Fracture toughness, critical size, residual strength, margin."""

from __future__ import annotations

import dataclasses
import math

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GEOMETRY,
    ILLUSTRATIVE_TOUGHNESS_DISCLAIMER,
    FractureBoundaryStatus,
    FractureToughness,
    StressCycle,
    ThroughCrackGeometry,
    assess_fracture,
    critical_crack_length,
    delta_stress_intensity,
    fracture_margin,
    fracture_utilization,
    k_max,
    mpa_sqrt_m_to_pa_sqrt_m,
    pa_sqrt_m_to_mpa_sqrt_m,
    residual_strength,
)

UNIT = ThroughCrackGeometry(geometry_factor=1.0)
K30 = FractureToughness(
    name="test", k_ic=30.0e6, source_note="ILLUSTRATIVE (test fixture)"
)


# --------------------------------------------------------------------------
# A. toughness validation
# --------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [0.0, -1.0e6])
def test_non_positive_toughness_rejected(bad):
    with pytest.raises(ValueError, match="strictly positive"):
        FractureToughness(name="x", k_ic=bad, source_note="ILLUSTRATIVE (test)")


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_non_finite_toughness_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        FractureToughness(name="x", k_ic=bad, source_note="ILLUSTRATIVE (test)")


def test_toughness_is_immutable():
    with pytest.raises(dataclasses.FrozenInstanceError):
        K30.k_ic = 1.0


def test_empty_name_rejected():
    with pytest.raises(ValueError, match="name"):
        FractureToughness(name="  ", k_ic=30.0e6, source_note="ILLUSTRATIVE (test)")


# --------------------------------------------------------------------------
# B. toughness provenance
# --------------------------------------------------------------------------
def test_empty_provenance_rejected():
    with pytest.raises(ValueError, match="source_note"):
        FractureToughness(name="x", k_ic=30.0e6, source_note="")


def test_canonical_toughness_is_flagged_illustrative():
    assert (
        ILLUSTRATIVE_TOUGHNESS_DISCLAIMER
        in CANONICAL_FRACTURE_TOUGHNESS.source_note
    )


def test_canonical_toughness_value():
    assert CANONICAL_FRACTURE_TOUGHNESS.k_ic == pytest.approx(30.0e6, rel=1e-12)
    assert CANONICAL_FRACTURE_TOUGHNESS.k_ic_in_mpa_sqrt_m == pytest.approx(
        30.0, rel=1e-12
    )


# --------------------------------------------------------------------------
# C. MPa*sqrt(m) -> Pa*sqrt(m) conversion
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "mpa, pa", [(1.0, 1.0e6), (30.0, 30.0e6), (0.5, 5.0e5), (72.5, 7.25e7)]
)
def test_toughness_unit_conversion_hand_calculation(mpa, pa):
    assert mpa_sqrt_m_to_pa_sqrt_m(mpa) == pytest.approx(pa, rel=1e-12)
    assert pa_sqrt_m_to_mpa_sqrt_m(pa) == pytest.approx(mpa, rel=1e-12)


def test_toughness_unit_conversion_round_trip():
    assert pa_sqrt_m_to_mpa_sqrt_m(mpa_sqrt_m_to_pa_sqrt_m(37.4)) == pytest.approx(
        37.4, rel=1e-12
    )


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_toughness_conversion_rejects_non_finite(bad):
    with pytest.raises(ValueError, match="finite"):
        mpa_sqrt_m_to_pa_sqrt_m(bad)


# --------------------------------------------------------------------------
# D. K_max hand calculation
# --------------------------------------------------------------------------
def test_k_max_hand_calculation():
    # K_max = 1.0 * 120e6 * sqrt(pi * 0.001) = 1.2e8 * 0.0560499121639793
    #       = 6.725989459677515e6 Pa*sqrt(m)
    assert k_max(1.0e-3, CANONICAL_CYCLE, UNIT) == pytest.approx(
        6.725989459677515e6, rel=1e-12
    )


def test_k_max_scales_with_sqrt_a_and_sigma_max():
    base = k_max(1.0e-3, CANONICAL_CYCLE, UNIT)
    assert k_max(4.0e-3, CANONICAL_CYCLE, UNIT) == pytest.approx(2.0 * base)


# --------------------------------------------------------------------------
# E. fracture uses sigma_max, NOT delta_sigma
# --------------------------------------------------------------------------
def test_k_max_uses_sigma_max_not_stress_range():
    """Two cycles with identical delta_sigma but different sigma_max must give
    different K_max -- and identical delta_K."""
    low = StressCycle(sigma_max=120.0e6, sigma_min=20.0e6)
    high = StressCycle(sigma_max=180.0e6, sigma_min=80.0e6)
    assert low.stress_range == high.stress_range
    a = 2.0e-3
    assert delta_stress_intensity(a, low, UNIT) == pytest.approx(
        delta_stress_intensity(a, high, UNIT)
    )
    assert k_max(a, high, UNIT) == pytest.approx(1.5 * k_max(a, low, UNIT))
    assert k_max(a, low, UNIT) != pytest.approx(delta_stress_intensity(a, low, UNIT))


def test_critical_size_uses_sigma_max_not_stress_range():
    """If delta_sigma were used as the fracture criterion, these two cycles
    would share a critical size. They must not."""
    low = StressCycle(sigma_max=120.0e6, sigma_min=20.0e6)
    high = StressCycle(sigma_max=180.0e6, sigma_min=80.0e6)
    a_low = critical_crack_length(low, UNIT, K30).critical_crack_length
    a_high = critical_crack_length(high, UNIT, K30).critical_crack_length
    assert a_high < a_low
    assert a_high == pytest.approx(a_low / 1.5**2, rel=1e-12)


def test_critical_size_is_not_the_delta_sigma_based_value():
    """Guard against silently substituting delta_sigma for sigma_max."""
    a_c = critical_crack_length(CANONICAL_CYCLE, UNIT, K30).critical_crack_length
    wrong = (K30.k_ic / (1.0 * CANONICAL_CYCLE.stress_range)) ** 2 / math.pi
    assert a_c != pytest.approx(wrong)
    assert a_c == pytest.approx(wrong * (100.0 / 120.0) ** 2, rel=1e-12)


# --------------------------------------------------------------------------
# F. critical-size hand calculation
# --------------------------------------------------------------------------
def test_critical_crack_length_hand_calculation():
    # a_c = (1/pi) * (30e6 / (1.0 * 120e6))**2 = (1/pi) * 0.0625
    #     = 0.019894367886486918 m = 19.894 mm
    result = critical_crack_length(CANONICAL_CYCLE, UNIT, K30)
    assert result.critical_crack_length == pytest.approx(
        0.019894367886486918, rel=1e-12
    )
    assert result.status is FractureBoundaryStatus.FINITE_CRITICAL_SIZE
    assert result.has_finite_boundary
    assert result.sigma_max == 120.0e6
    assert result.geometry_factor == 1.0
    assert result.k_ic == 30.0e6


# --------------------------------------------------------------------------
# G/H/I. exact scaling of the critical size
# --------------------------------------------------------------------------
@pytest.mark.parametrize("factor", [0.5, 2.0, 3.0])
def test_critical_size_scales_as_k_ic_squared(factor):
    base = critical_crack_length(CANONICAL_CYCLE, UNIT, K30).critical_crack_length
    scaled = critical_crack_length(
        CANONICAL_CYCLE,
        UNIT,
        FractureToughness(
            name="s", k_ic=30.0e6 * factor, source_note="ILLUSTRATIVE (test)"
        ),
    ).critical_crack_length
    assert scaled == pytest.approx(base * factor**2, rel=1e-12)


@pytest.mark.parametrize("factor", [0.5, 1.5, 2.0])
def test_critical_size_scales_as_sigma_max_to_the_minus_two(factor):
    base = critical_crack_length(
        StressCycle(sigma_max=120.0e6, sigma_min=0.0), UNIT, K30
    ).critical_crack_length
    scaled = critical_crack_length(
        StressCycle(sigma_max=120.0e6 * factor, sigma_min=0.0), UNIT, K30
    ).critical_crack_length
    assert scaled == pytest.approx(base * factor**-2, rel=1e-12)


@pytest.mark.parametrize("factor", [0.8, 1.2, 2.0])
def test_critical_size_scales_as_y_to_the_minus_two(factor):
    base = critical_crack_length(CANONICAL_CYCLE, UNIT, K30).critical_crack_length
    scaled = critical_crack_length(
        CANONICAL_CYCLE, ThroughCrackGeometry(factor), K30
    ).critical_crack_length
    assert scaled == pytest.approx(base * factor**-2, rel=1e-12)


# --------------------------------------------------------------------------
# J. compression-only: no tensile mode-I fracture boundary
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "cycle",
    [
        StressCycle(sigma_max=-10.0e6, sigma_min=-100.0e6),
        StressCycle(sigma_max=0.0, sigma_min=-100.0e6),
        StressCycle(sigma_max=0.0, sigma_min=0.0),
    ],
)
def test_no_tensile_boundary_for_non_opening_cycle(cycle):
    result = critical_crack_length(cycle, UNIT, K30)
    assert result.status is FractureBoundaryStatus.NO_TENSILE_BOUNDARY
    assert result.critical_crack_length == math.inf
    assert not result.has_finite_boundary


def test_compression_does_not_produce_a_mirrored_critical_size():
    """abs(sigma_max) must never be taken: -120 MPa is not equivalent to +120."""
    tensile = critical_crack_length(
        StressCycle(sigma_max=120.0e6, sigma_min=0.0), UNIT, K30
    )
    compressive = critical_crack_length(
        StressCycle(sigma_max=-120.0e6, sigma_min=-200.0e6), UNIT, K30
    )
    assert math.isfinite(tensile.critical_crack_length)
    assert compressive.critical_crack_length == math.inf


def test_compression_only_assessment_is_non_governing():
    cycle = StressCycle(sigma_max=-10.0e6, sigma_min=-100.0e6)
    assessment = assess_fracture(1.0e-3, cycle, UNIT, K30)
    assert assessment.k_max < 0.0
    assert assessment.utilization == 0.0
    assert assessment.margin == math.inf
    assert assessment.passes
    assert assessment.status is FractureBoundaryStatus.NO_TENSILE_BOUNDARY


# --------------------------------------------------------------------------
# K. residual-strength hand calculation
# --------------------------------------------------------------------------
def test_residual_strength_hand_calculation():
    # sigma_res(1 mm) = 30e6 / (1.0 * sqrt(pi * 0.001))
    #                 = 30e6 / 0.0560499121639793 = 5.3523723484583133e8 Pa
    assert residual_strength(1.0e-3, UNIT, K30) == pytest.approx(
        5.3523723484583133e8, rel=1e-12
    )


# --------------------------------------------------------------------------
# L/M/N. residual-strength scaling
# --------------------------------------------------------------------------
def test_residual_strength_scales_as_inverse_sqrt_a():
    base = residual_strength(1.0e-3, UNIT, K30)
    assert residual_strength(4.0e-3, UNIT, K30) == pytest.approx(base / 2.0, rel=1e-12)
    assert residual_strength(9.0e-3, UNIT, K30) == pytest.approx(base / 3.0, rel=1e-12)


def test_residual_strength_decreases_monotonically_with_crack_size():
    values = [
        residual_strength(a, UNIT, K30)
        for a in (0.5e-3, 1.0e-3, 2.0e-3, 5.0e-3, 10.0e-3, 20.0e-3)
    ]
    assert all(later < earlier for earlier, later in zip(values, values[1:]))


def test_residual_strength_scales_linearly_with_toughness():
    base = residual_strength(2.0e-3, UNIT, K30)
    doubled = residual_strength(
        2.0e-3,
        UNIT,
        FractureToughness(name="s", k_ic=60.0e6, source_note="ILLUSTRATIVE (test)"),
    )
    assert doubled == pytest.approx(2.0 * base, rel=1e-12)


def test_residual_strength_scales_inversely_with_geometry_factor():
    base = residual_strength(2.0e-3, UNIT, K30)
    stiffer = residual_strength(2.0e-3, ThroughCrackGeometry(2.0), K30)
    assert stiffer == pytest.approx(base / 2.0, rel=1e-12)


# --------------------------------------------------------------------------
# O. critical-size / residual-strength round trip
# --------------------------------------------------------------------------
@pytest.mark.parametrize("sigma_max", [40.0e6, 120.0e6, 250.0e6])
@pytest.mark.parametrize("k_ic_mpa", [20.0, 30.0, 55.0])
@pytest.mark.parametrize("y", [0.8, 1.0, 1.5])
def test_residual_strength_at_critical_size_returns_sigma_max(sigma_max, k_ic_mpa, y):
    """The single most important consistency check of this module: the critical
    size and the residual strength are the same equation solved two ways."""
    cycle = StressCycle(sigma_max=sigma_max, sigma_min=0.0)
    geometry = ThroughCrackGeometry(y)
    toughness = FractureToughness(
        name="rt",
        k_ic=mpa_sqrt_m_to_pa_sqrt_m(k_ic_mpa),
        source_note="ILLUSTRATIVE (test)",
    )
    a_c = critical_crack_length(cycle, geometry, toughness).critical_crack_length
    assert residual_strength(a_c, geometry, toughness) == pytest.approx(
        sigma_max, rel=1e-12
    )


# --------------------------------------------------------------------------
# AB. K_max at the critical size equals K_IC
# --------------------------------------------------------------------------
@pytest.mark.parametrize("k_ic_mpa", [20.0, 30.0, 45.0])
def test_k_max_at_critical_size_equals_toughness(k_ic_mpa):
    toughness = FractureToughness(
        name="b",
        k_ic=mpa_sqrt_m_to_pa_sqrt_m(k_ic_mpa),
        source_note="ILLUSTRATIVE (test)",
    )
    a_c = critical_crack_length(
        CANONICAL_CYCLE, UNIT, toughness
    ).critical_crack_length
    assert k_max(a_c, CANONICAL_CYCLE, UNIT) == pytest.approx(
        toughness.k_ic, rel=1e-12
    )


# --------------------------------------------------------------------------
# AA. delta_K(a_c) / K_IC = delta_sigma / sigma_max
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "sigma_max, sigma_min", [(120.0e6, 20.0e6), (200.0e6, 0.0), (150.0e6, -50.0e6)]
)
@pytest.mark.parametrize("y", [0.9, 1.0, 1.3])
def test_delta_k_at_critical_is_the_stress_ratio_share_of_toughness(
    sigma_max, sigma_min, y
):
    cycle = StressCycle(sigma_max=sigma_max, sigma_min=sigma_min)
    geometry = ThroughCrackGeometry(y)
    a_c = critical_crack_length(cycle, geometry, K30).critical_crack_length
    expected = (cycle.stress_range / cycle.sigma_max) * K30.k_ic
    assert delta_stress_intensity(a_c, cycle, geometry) == pytest.approx(
        expected, rel=1e-12
    )


def test_canonical_delta_k_at_critical_is_five_sixths_of_toughness():
    """delta_sigma / sigma_max = 100/120 = 5/6, so delta_K(a_c) = 25 MPa*sqrt(m)."""
    a_c = critical_crack_length(
        CANONICAL_CYCLE, UNIT, K30
    ).critical_crack_length
    dk = delta_stress_intensity(a_c, CANONICAL_CYCLE, UNIT)
    assert dk / K30.k_ic == pytest.approx(5.0 / 6.0, rel=1e-12)
    assert dk == pytest.approx(25.0e6, rel=1e-12)
    assert dk < K30.k_ic


# --------------------------------------------------------------------------
# P. fracture utilization hand calculation
# --------------------------------------------------------------------------
def test_fracture_utilization_hand_calculation():
    # K_max(1 mm) / K_IC = 6.725989459677515e6 / 30e6 = 0.22419964865591716
    k = k_max(1.0e-3, CANONICAL_CYCLE, UNIT)
    assert fracture_utilization(k, K30) == pytest.approx(
        0.22419964865591716, rel=1e-12
    )


def test_utilization_is_one_at_the_boundary():
    assert fracture_utilization(K30.k_ic, K30) == pytest.approx(1.0, rel=1e-12)


# --------------------------------------------------------------------------
# Q. fracture margin hand calculation
# --------------------------------------------------------------------------
def test_fracture_margin_hand_calculation():
    # MS_K = 30e6 / 6.725989459677515e6 - 1 = 3.460310290381928
    k = k_max(1.0e-3, CANONICAL_CYCLE, UNIT)
    assert fracture_margin(k, K30) == pytest.approx(3.460310290381928, rel=1e-12)


def test_margin_and_utilization_are_consistent():
    k = k_max(3.0e-3, CANONICAL_CYCLE, UNIT)
    assert fracture_margin(k, K30) == pytest.approx(
        1.0 / fracture_utilization(k, K30) - 1.0, rel=1e-12
    )


# --------------------------------------------------------------------------
# R/S/T. behaviour at, just below, and just above the boundary
# --------------------------------------------------------------------------
def test_exact_boundary_passes_with_zero_margin():
    a_c = critical_crack_length(CANONICAL_CYCLE, UNIT, K30).critical_crack_length
    assessment = assess_fracture(a_c, CANONICAL_CYCLE, UNIT, K30)
    assert assessment.passes
    assert assessment.utilization == pytest.approx(1.0, rel=1e-12)
    assert assessment.margin == pytest.approx(0.0, abs=1e-12)
    assert assessment.residual_strength == pytest.approx(120.0e6, rel=1e-12)


def test_just_below_boundary_passes():
    a_c = critical_crack_length(CANONICAL_CYCLE, UNIT, K30).critical_crack_length
    assessment = assess_fracture(a_c * 0.999, CANONICAL_CYCLE, UNIT, K30)
    assert assessment.passes
    assert assessment.utilization < 1.0
    assert assessment.margin > 0.0
    assert assessment.residual_strength > 120.0e6


def test_just_above_boundary_fails():
    a_c = critical_crack_length(CANONICAL_CYCLE, UNIT, K30).critical_crack_length
    assessment = assess_fracture(a_c * 1.001, CANONICAL_CYCLE, UNIT, K30)
    assert not assessment.passes
    assert assessment.utilization > 1.0
    assert assessment.margin < 0.0
    assert assessment.residual_strength < 120.0e6


# --------------------------------------------------------------------------
# AN. assessment fields independently verified
# --------------------------------------------------------------------------
def test_assessment_fields_are_independently_correct():
    a = 4.0e-3
    assessment = assess_fracture(a, CANONICAL_CYCLE, UNIT, K30)
    assert assessment.crack_length == a
    assert assessment.applied_sigma_max == 120.0e6
    assert assessment.k_ic == 30.0e6
    # K_max = 120e6 * sqrt(pi * 0.004) = 1.2e8 * 0.11209982432795862
    assert assessment.k_max == pytest.approx(1.345197891935503e7, rel=1e-12)
    # sigma_res = 30e6 / sqrt(pi * 0.004) = 2.6761861742291567e8
    assert assessment.residual_strength == pytest.approx(
        2.6761861742291567e8, rel=1e-12
    )
    assert assessment.utilization == pytest.approx(
        assessment.k_max / assessment.k_ic, rel=1e-12
    )
    assert assessment.margin == pytest.approx(
        assessment.k_ic / assessment.k_max - 1.0, rel=1e-12
    )
    assert assessment.passes is True
    # utilization and the applied/residual stress ratio must agree
    assert assessment.utilization == pytest.approx(
        assessment.applied_sigma_max / assessment.residual_strength, rel=1e-12
    )


# --------------------------------------------------------------------------
# AK-AM. input rejection
# --------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [0.0, -1.0e-3])
def test_non_positive_crack_length_rejected(bad):
    with pytest.raises(ValueError, match="strictly positive"):
        residual_strength(bad, UNIT, K30)
    with pytest.raises(ValueError, match="strictly positive"):
        assess_fracture(bad, CANONICAL_CYCLE, UNIT, K30)
    with pytest.raises(ValueError, match="strictly positive"):
        k_max(bad, CANONICAL_CYCLE, UNIT)


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_non_finite_crack_length_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        residual_strength(bad, UNIT, K30)
    with pytest.raises(ValueError, match="finite"):
        assess_fracture(bad, CANONICAL_CYCLE, UNIT, K30)


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_non_finite_k_max_rejected_by_margin_helpers(bad):
    with pytest.raises(ValueError, match="finite"):
        fracture_utilization(bad, K30)
    with pytest.raises(ValueError, match="finite"):
        fracture_margin(bad, K30)


def test_wrong_types_rejected():
    with pytest.raises(TypeError):
        critical_crack_length("cycle", UNIT, K30)
    with pytest.raises(TypeError):
        critical_crack_length(CANONICAL_CYCLE, "geometry", K30)
    with pytest.raises(TypeError):
        critical_crack_length(CANONICAL_CYCLE, UNIT, "toughness")
    with pytest.raises(TypeError):
        residual_strength(1.0e-3, UNIT, "toughness")
    with pytest.raises(TypeError):
        k_max(1.0e-3, "cycle", UNIT)
    with pytest.raises(TypeError):
        fracture_margin(1.0e6, "toughness")


def test_canonical_geometry_and_cycle_unchanged_by_milestone_two():
    assert CANONICAL_GEOMETRY.geometry_factor == 1.0
    assert CANONICAL_CYCLE.sigma_max == 120.0e6
    assert CANONICAL_CYCLE.sigma_min == 20.0e6
