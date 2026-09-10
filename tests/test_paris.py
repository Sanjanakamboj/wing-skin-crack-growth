"""M-S. Paris-law validation, unit conversion, and growth-rate behaviour."""

from __future__ import annotations

import dataclasses
import math

import pytest

from crackgrowth import (
    ILLUSTRATIVE_ALUMINIUM_LIKE_PARIS,
    ILLUSTRATIVE_DISCLAIMER,
    ParisLaw,
    StressCycle,
    ThroughCrackGeometry,
    crack_growth_rate,
    crack_growth_rate_point,
    delta_stress_intensity,
    paris_c_from_mpa_basis,
    paris_c_to_mpa_basis,
    paris_law_from_mpa_basis,
)


# M. Paris-law validation
@pytest.mark.parametrize("bad_c", [0.0, -1.0e-30])
def test_non_positive_c_rejected(bad_c):
    with pytest.raises(ValueError, match="strictly positive"):
        ParisLaw(C=bad_c, m=3.0, name="x", source_note="test")


@pytest.mark.parametrize("bad_m", [0.0, -1.0])
def test_non_positive_m_rejected(bad_m):
    with pytest.raises(ValueError, match="strictly positive"):
        ParisLaw(C=1.0e-29, m=bad_m, name="x", source_note="test")


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_non_finite_paris_parameters_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        ParisLaw(C=bad, m=3.0, name="x", source_note="test")
    with pytest.raises(ValueError, match="finite"):
        ParisLaw(C=1.0e-29, m=bad, name="x", source_note="test")


def test_empty_name_or_source_note_rejected():
    with pytest.raises(ValueError, match="name"):
        ParisLaw(C=1.0e-29, m=3.0, name="  ", source_note="test")
    with pytest.raises(ValueError, match="source_note"):
        ParisLaw(C=1.0e-29, m=3.0, name="x", source_note="")


def test_paris_law_is_immutable(simple_paris):
    with pytest.raises(dataclasses.FrozenInstanceError):
        simple_paris.C = 1.0


def test_canonical_curve_is_flagged_illustrative():
    assert ILLUSTRATIVE_DISCLAIMER in ILLUSTRATIVE_ALUMINIUM_LIKE_PARIS.source_note


# N. MPa*sqrt(m) -> Pa*sqrt(m) conversion, verified independently
@pytest.mark.parametrize(
    "c_mpa, m, expected_c_si",
    [
        (1.0e-11, 3.0, 1.0e-11 / 1.0e18),  # (1e6)**3 = 1e18 -> 1e-29
        (1.0e-12, 4.0, 1.0e-12 / 1.0e24),  # (1e6)**4 = 1e24 -> 1e-36
        (2.0e-10, 2.0, 2.0e-10 / 1.0e12),  # (1e6)**2 = 1e12 -> 2e-22
    ],
)
def test_c_conversion_hand_calculation(c_mpa, m, expected_c_si):
    assert paris_c_from_mpa_basis(c_mpa, m) == pytest.approx(expected_c_si, rel=1e-12)


def test_c_conversion_for_non_integer_exponent():
    # C_SI = 1e-12 / (1e6)**3.5 = 1e-12 / 1e21 = 1e-33
    assert paris_c_from_mpa_basis(1.0e-12, 3.5) == pytest.approx(1.0e-33, rel=1e-12)


def test_c_conversion_round_trip():
    c_mpa, m = 3.7e-11, 2.9
    c_si = paris_c_from_mpa_basis(c_mpa, m)
    assert paris_c_to_mpa_basis(c_si, m) == pytest.approx(c_mpa, rel=1e-12)


def test_canonical_curve_reports_its_mpa_basis_coefficient():
    assert ILLUSTRATIVE_ALUMINIUM_LIKE_PARIS.C == pytest.approx(1.0e-29, rel=1e-12)
    assert ILLUSTRATIVE_ALUMINIUM_LIKE_PARIS.c_on_mpa_basis == pytest.approx(
        1.0e-11, rel=1e-12
    )


def test_conversion_rejects_invalid_input():
    with pytest.raises(ValueError):
        paris_c_from_mpa_basis(0.0, 3.0)
    with pytest.raises(ValueError):
        paris_c_from_mpa_basis(1.0e-11, -3.0)


# O. converted coefficient reproduces an identical da/dN
@pytest.mark.parametrize("m", [2.0, 3.0, 3.5, 4.2])
@pytest.mark.parametrize("delta_k_mpa", [3.0, 5.604991216, 17.724538509, 25.0])
def test_converted_coefficient_reproduces_handbook_growth_rate(m, delta_k_mpa):
    """da/dN from the SI curve must equal C_MPa * (delta_K in MPa*sqrt(m))**m."""
    c_mpa = 1.0e-11
    law = paris_law_from_mpa_basis(
        c_mpa, m, name="converted", source_note="ILLUSTRATIVE (test)"
    )
    handbook_rate = c_mpa * delta_k_mpa**m
    si_rate = law.growth_rate(delta_k_mpa * 1.0e6)
    assert si_rate == pytest.approx(handbook_rate, rel=1e-12)


def test_mismatched_unit_basis_is_grossly_wrong():
    """Guard: using an MPa-basis C with an SI delta_K is off by (1e6)**m."""
    m = 3.0
    c_mpa = 1.0e-11
    wrong = ParisLaw(C=c_mpa, m=m, name="wrong", source_note="ILLUSTRATIVE (test)")
    right = paris_law_from_mpa_basis(
        c_mpa, m, name="right", source_note="ILLUSTRATIVE (test)"
    )
    dk_si = 5.0e6
    assert wrong.growth_rate(dk_si) == pytest.approx(
        right.growth_rate(dk_si) * (1.0e6) ** m, rel=1e-12
    )


# P. da/dN hand calculation
def test_growth_rate_hand_calculation(unit_geometry, cycle_100mpa_range, simple_paris):
    # delta_K(1 mm) = 5.604991216e6 Pa*sqrt(m)
    # da/dN = 1e-29 * (5.604991216e6)**3
    #       = 1e-29 * 1.76085992e20 = 1.76085992e-9 m/cycle
    rate = crack_growth_rate(1.0e-3, cycle_100mpa_range, unit_geometry, simple_paris)
    assert rate == pytest.approx(1.7608599228871055e-09, rel=1e-10)


def test_growth_rate_at_target_hand_calculation(
    unit_geometry, cycle_100mpa_range, simple_paris
):
    # delta_K(10 mm) = 1.7724538509e7; da/dN = 1e-29 * (1.7724538509e7)**3
    rate = crack_growth_rate(10.0e-3, cycle_100mpa_range, unit_geometry, simple_paris)
    assert rate == pytest.approx(5.568327996831707e-08, rel=1e-10)


# Q. monotonic with a
def test_growth_rate_increases_monotonically_with_crack_length(
    unit_geometry, cycle_100mpa_range, simple_paris
):
    lengths = [0.5e-3, 1.0e-3, 2.0e-3, 5.0e-3, 10.0e-3, 20.0e-3]
    rates = [
        crack_growth_rate(a, cycle_100mpa_range, unit_geometry, simple_paris)
        for a in lengths
    ]
    assert all(later > earlier for earlier, later in zip(rates, rates[1:]))


# R. scaling with C
def test_growth_rate_scales_linearly_with_c(unit_geometry, cycle_100mpa_range):
    base = ParisLaw(C=1.0e-29, m=3.0, name="a", source_note="ILLUSTRATIVE (test)")
    tenfold = ParisLaw(C=1.0e-28, m=3.0, name="b", source_note="ILLUSTRATIVE (test)")
    a = 2.0e-3
    assert crack_growth_rate(
        a, cycle_100mpa_range, unit_geometry, tenfold
    ) == pytest.approx(
        10.0 * crack_growth_rate(a, cycle_100mpa_range, unit_geometry, base)
    )


# S. scaling with delta_K ** m
@pytest.mark.parametrize("m", [1.5, 2.0, 3.0, 4.0])
def test_growth_rate_scales_as_delta_k_to_the_m(unit_geometry, m):
    law = ParisLaw(C=1.0e-30, m=m, name="c", source_note="ILLUSTRATIVE (test)")
    small = StressCycle(sigma_max=50.0e6, sigma_min=0.0)
    large = StressCycle(sigma_max=100.0e6, sigma_min=0.0)
    a = 3.0e-3
    ratio = crack_growth_rate(a, large, unit_geometry, law) / crack_growth_rate(
        a, small, unit_geometry, law
    )
    assert ratio == pytest.approx(2.0**m, rel=1e-12)


def test_growth_rate_scales_with_geometry_factor_to_the_m(
    cycle_100mpa_range, simple_paris
):
    a = 3.0e-3
    base = crack_growth_rate(
        a, cycle_100mpa_range, ThroughCrackGeometry(1.0), simple_paris
    )
    doubled = crack_growth_rate(
        a, cycle_100mpa_range, ThroughCrackGeometry(2.0), simple_paris
    )
    assert doubled == pytest.approx(base * 2.0**simple_paris.m, rel=1e-12)


# Zero-range policy at the rate level
def test_zero_stress_range_gives_zero_growth_rate(unit_geometry, simple_paris):
    cycle = StressCycle(sigma_max=90.0e6, sigma_min=90.0e6)
    assert crack_growth_rate(1.0e-3, cycle, unit_geometry, simple_paris) == 0.0


def test_growth_rate_rejects_negative_delta_k(simple_paris):
    with pytest.raises(ValueError, match="non-negative"):
        simple_paris.growth_rate(-1.0e6)


def test_growth_rate_rejects_bad_inputs(
    unit_geometry, cycle_100mpa_range, simple_paris
):
    with pytest.raises(TypeError):
        crack_growth_rate(1.0e-3, cycle_100mpa_range, unit_geometry, "not a law")
    with pytest.raises(ValueError, match="strictly positive"):
        crack_growth_rate(0.0, cycle_100mpa_range, unit_geometry, simple_paris)
    with pytest.raises(ValueError, match="finite"):
        crack_growth_rate(math.nan, cycle_100mpa_range, unit_geometry, simple_paris)


# Structured, auditable rate point
def test_rate_point_exposes_consistent_intermediates(
    unit_geometry, cycle_100mpa_range, simple_paris
):
    a = 4.0e-3
    point = crack_growth_rate_point(
        a, cycle_100mpa_range, unit_geometry, simple_paris
    )
    assert point.crack_length == a
    assert point.delta_stress == pytest.approx(100.0e6)
    assert point.geometry_factor == 1.0
    assert point.delta_k == pytest.approx(
        delta_stress_intensity(a, cycle_100mpa_range, unit_geometry)
    )
    assert point.growth_rate == pytest.approx(
        crack_growth_rate(a, cycle_100mpa_range, unit_geometry, simple_paris)
    )
    assert point.growth_rate == pytest.approx(
        simple_paris.C * point.delta_k**simple_paris.m
    )
