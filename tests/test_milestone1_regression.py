"""AI/AJ. Milestone 1 behaviour must be preserved exactly by Milestone 2.

These are deliberate regression locks. Milestone 2 adds a fracture boundary but
must not perturb any Milestone 1 number, API, or convention. If a change to the
growth mechanics is ever made, these tests should be the first to complain.
"""

from __future__ import annotations

import math

import pytest

import crackgrowth as cg
from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_GEOMETRY,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CANONICAL_TARGET_CRACK_LENGTH,
    analytical_cycles_to_crack_length,
    crack_growth_rate,
    cycles_to_crack_length,
    delta_stress_intensity,
    stress_intensity,
)


# --------------------------------------------------------------------------
# AI. the Milestone 1 life-to-10-mm result is unchanged
# --------------------------------------------------------------------------
def test_milestone1_canonical_life_to_imposed_target_is_unchanged():
    result = cycles_to_crack_length(
        CANONICAL_INITIAL_CRACK_LENGTH,
        CANONICAL_TARGET_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    )
    analytical = analytical_cycles_to_crack_length(
        CANONICAL_INITIAL_CRACK_LENGTH,
        CANONICAL_TARGET_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    )
    assert result.predicted_cycles == pytest.approx(776634.4447, abs=1e-3)
    assert analytical == pytest.approx(776634.4444503565, rel=1e-12)
    assert result.integration_intervals == 1000


def test_milestone1_imposed_target_is_still_ten_millimetres():
    """Milestone 2 must not repoint the M1 constant at the critical size."""
    assert CANONICAL_TARGET_CRACK_LENGTH == 10.0e-3


def test_milestone1_canonical_inputs_are_unchanged():
    assert CANONICAL_GEOMETRY.geometry_factor == 1.0
    assert CANONICAL_CYCLE.sigma_max == 120.0e6
    assert CANONICAL_CYCLE.sigma_min == 20.0e6
    assert CANONICAL_CYCLE.stress_range == pytest.approx(100.0e6)
    assert CANONICAL_CYCLE.stress_ratio == pytest.approx(1.0 / 6.0)
    assert CANONICAL_INITIAL_CRACK_LENGTH == 1.0e-3
    assert CANONICAL_PARIS_LAW.C == pytest.approx(1.0e-29, rel=1e-12)
    assert CANONICAL_PARIS_LAW.m == 3.0
    assert CANONICAL_PARIS_LAW.c_on_mpa_basis == pytest.approx(1.0e-11, rel=1e-12)


def test_milestone1_endpoint_values_are_unchanged():
    a0, af = CANONICAL_INITIAL_CRACK_LENGTH, CANONICAL_TARGET_CRACK_LENGTH
    assert delta_stress_intensity(
        a0, CANONICAL_CYCLE, CANONICAL_GEOMETRY
    ) == pytest.approx(5.604991216397929e6, rel=1e-12)
    assert delta_stress_intensity(
        af, CANONICAL_CYCLE, CANONICAL_GEOMETRY
    ) == pytest.approx(1.7724538509055160e7, rel=1e-12)
    assert crack_growth_rate(
        a0, CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_PARIS_LAW
    ) == pytest.approx(1.7608599228871055e-09, rel=1e-10)
    assert crack_growth_rate(
        af, CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_PARIS_LAW
    ) == pytest.approx(5.568327996831707e-08, rel=1e-10)


# --------------------------------------------------------------------------
# AJ. the Milestone 1 API and conventions still behave identically
# --------------------------------------------------------------------------
def test_milestone1_public_api_is_still_exported():
    for name in (
        "ThroughCrackGeometry",
        "INFINITE_PLATE_THROUGH_CRACK",
        "StressCycle",
        "stress_intensity",
        "delta_stress_intensity",
        "ParisLaw",
        "crack_growth_rate",
        "crack_growth_rate_point",
        "paris_c_from_mpa_basis",
        "paris_c_to_mpa_basis",
        "paris_law_from_mpa_basis",
        "MPA_ROOT_M_IN_PA_ROOT_M",
        "ILLUSTRATIVE_DISCLAIMER",
        "ILLUSTRATIVE_ALUMINIUM_LIKE_PARIS",
        "CrackGrowthResult",
        "cycles_to_crack_length",
        "analytical_cycles_to_crack_length",
        "DEFAULT_INTEGRATION_INTERVALS",
        "stress_range_sensitivity",
        "initial_crack_length_sensitivity",
    ):
        assert name in cg.__all__, name
        assert hasattr(cg, name), name


def test_milestone1_signed_stress_intensity_convention_is_unchanged():
    tensile = stress_intensity(100.0e6, 2.0e-3, CANONICAL_GEOMETRY)
    compressive = stress_intensity(-100.0e6, 2.0e-3, CANONICAL_GEOMETRY)
    assert tensile > 0.0 > compressive
    assert compressive == pytest.approx(-tensile)


def test_milestone1_zero_range_no_growth_policy_is_unchanged():
    cycle = cg.StressCycle(sigma_max=80.0e6, sigma_min=80.0e6)
    result = cycles_to_crack_length(
        1.0e-3, 10.0e-3, cycle, CANONICAL_GEOMETRY, CANONICAL_PARIS_LAW
    )
    assert result.predicted_cycles == math.inf
    assert (
        analytical_cycles_to_crack_length(
            1.0e-3, 10.0e-3, cycle, CANONICAL_GEOMETRY, CANONICAL_PARIS_LAW
        )
        == math.inf
    )


def test_milestone1_integrator_still_requires_increasing_limits():
    with pytest.raises(ValueError, match="a_final must be strictly greater"):
        cycles_to_crack_length(
            10.0e-3, 1.0e-3, CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_PARIS_LAW
        )


def test_milestone1_integrator_still_requires_even_intervals():
    with pytest.raises(ValueError, match="even"):
        cycles_to_crack_length(
            1.0e-3,
            10.0e-3,
            CANONICAL_CYCLE,
            CANONICAL_GEOMETRY,
            CANONICAL_PARIS_LAW,
            intervals=101,
        )


def test_milestone1_stress_range_sensitivity_is_unchanged():
    points = cg.stress_range_sensitivity(
        (0.5, 0.75, 1.0, 1.25, 1.5),
        1.0e-3,
        CANONICAL_TARGET_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    )
    lives = [p.numerical_cycles for p in points]
    assert lives[2] == pytest.approx(776634.4447, abs=1e-3)
    reference = lives[2]
    for point, life in zip(points, lives):
        assert life == pytest.approx(
            reference * point.parameter_value**-CANONICAL_PARIS_LAW.m, rel=1e-9
        )


def test_milestone2_does_not_change_the_growth_model_used_at_ten_millimetres():
    """The M2 helper integrating past 10 mm must pass through exactly the same
    growth rates the M1 model gives there."""
    a = 10.0e-3
    m1_rate = crack_growth_rate(
        a, CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_PARIS_LAW
    )
    result = cg.cycles_to_critical_crack(
        a,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
        cg.CANONICAL_FRACTURE_TOUGHNESS,
    )
    assert result.growth_result is not None
    assert result.growth_result.initial_growth_rate == pytest.approx(m1_rate)


def test_milestone2_life_to_critical_exceeds_milestone1_life_to_ten_mm():
    """a_c = 19.894 mm lies beyond the imposed 10 mm target, so the
    toughness-derived life must be the larger of the two, and the difference is
    exactly the M1 life from 10 mm to a_c."""
    m1 = analytical_cycles_to_crack_length(
        1.0e-3,
        CANONICAL_TARGET_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    )
    m2 = cg.cycles_to_critical_crack(
        1.0e-3,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
        cg.CANONICAL_FRACTURE_TOUGHNESS,
    )
    assert m2.analytical_cycles > m1
    remainder = analytical_cycles_to_crack_length(
        CANONICAL_TARGET_CRACK_LENGTH,
        m2.critical_crack_length,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    )
    assert m1 + remainder == pytest.approx(m2.analytical_cycles, rel=1e-12)
