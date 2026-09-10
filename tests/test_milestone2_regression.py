"""AT/AU/AV. Milestone 1 and 2 behaviour must survive Milestone 3 untouched.

The Milestone 1 lock lives in ``test_milestone1_regression.py``; this module
locks the Milestone 2 fracture results and re-checks the Milestone 1 headline
numbers now that a second geometry and a second integrator exist.
"""

from __future__ import annotations

import math

import pytest

import crackgrowth as cg
from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GEOMETRY,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CANONICAL_TARGET_CRACK_LENGTH,
    FlawAdmissibility,
    analytical_cycles_to_crack_length,
    assess_fracture,
    critical_crack_length,
    cycles_to_crack_length,
    cycles_to_critical_crack,
    delta_stress_intensity,
    k_max,
    residual_strength,
)


# --------------------------------------------------------------------------
# AT. Milestone 1 headline numbers
# --------------------------------------------------------------------------
def test_milestone1_canonical_values_unchanged():
    a0, af = CANONICAL_INITIAL_CRACK_LENGTH, CANONICAL_TARGET_CRACK_LENGTH
    assert a0 == 1.0e-3
    assert af == 10.0e-3
    assert CANONICAL_GEOMETRY.geometry_factor == 1.0
    assert CANONICAL_CYCLE.sigma_max == 120.0e6
    assert CANONICAL_CYCLE.sigma_min == 20.0e6
    assert CANONICAL_CYCLE.stress_range == pytest.approx(100.0e6)
    assert delta_stress_intensity(
        a0, CANONICAL_CYCLE, CANONICAL_GEOMETRY
    ) == pytest.approx(5.604991216397929e6, rel=1e-12)
    assert delta_stress_intensity(
        af, CANONICAL_CYCLE, CANONICAL_GEOMETRY
    ) == pytest.approx(1.7724538509055160e7, rel=1e-12)
    assert analytical_cycles_to_crack_length(
        a0, af, CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_PARIS_LAW
    ) == pytest.approx(776634.4444503565, rel=1e-12)
    assert cycles_to_crack_length(
        a0, af, CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_PARIS_LAW
    ).predicted_cycles == pytest.approx(776634.4447, abs=1e-3)


# --------------------------------------------------------------------------
# AU. Milestone 2 headline numbers
# --------------------------------------------------------------------------
def test_milestone2_canonical_values_unchanged():
    assert CANONICAL_FRACTURE_TOUGHNESS.k_ic == pytest.approx(30.0e6, rel=1e-12)

    a0 = CANONICAL_INITIAL_CRACK_LENGTH
    assert k_max(a0, CANONICAL_CYCLE, CANONICAL_GEOMETRY) == pytest.approx(
        6.725989459677515e6, rel=1e-12
    )
    assert residual_strength(
        a0, CANONICAL_GEOMETRY, CANONICAL_FRACTURE_TOUGHNESS
    ) == pytest.approx(5.3523723484583133e8, rel=1e-12)

    critical = critical_crack_length(
        CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_FRACTURE_TOUGHNESS
    )
    assert critical.critical_crack_length == pytest.approx(
        0.019894367886486918, rel=1e-12
    )

    assessment = assess_fracture(
        a0, CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_FRACTURE_TOUGHNESS
    )
    assert assessment.margin == pytest.approx(3.460310290381928, rel=1e-12)
    assert assessment.utilization == pytest.approx(0.22419964865591716, rel=1e-12)

    life = cycles_to_critical_crack(
        a0,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
        CANONICAL_FRACTURE_TOUGHNESS,
    )
    assert life.admissibility is FlawAdmissibility.BELOW_CRITICAL
    assert life.analytical_cycles == pytest.approx(881160.779753657, rel=1e-9)
    assert life.predicted_cycles == pytest.approx(881160.7850256478, rel=1e-9)
    assert life.delta_k_at_critical == pytest.approx(25.0e6, rel=1e-12)


def test_milestone2_closed_form_critical_size_still_exact():
    """The constant-Y closed form remains available and exact as the
    infinite-plate regression reference."""
    critical = critical_crack_length(
        CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_FRACTURE_TOUGHNESS
    )
    expected = (
        CANONICAL_FRACTURE_TOUGHNESS.k_ic
        / (CANONICAL_GEOMETRY.geometry_factor * CANONICAL_CYCLE.sigma_max)
    ) ** 2 / math.pi
    assert critical.critical_crack_length == pytest.approx(expected, rel=1e-15)


# --------------------------------------------------------------------------
# AV. the prior public API is intact
# --------------------------------------------------------------------------
def test_prior_public_api_is_still_exported():
    for name in (
        # Milestone 1
        "ThroughCrackGeometry",
        "INFINITE_PLATE_THROUGH_CRACK",
        "StressCycle",
        "stress_intensity",
        "delta_stress_intensity",
        "ParisLaw",
        "crack_growth_rate",
        "paris_c_from_mpa_basis",
        "CrackGrowthResult",
        "cycles_to_crack_length",
        "analytical_cycles_to_crack_length",
        "DEFAULT_INTEGRATION_INTERVALS",
        "stress_range_sensitivity",
        "initial_crack_length_sensitivity",
        # Milestone 2
        "FractureToughness",
        "FractureBoundaryStatus",
        "CriticalCrackResult",
        "FractureAssessment",
        "mpa_sqrt_m_to_pa_sqrt_m",
        "k_max",
        "critical_crack_length",
        "residual_strength",
        "fracture_utilization",
        "fracture_margin",
        "assess_fracture",
        "FlawAdmissibility",
        "InitialFlawBeyondCriticalError",
        "cycles_to_critical_crack",
        "toughness_sensitivity",
        "residual_strength_table",
    ):
        assert name in cg.__all__, name
        assert hasattr(cg, name), name


def test_constant_geometry_dataclass_shape_is_unchanged():
    """geometry_factor_at was ADDED; the scalar attribute still works."""
    geometry = cg.ThroughCrackGeometry(geometry_factor=1.12)
    assert geometry.geometry_factor == 1.12
    assert geometry.geometry_factor_at(1.0e-3) == 1.12
    assert geometry.geometry_factor_at(50.0e-3) == 1.12


def test_milestone1_integrator_defaults_are_not_repointed():
    """Milestone 3 must not silently re-default M1 to the log grid."""
    assert cg.DEFAULT_INTEGRATION_INTERVALS == 1000
    result = cycles_to_crack_length(
        1.0e-3,
        10.0e-3,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    )
    assert type(result).__name__ == "CrackGrowthResult"
    assert result.integration_intervals == 1000


def test_milestone3_did_not_repoint_the_canonical_constants():
    assert CANONICAL_TARGET_CRACK_LENGTH == 10.0e-3
    assert CANONICAL_INITIAL_CRACK_LENGTH == 1.0e-3
    assert CANONICAL_GEOMETRY.geometry_factor == 1.0
    assert cg.CANONICAL_PLATE_WIDTH == 100.0e-3
