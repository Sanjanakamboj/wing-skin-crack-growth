"""AT/AU/AV/AW. Milestones 1-3 must survive Milestone 4 untouched.

M1 is locked in ``test_milestone1_regression.py`` and M2 in
``test_milestone2_regression.py``; this module locks the Milestone 3
finite-width results and re-checks the earlier headline numbers now that a
threshold path exists alongside them.
"""

from __future__ import annotations

import pytest

import crackgrowth as cg
from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GEOMETRY,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CANONICAL_PLATE_WIDTH,
    CANONICAL_TARGET_CRACK_LENGTH,
    analytical_cycles_to_crack_length,
    crack_growth_rate,
    crack_growth_rate_for_geometry,
    critical_crack_length,
    cycles_to_critical_crack,
    cycles_to_finite_width_fracture,
    delta_k_for_geometry,
    delta_stress_intensity,
    finite_width_critical_crack_length,
    frozen_geometry_factor_life,
    integrate_crack_growth_log_grid,
)

PANEL = CANONICAL_FINITE_WIDTH_GEOMETRY
A0 = CANONICAL_INITIAL_CRACK_LENGTH


# --------------------------------------------------------------------------
# AT/AU. Milestone 1 and 2 headline numbers
# --------------------------------------------------------------------------
def test_milestone1_canonical_values_unchanged():
    assert analytical_cycles_to_crack_length(
        A0,
        CANONICAL_TARGET_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    ) == pytest.approx(776634.4444503565, rel=1e-12)
    assert delta_stress_intensity(
        A0, CANONICAL_CYCLE, CANONICAL_GEOMETRY
    ) == pytest.approx(5.604991216397929e6, rel=1e-12)
    assert CANONICAL_TARGET_CRACK_LENGTH == 10.0e-3


def test_milestone2_canonical_values_unchanged():
    assert critical_crack_length(
        CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_FRACTURE_TOUGHNESS
    ).critical_crack_length == pytest.approx(0.019894367886486918, rel=1e-12)
    life = cycles_to_critical_crack(
        A0,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
        CANONICAL_FRACTURE_TOUGHNESS,
    )
    assert life.analytical_cycles == pytest.approx(881160.779753657, rel=1e-9)
    assert life.delta_k_at_critical == pytest.approx(25.0e6, rel=1e-12)


# --------------------------------------------------------------------------
# AV. Milestone 3 headline numbers
# --------------------------------------------------------------------------
def test_milestone3_canonical_values_unchanged():
    assert CANONICAL_PLATE_WIDTH == 100.0e-3
    assert PANEL.geometry_factor_at(A0) == pytest.approx(
        1.0002468111606977, rel=1e-12
    )
    assert delta_k_for_geometry(A0, CANONICAL_CYCLE, PANEL) == pytest.approx(
        5606374.590785748, rel=1e-9
    )

    critical = finite_width_critical_crack_length(
        CANONICAL_CYCLE, PANEL, CANONICAL_FRACTURE_TOUGHNESS
    )
    assert critical.critical_crack_length == pytest.approx(
        0.017093952819407475, rel=1e-9
    )
    assert critical.geometry_factor_at_critical == pytest.approx(
        1.078807164149301, rel=1e-9
    )
    assert critical.k_max_at_critical == pytest.approx(30.0e6, rel=1e-9)
    assert critical.ligament_fraction == pytest.approx(0.6581209436118505, rel=1e-9)

    life = cycles_to_finite_width_fracture(
        A0, CANONICAL_CYCLE, PANEL, CANONICAL_PARIS_LAW, CANONICAL_FRACTURE_TOUGHNESS
    )
    assert life.predicted_cycles == pytest.approx(842070.5528120178, rel=1e-9)
    assert life.delta_k_at_critical == pytest.approx(25.0e6, rel=1e-9)


def test_milestone3_frozen_y_diagnostic_unchanged():
    a_c = finite_width_critical_crack_length(
        CANONICAL_CYCLE, PANEL, CANONICAL_FRACTURE_TOUGHNESS
    ).critical_crack_length
    assert frozen_geometry_factor_life(
        A0, a_c, CANONICAL_CYCLE, PANEL, CANONICAL_PARIS_LAW
    ) == pytest.approx(860455.4329222413, rel=1e-9)


def test_milestone3_log_grid_integrator_unchanged():
    analytical = analytical_cycles_to_crack_length(
        A0, 10.0e-3, CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_PARIS_LAW
    )
    result = integrate_crack_growth_log_grid(
        A0, 10.0e-3, CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_PARIS_LAW
    )
    assert result.predicted_cycles == pytest.approx(analytical, rel=1e-12)
    assert result.integration_method == "composite Simpson on a log(a) grid"


# --------------------------------------------------------------------------
# The no-threshold growth rate is still threshold-free
# --------------------------------------------------------------------------
def test_original_growth_rate_apis_apply_no_threshold():
    """Milestone 4 must NOT have silently added a cutoff to the original path."""
    tiny = 0.05e-3
    assert delta_k_for_geometry(tiny, CANONICAL_CYCLE, PANEL) < (
        cg.CANONICAL_GROWTH_THRESHOLD.delta_k_threshold
    )
    assert crack_growth_rate_for_geometry(
        tiny, CANONICAL_CYCLE, PANEL, CANONICAL_PARIS_LAW
    ) > 0.0
    assert crack_growth_rate(
        tiny, CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_PARIS_LAW
    ) > 0.0


def test_milestone3_life_helper_still_ignores_the_threshold():
    """cycles_to_finite_width_fracture takes no threshold and applies none."""
    result = cycles_to_finite_width_fracture(
        0.1e-3,
        CANONICAL_CYCLE,
        PANEL,
        CANONICAL_PARIS_LAW,
        CANONICAL_FRACTURE_TOUGHNESS,
    )
    assert result.predicted_cycles > 0.0
    assert result.growth_result is not None


# --------------------------------------------------------------------------
# AW. prior public API intact
# --------------------------------------------------------------------------
def test_prior_public_api_is_still_exported():
    for name in (
        # M1
        "ThroughCrackGeometry", "StressCycle", "stress_intensity",
        "delta_stress_intensity", "ParisLaw", "crack_growth_rate",
        "cycles_to_crack_length", "analytical_cycles_to_crack_length",
        "DEFAULT_INTEGRATION_INTERVALS",
        # M2
        "FractureToughness", "critical_crack_length", "residual_strength",
        "assess_fracture", "FlawAdmissibility", "cycles_to_critical_crack",
        "mpa_sqrt_m_to_pa_sqrt_m",
        # M3
        "FiniteWidthCenterCrack", "CrackGeometry", "geometry_factor_at",
        "delta_k_for_geometry", "crack_growth_rate_for_geometry",
        "finite_width_critical_crack_length", "residual_strength_for_geometry",
        "integrate_crack_growth_log_grid", "cycles_to_finite_width_fracture",
        "frozen_geometry_factor_life", "width_sensitivity",
        "geometry_amplification_table", "CriticalCrackStatus",
    ):
        assert name in cg.__all__, name
        assert hasattr(cg, name), name


def test_canonical_constants_are_not_repointed():
    assert CANONICAL_INITIAL_CRACK_LENGTH == 1.0e-3
    assert CANONICAL_TARGET_CRACK_LENGTH == 10.0e-3
    assert CANONICAL_GEOMETRY.geometry_factor == 1.0
    assert CANONICAL_PLATE_WIDTH == 100.0e-3
    assert CANONICAL_FRACTURE_TOUGHNESS.k_ic == pytest.approx(30.0e6, rel=1e-12)
    assert cg.CANONICAL_GROWTH_THRESHOLD.delta_k_threshold == pytest.approx(
        4.0e6, rel=1e-12
    )
