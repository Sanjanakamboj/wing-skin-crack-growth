"""AZ-BD. Milestones 1-4 must survive Milestone 5 untouched."""

from __future__ import annotations

import pytest

import crackgrowth as cg
from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GEOMETRY,
    CANONICAL_GROWTH_THRESHOLD,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CANONICAL_PLATE_WIDTH,
    CANONICAL_TARGET_CRACK_LENGTH,
    GrowthState,
    analytical_cycles_to_crack_length,
    crack_growth_rate_for_geometry,
    critical_crack_length,
    cycles_to_critical_crack,
    cycles_to_finite_width_fracture,
    cycles_to_fracture_with_threshold,
    delta_k_for_geometry,
    delta_stress_intensity,
    finite_width_critical_crack_length,
    finite_width_threshold_crack_length,
    thresholded_crack_growth_rate,
)

G = CANONICAL_FINITE_WIDTH_GEOMETRY
A0 = CANONICAL_INITIAL_CRACK_LENGTH


# --------------------------------------------------------------------------
# AZ/BA. Milestones 1 and 2
# --------------------------------------------------------------------------
def test_milestone1_canonical_values_unchanged():
    assert CANONICAL_TARGET_CRACK_LENGTH == 10.0e-3
    assert delta_stress_intensity(
        A0, CANONICAL_CYCLE, CANONICAL_GEOMETRY
    ) == pytest.approx(5.604991216397929e6, rel=1e-12)
    assert analytical_cycles_to_crack_length(
        A0,
        CANONICAL_TARGET_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    ) == pytest.approx(776634.4444503565, rel=1e-12)


def test_milestone2_canonical_values_unchanged():
    assert critical_crack_length(
        CANONICAL_CYCLE, CANONICAL_GEOMETRY, CANONICAL_FRACTURE_TOUGHNESS
    ).critical_crack_length == pytest.approx(0.019894367886486918, rel=1e-12)
    assert cycles_to_critical_crack(
        A0,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
        CANONICAL_FRACTURE_TOUGHNESS,
    ).analytical_cycles == pytest.approx(881160.779753657, rel=1e-9)


# --------------------------------------------------------------------------
# BB. Milestone 3
# --------------------------------------------------------------------------
def test_milestone3_canonical_values_unchanged():
    assert CANONICAL_PLATE_WIDTH == 100.0e-3
    assert G.geometry_factor_at(A0) == pytest.approx(1.0002468111606977, rel=1e-12)
    critical = finite_width_critical_crack_length(
        CANONICAL_CYCLE, G, CANONICAL_FRACTURE_TOUGHNESS
    )
    assert critical.critical_crack_length == pytest.approx(
        0.017093952819407475, rel=1e-9
    )
    assert critical.ligament_fraction == pytest.approx(0.6581209436118505, rel=1e-9)
    assert cycles_to_finite_width_fracture(
        A0, CANONICAL_CYCLE, G, CANONICAL_PARIS_LAW, CANONICAL_FRACTURE_TOUGHNESS
    ).predicted_cycles == pytest.approx(842070.5528120178, rel=1e-9)


# --------------------------------------------------------------------------
# BC. Milestone 4
# --------------------------------------------------------------------------
def test_milestone4_canonical_values_unchanged():
    assert CANONICAL_GROWTH_THRESHOLD.delta_k_threshold == pytest.approx(
        4.0e6, rel=1e-12
    )
    assert delta_k_for_geometry(A0, CANONICAL_CYCLE, G) == pytest.approx(
        5606374.590785748, rel=1e-9
    )
    assert finite_width_threshold_crack_length(
        CANONICAL_CYCLE, G, CANONICAL_GROWTH_THRESHOLD
    ).threshold_crack_length == pytest.approx(0.0005092306461735948, rel=1e-9)

    result = cycles_to_fracture_with_threshold(
        A0,
        CANONICAL_CYCLE,
        G,
        CANONICAL_PARIS_LAW,
        CANONICAL_FRACTURE_TOUGHNESS,
        CANONICAL_GROWTH_THRESHOLD,
    )
    assert result.state is GrowthState.ACTIVE_GROWTH
    assert result.threshold_ratio == pytest.approx(1.401594, rel=1e-6)
    assert result.predicted_cycles == pytest.approx(842070.5528120178, rel=1e-9)


def test_milestone4_hard_cutoff_still_applies_only_on_its_own_path():
    """Milestone 5 must not have leaked a threshold into the plain Paris rate,
    nor removed it from the threshold-aware helper."""
    tiny = 0.1e-3
    assert delta_k_for_geometry(tiny, CANONICAL_CYCLE, G) < (
        CANONICAL_GROWTH_THRESHOLD.delta_k_threshold
    )
    assert crack_growth_rate_for_geometry(
        tiny, CANONICAL_CYCLE, G, CANONICAL_PARIS_LAW
    ) > 0.0
    assert thresholded_crack_growth_rate(
        tiny,
        CANONICAL_CYCLE,
        G,
        CANONICAL_PARIS_LAW,
        CANONICAL_GROWTH_THRESHOLD,
    ).growth_rate == 0.0


# --------------------------------------------------------------------------
# BD. prior public API intact
# --------------------------------------------------------------------------
def test_prior_public_api_is_still_exported():
    for name in (
        # M1
        "ThroughCrackGeometry", "StressCycle", "stress_intensity",
        "delta_stress_intensity", "ParisLaw", "crack_growth_rate",
        "cycles_to_crack_length", "analytical_cycles_to_crack_length",
        # M2
        "FractureToughness", "critical_crack_length", "residual_strength",
        "assess_fracture", "cycles_to_critical_crack", "mpa_sqrt_m_to_pa_sqrt_m",
        # M3
        "FiniteWidthCenterCrack", "CrackGeometry", "delta_k_for_geometry",
        "finite_width_critical_crack_length", "integrate_crack_growth_log_grid",
        "cycles_to_finite_width_fracture", "frozen_geometry_factor_life",
        # M4
        "CrackGrowthThreshold", "GrowthState", "BoundaryOrdering",
        "thresholded_crack_growth_rate", "cycles_to_fracture_with_threshold",
        "finite_width_threshold_crack_length", "threshold_crack_length_constant_y",
    ):
        assert name in cg.__all__, name
        assert hasattr(cg, name), name


def test_canonical_constants_are_not_repointed():
    assert CANONICAL_INITIAL_CRACK_LENGTH == 1.0e-3
    assert CANONICAL_TARGET_CRACK_LENGTH == 10.0e-3
    assert CANONICAL_GEOMETRY.geometry_factor == 1.0
    assert CANONICAL_PLATE_WIDTH == 100.0e-3
    assert CANONICAL_CYCLE.sigma_max == 120.0e6
    assert CANONICAL_FRACTURE_TOUGHNESS.k_ic == pytest.approx(30.0e6, rel=1e-12)
    assert CANONICAL_GROWTH_THRESHOLD.delta_k_threshold == pytest.approx(
        4.0e6, rel=1e-12
    )


def test_milestone5_did_not_change_prior_defaults():
    assert cg.DEFAULT_INTEGRATION_INTERVALS == 1000
    assert cg.DEFAULT_LOG_INTERVALS == 1000
    assert cg.DEFAULT_ROOT_TOLERANCE == 1.0e-12
