"""O-Y, AQ. Finite-width fracture boundary, solver, residual strength."""

from __future__ import annotations

import math

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CriticalCrackStatus,
    FiniteWidthCenterCrack,
    FractureToughness,
    StressCycle,
    ThroughCrackGeometry,
    assess_fracture_for_geometry,
    critical_crack_length,
    delta_k_for_geometry,
    finite_width_critical_crack_length,
    k_max_for_geometry,
    mpa_sqrt_m_to_pa_sqrt_m,
    residual_strength,
    residual_strength_for_geometry,
)

W = 0.1
PANEL = CANONICAL_FINITE_WIDTH_GEOMETRY
UNIT = ThroughCrackGeometry(geometry_factor=1.0)
K30 = CANONICAL_FRACTURE_TOUGHNESS


def _solve(cycle=CANONICAL_CYCLE, geometry=PANEL, toughness=K30, **kwargs):
    return finite_width_critical_crack_length(cycle, geometry, toughness, **kwargs)


# --------------------------------------------------------------------------
# O. root is correctly bracketed and located
# --------------------------------------------------------------------------
def test_critical_crack_length_canonical_value():
    result = _solve()
    assert result.status is CriticalCrackStatus.ROOT_FOUND
    assert result.has_finite_boundary
    assert result.critical_crack_length == pytest.approx(
        0.017093952819407475, rel=1e-9
    )


def test_root_is_bracketed_by_a_sign_change():
    """Independent bracket check: f < 0 just below a_c and f > 0 just above."""
    a_c = _solve().critical_crack_length
    f = lambda a: k_max_for_geometry(a, CANONICAL_CYCLE, PANEL) - K30.k_ic
    assert f(a_c * 0.99) < 0.0
    assert f(a_c * 1.01) > 0.0


def test_critical_result_diagnostics_are_consistent():
    result = _solve()
    a_c = result.critical_crack_length
    assert result.plate_width == W
    assert result.total_crack_length_at_critical == pytest.approx(2.0 * a_c)
    assert result.residual_ligament == pytest.approx(W - 2.0 * a_c)
    assert result.ligament_fraction == pytest.approx(1.0 - 2.0 * a_c / W)
    assert result.ligament_fraction == pytest.approx(0.6581209436118505, rel=1e-9)
    assert result.geometry_factor_at_critical == pytest.approx(
        1.078807164149301, rel=1e-9
    )
    assert result.iterations > 0
    assert result.tolerance > 0.0
    assert 0.0 < result.ligament_fraction < 1.0


# --------------------------------------------------------------------------
# P. K_max at the critical size equals K_IC
# --------------------------------------------------------------------------
@pytest.mark.parametrize("k_ic_mpa", [20.0, 30.0, 45.0])
@pytest.mark.parametrize("width", [0.05, 0.1, 0.2])
def test_k_max_at_critical_equals_toughness(k_ic_mpa, width):
    toughness = FractureToughness(
        name="t",
        k_ic=mpa_sqrt_m_to_pa_sqrt_m(k_ic_mpa),
        source_note="ILLUSTRATIVE (test)",
    )
    panel = FiniteWidthCenterCrack(plate_width=width)
    result = _solve(geometry=panel, toughness=toughness)
    assert result.status is CriticalCrackStatus.ROOT_FOUND
    assert k_max_for_geometry(
        result.critical_crack_length, CANONICAL_CYCLE, panel
    ) == pytest.approx(toughness.k_ic, rel=1e-9)
    assert result.k_max_at_critical == pytest.approx(toughness.k_ic, rel=1e-9)


def test_constant_y_formula_would_give_the_wrong_answer():
    """Guard: the M2 closed form must NOT be used for a finite-width panel."""
    a_c = _solve().critical_crack_length
    m2 = critical_crack_length(CANONICAL_CYCLE, UNIT, K30).critical_crack_length
    assert a_c < m2
    assert a_c != pytest.approx(m2, rel=1e-3)


# --------------------------------------------------------------------------
# Q. residual strength at a_c returns sigma_max
# --------------------------------------------------------------------------
@pytest.mark.parametrize("width", [0.05, 0.1, 0.5])
@pytest.mark.parametrize("sigma_max", [80.0e6, 120.0e6, 200.0e6])
def test_residual_strength_at_critical_returns_sigma_max(width, sigma_max):
    panel = FiniteWidthCenterCrack(plate_width=width)
    cycle = StressCycle(sigma_max=sigma_max, sigma_min=0.0)
    result = _solve(cycle=cycle, geometry=panel)
    assert result.status is CriticalCrackStatus.ROOT_FOUND
    assert residual_strength_for_geometry(
        result.critical_crack_length, panel, K30
    ) == pytest.approx(sigma_max, rel=1e-9)


# --------------------------------------------------------------------------
# R. ΔK(a_c)/K_IC = Δσ/σ_max survives Y(a)
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "sigma_max, sigma_min", [(120.0e6, 20.0e6), (200.0e6, 0.0), (150.0e6, -50.0e6)]
)
@pytest.mark.parametrize("width", [0.05, 0.1, 0.3])
def test_delta_k_at_critical_is_the_stress_share_of_toughness(
    sigma_max, sigma_min, width
):
    """Y(a_c) multiplies both delta_K and K_max, so it cancels in the ratio."""
    cycle = StressCycle(sigma_max=sigma_max, sigma_min=sigma_min)
    panel = FiniteWidthCenterCrack(plate_width=width)
    result = _solve(cycle=cycle, geometry=panel)
    dk = delta_k_for_geometry(result.critical_crack_length, cycle, panel)
    assert dk / K30.k_ic == pytest.approx(
        cycle.stress_range / cycle.sigma_max, rel=1e-9
    )


def test_canonical_delta_k_at_critical_is_five_sixths_of_toughness():
    a_c = _solve().critical_crack_length
    dk = delta_k_for_geometry(a_c, CANONICAL_CYCLE, PANEL)
    assert dk / K30.k_ic == pytest.approx(5.0 / 6.0, rel=1e-9)
    assert dk == pytest.approx(25.0e6, rel=1e-9)


# --------------------------------------------------------------------------
# S/T. finite width is more severe than the infinite plate
# --------------------------------------------------------------------------
def test_finite_width_critical_size_is_smaller_than_infinite_plate():
    a_c = _solve().critical_crack_length
    m2 = critical_crack_length(CANONICAL_CYCLE, UNIT, K30).critical_crack_length
    assert a_c < m2
    assert a_c / m2 == pytest.approx(0.8592, rel=1e-3)  # ~14 % reduction


@pytest.mark.parametrize("a", [1.0e-3, 5.0e-3, 20.0e-3, 40.0e-3])
def test_finite_width_residual_strength_is_below_infinite_plate(a):
    assert residual_strength_for_geometry(a, PANEL, K30) < residual_strength(
        a, UNIT, K30
    )


def test_residual_strength_hand_calculation_and_monotonicity():
    # sigma_res(20 mm) = 30e6 / (sqrt(pi*0.02) * 1.1117859405028423)
    #                  = 1.0764903544858584e8 Pa
    assert residual_strength_for_geometry(20.0e-3, PANEL, K30) == pytest.approx(
        1.0764903544858584e8, rel=1e-12
    )
    values = [
        residual_strength_for_geometry(a, PANEL, K30)
        for a in (0.5e-3, 1.0e-3, 5.0e-3, 20.0e-3, 40.0e-3, 49.0e-3)
    ]
    assert all(later < earlier for earlier, later in zip(values, values[1:]))


def test_residual_strength_falls_toward_zero_near_half_width():
    assert residual_strength_for_geometry(0.4999 * W, PANEL, K30) < 2.0e6
    assert residual_strength_for_geometry(0.499999 * W, PANEL, K30) < 3.0e5


# --------------------------------------------------------------------------
# W/X/Y. solver determinism, bounds, and never touching W/2
# --------------------------------------------------------------------------
def test_solver_is_deterministic():
    first = _solve()
    for _ in range(5):
        assert _solve() == first


def test_solver_respects_a_tighter_tolerance():
    coarse = _solve(tolerance=1.0e-6)
    fine = _solve(tolerance=1.0e-14)
    assert fine.iterations > coarse.iterations
    assert fine.critical_crack_length == pytest.approx(
        coarse.critical_crack_length, abs=1.0e-6
    )
    assert fine.tolerance == 1.0e-14


def test_solver_never_evaluates_at_or_beyond_half_width():
    """Y diverges at W/2. The panel rejects a >= W/2, so if the solver ever
    probed there the solve would raise rather than return."""
    result = _solve(tolerance=1.0e-15, max_iterations=400)
    assert result.status is CriticalCrackStatus.ROOT_FOUND
    assert result.critical_crack_length < PANEL.max_admissible_crack_length
    assert math.isfinite(result.geometry_factor_at_critical)


def test_solver_result_stays_inside_the_bracket():
    result = _solve(lower_bound=1.0e-6)
    assert 1.0e-6 <= result.critical_crack_length < 0.5 * W


def test_solver_rejects_a_lower_bound_above_the_upper_bracket():
    with pytest.raises(ValueError, match="not below the admissible upper bracket"):
        _solve(lower_bound=0.06)


def test_solver_reports_when_lower_bound_is_already_critical():
    """A very tough-limited case: K_max already exceeds K_IC at the bound."""
    weak = FractureToughness(
        name="weak", k_ic=1.0e5, source_note="ILLUSTRATIVE (test)"
    )
    result = _solve(toughness=weak, lower_bound=1.0e-3)
    assert result.status is CriticalCrackStatus.ALREADY_CRITICAL_AT_LOWER_BOUND
    assert result.critical_crack_length == math.inf
    assert not result.has_finite_boundary


def test_solver_reports_no_root_when_toughness_is_unreachable():
    """With Y -> infinity at W/2 a root normally exists, but the status is
    implemented honestly rather than assumed impossible."""
    enormous = FractureToughness(
        name="huge", k_ic=1.0e15, source_note="ILLUSTRATIVE (test)"
    )
    result = _solve(toughness=enormous)
    assert result.status is CriticalCrackStatus.NO_ROOT_IN_BRACKET
    assert result.critical_crack_length == math.inf


# --------------------------------------------------------------------------
# AQ. non-tensile sigma_max
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "cycle",
    [
        StressCycle(sigma_max=-10.0e6, sigma_min=-100.0e6),
        StressCycle(sigma_max=0.0, sigma_min=-100.0e6),
    ],
)
def test_no_tensile_boundary_for_non_opening_cycle(cycle):
    result = _solve(cycle=cycle)
    assert result.status is CriticalCrackStatus.NO_TENSILE_BOUNDARY
    assert result.critical_crack_length == math.inf
    assert result.iterations == 0


def test_compression_is_not_mirrored_into_a_boundary():
    tensile = _solve(cycle=StressCycle(sigma_max=120.0e6, sigma_min=0.0))
    compressive = _solve(
        cycle=StressCycle(sigma_max=-120.0e6, sigma_min=-200.0e6)
    )
    assert tensile.has_finite_boundary
    assert compressive.critical_crack_length == math.inf


# --------------------------------------------------------------------------
# Geometry-aware assessment
# --------------------------------------------------------------------------
def test_assessment_fields_are_consistent():
    a = 20.0e-3
    assessment = assess_fracture_for_geometry(a, CANONICAL_CYCLE, PANEL, K30)
    assert assessment.geometry_factor == pytest.approx(
        1.1117859405028423, rel=1e-12
    )
    assert assessment.k_max == pytest.approx(3.3442008885619722e7, rel=1e-12)
    assert assessment.utilization == pytest.approx(
        assessment.k_max / assessment.k_ic, rel=1e-12
    )
    assert assessment.margin == pytest.approx(
        assessment.k_ic / assessment.k_max - 1.0, rel=1e-12
    )
    assert assessment.utilization == pytest.approx(
        assessment.applied_sigma_max / assessment.residual_strength, rel=1e-12
    )
    assert assessment.ligament_fraction == pytest.approx(0.6)
    # a = 20 mm is beyond the 17.09 mm boundary
    assert not assessment.passes


def test_assessment_passes_below_and_fails_above_the_boundary():
    a_c = _solve().critical_crack_length
    below = assess_fracture_for_geometry(a_c * 0.99, CANONICAL_CYCLE, PANEL, K30)
    above = assess_fracture_for_geometry(a_c * 1.01, CANONICAL_CYCLE, PANEL, K30)
    assert below.passes and below.margin > 0.0 and below.utilization < 1.0
    assert not above.passes and above.margin < 0.0 and above.utilization > 1.0


def test_assessment_ligament_fraction_is_none_for_constant_geometry():
    assessment = assess_fracture_for_geometry(10.0e-3, CANONICAL_CYCLE, UNIT, K30)
    assert assessment.ligament_fraction is None
    assert assessment.geometry_factor == 1.0


def test_ligament_fraction_is_not_blended_into_the_margin():
    """The margin depends on K_max/K_IC alone, never on the ligament."""
    a = 10.0e-3
    narrow = assess_fracture_for_geometry(
        a, CANONICAL_CYCLE, FiniteWidthCenterCrack(0.05), K30
    )
    wide = assess_fracture_for_geometry(
        a, CANONICAL_CYCLE, FiniteWidthCenterCrack(0.5), K30
    )
    assert narrow.ligament_fraction != pytest.approx(wide.ligament_fraction)
    for assessment in (narrow, wide):
        assert assessment.margin == pytest.approx(
            assessment.k_ic / assessment.k_max - 1.0, rel=1e-12
        )


def test_non_opening_assessment_is_non_governing():
    cycle = StressCycle(sigma_max=-10.0e6, sigma_min=-100.0e6)
    assessment = assess_fracture_for_geometry(10.0e-3, cycle, PANEL, K30)
    assert assessment.k_max < 0.0
    assert assessment.utilization == 0.0
    assert assessment.margin == math.inf
    assert assessment.passes


def test_solver_rejects_bad_types_and_values():
    with pytest.raises(TypeError):
        finite_width_critical_crack_length("cycle", PANEL, K30)
    with pytest.raises(TypeError):
        finite_width_critical_crack_length(CANONICAL_CYCLE, UNIT, K30)
    with pytest.raises(TypeError):
        finite_width_critical_crack_length(CANONICAL_CYCLE, PANEL, "toughness")
    with pytest.raises(ValueError, match="strictly positive"):
        _solve(tolerance=0.0)
    with pytest.raises(ValueError, match="at least 1"):
        _solve(max_iterations=0)
