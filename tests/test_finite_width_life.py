"""U, V, AF-AP. Life to the finite-width fracture boundary."""

from __future__ import annotations

import math

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GEOMETRY,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CriticalCrackStatus,
    FiniteWidthCenterCrack,
    FlawAdmissibility,
    FractureToughness,
    InitialFlawBeyondCriticalError,
    StressCycle,
    cycles_to_critical_crack,
    cycles_to_finite_width_fracture,
    delta_k_for_geometry,
    finite_width_critical_crack_length,
    frozen_geometry_factor_life,
    mpa_sqrt_m_to_pa_sqrt_m,
)

PANEL = CANONICAL_FINITE_WIDTH_GEOMETRY
LAW = CANONICAL_PARIS_LAW
K30 = CANONICAL_FRACTURE_TOUGHNESS


def _life(a0=1.0e-3, cycle=CANONICAL_CYCLE, geometry=PANEL, toughness=K30, **kwargs):
    return cycles_to_finite_width_fracture(
        a0, cycle, geometry, LAW, toughness, **kwargs
    )


# --------------------------------------------------------------------------
# Canonical result, fields independently verified
# --------------------------------------------------------------------------
def test_canonical_finite_width_result():
    result = _life()
    assert result.admissibility is FlawAdmissibility.BELOW_CRITICAL
    assert result.plate_width == 100.0e-3
    assert result.critical_crack_length == pytest.approx(
        0.017093952819407475, rel=1e-9
    )
    assert result.predicted_cycles == pytest.approx(842070.5528120178, rel=1e-9)
    assert result.initial_geometry_factor == pytest.approx(
        1.0002468111606977, rel=1e-12
    )
    assert result.ligament_fraction_at_critical == pytest.approx(
        0.6581209436118505, rel=1e-9
    )
    assert result.integration_method == "composite Simpson on a log(a) grid"

    growth = result.growth_result
    assert growth is not None
    assert growth.final_geometry_factor == pytest.approx(1.078807164149301, rel=1e-9)
    assert growth.initial_delta_k == pytest.approx(5606374.590785748, rel=1e-9)
    assert growth.initial_growth_rate == pytest.approx(1.762164044350064e-09, rel=1e-9)
    assert growth.final_growth_rate == pytest.approx(1.5625e-07, rel=1e-9)
    assert result.delta_k_at_critical == pytest.approx(25.0e6, rel=1e-9)

    assert result.critical_assessment is not None
    assert result.critical_assessment.k_max == pytest.approx(30.0e6, rel=1e-9)
    assert result.critical_assessment.utilization == pytest.approx(1.0, rel=1e-9)
    assert result.initial_assessment.passes


# --------------------------------------------------------------------------
# U. finite-width life is shorter than the infinite-plate life
# --------------------------------------------------------------------------
def test_finite_width_life_is_shorter_than_infinite_plate():
    finite = _life().predicted_cycles
    infinite = cycles_to_critical_crack(
        CANONICAL_INITIAL_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        LAW,
        K30,
    ).analytical_cycles
    assert finite < infinite
    assert finite / infinite == pytest.approx(0.9556, rel=1e-3)  # ~4.4 % shorter


def test_canonical_reductions_relative_to_milestone_two():
    finite = _life()
    infinite = cycles_to_critical_crack(
        CANONICAL_INITIAL_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        LAW,
        K30,
    )
    size_reduction = 1.0 - finite.critical_crack_length / infinite.critical_crack_length
    life_reduction = 1.0 - finite.predicted_cycles / infinite.analytical_cycles
    assert size_reduction == pytest.approx(0.1408, rel=1e-3)
    assert life_reduction == pytest.approx(0.0444, rel=1e-2)
    # the critical size moves much more than the life does, because most of the
    # life is spent while the crack is small and Y is still close to 1
    assert size_reduction > 3.0 * life_reduction


# --------------------------------------------------------------------------
# V. constant-Y limit as the panel widens
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "width, size_tol, life_tol",
    [(1.0, 1e-2, 1e-3), (10.0, 1e-4, 1e-5), (100.0, 1e-6, 1e-7)],
)
def test_wide_panel_converges_to_the_milestone_two_result(width, size_tol, life_tol):
    infinite = cycles_to_critical_crack(
        1.0e-3, CANONICAL_CYCLE, CANONICAL_GEOMETRY, LAW, K30
    )
    finite = _life(geometry=FiniteWidthCenterCrack(plate_width=width))
    assert finite.critical_crack_length == pytest.approx(
        infinite.critical_crack_length, rel=size_tol
    )
    assert finite.predicted_cycles == pytest.approx(
        infinite.analytical_cycles, rel=life_tol
    )
    assert finite.initial_geometry_factor == pytest.approx(1.0, rel=size_tol)


def test_convergence_to_infinite_plate_is_monotonic_in_width():
    infinite = cycles_to_critical_crack(
        1.0e-3, CANONICAL_CYCLE, CANONICAL_GEOMETRY, LAW, K30
    ).analytical_cycles
    errors = [
        abs(_life(geometry=FiniteWidthCenterCrack(plate_width=w)).predicted_cycles
            - infinite)
        for w in (0.1, 0.2, 0.5, 1.0, 5.0)
    ]
    assert all(later < earlier for earlier, later in zip(errors, errors[1:]))


# --------------------------------------------------------------------------
# AO/AP. admissibility edge cases
# --------------------------------------------------------------------------
def test_initial_flaw_at_critical_gives_zero_life():
    a_c = finite_width_critical_crack_length(
        CANONICAL_CYCLE, PANEL, K30
    ).critical_crack_length
    result = _life(a0=a_c)
    assert result.admissibility is FlawAdmissibility.AT_CRITICAL
    assert result.predicted_cycles == 0.0
    assert result.growth_result is None


def test_initial_flaw_above_critical_is_reported_honestly():
    a_c = finite_width_critical_crack_length(
        CANONICAL_CYCLE, PANEL, K30
    ).critical_crack_length
    result = _life(a0=a_c * 1.2)
    assert result.admissibility is FlawAdmissibility.ABOVE_CRITICAL
    assert result.predicted_cycles == 0.0
    assert result.growth_result is None
    assert not result.initial_assessment.passes


def test_initial_flaw_above_critical_raises_in_strict_mode():
    a_c = finite_width_critical_crack_length(
        CANONICAL_CYCLE, PANEL, K30
    ).critical_crack_length
    with pytest.raises(InitialFlawBeyondCriticalError, match="finite-width"):
        _life(a0=a_c * 1.2, strict=True)


def test_marginally_admissible_flaw_gives_small_positive_life():
    a_c = finite_width_critical_crack_length(
        CANONICAL_CYCLE, PANEL, K30
    ).critical_crack_length
    lives = [_life(a0=a_c * f).predicted_cycles for f in (1 - 1e-3, 1 - 1e-6, 1 - 1e-9)]
    assert all(v > 0.0 and math.isfinite(v) for v in lives)
    assert all(later < earlier for earlier, later in zip(lives, lives[1:]))


def test_initial_flaw_beyond_half_width_rejected():
    with pytest.raises(ValueError, match="not admissible"):
        _life(a0=0.06)


# --------------------------------------------------------------------------
# AN. zero stress range and non-tensile cycles
# --------------------------------------------------------------------------
def test_zero_stress_range_gives_infinite_life():
    cycle = StressCycle(sigma_max=120.0e6, sigma_min=120.0e6)
    result = _life(cycle=cycle)
    assert result.admissibility is FlawAdmissibility.BELOW_CRITICAL
    assert math.isfinite(result.critical_crack_length)
    assert result.predicted_cycles == math.inf


def test_non_tensile_cycle_has_no_boundary_and_infinite_life():
    cycle = StressCycle(sigma_max=-10.0e6, sigma_min=-100.0e6)
    result = _life(cycle=cycle)
    assert result.admissibility is FlawAdmissibility.NO_TENSILE_BOUNDARY
    assert result.critical.status is CriticalCrackStatus.NO_TENSILE_BOUNDARY
    assert result.predicted_cycles == math.inf
    assert result.growth_result is None
    assert result.delta_k_at_critical is None


# --------------------------------------------------------------------------
# AF. frozen-Y diagnostic
# --------------------------------------------------------------------------
def test_frozen_geometry_factor_overestimates_life():
    result = _life()
    frozen = frozen_geometry_factor_life(
        1.0e-3, result.critical_crack_length, CANONICAL_CYCLE, PANEL, LAW
    )
    assert frozen > result.predicted_cycles
    assert frozen / result.predicted_cycles == pytest.approx(1.0218, rel=1e-3)


def test_frozen_geometry_factor_error_grows_for_a_narrower_panel():
    """The approximation gets worse where finite-width effects are stronger."""
    errors = []
    for width in (0.2, 0.1, 0.05):
        panel = FiniteWidthCenterCrack(plate_width=width)
        result = _life(geometry=panel)
        frozen = frozen_geometry_factor_life(
            1.0e-3, result.critical_crack_length, CANONICAL_CYCLE, panel, LAW
        )
        errors.append(frozen / result.predicted_cycles - 1.0)
    assert all(later > earlier for earlier, later in zip(errors, errors[1:]))


# --------------------------------------------------------------------------
# Determinism and convergence
# --------------------------------------------------------------------------
def test_repeated_calculation_is_bit_for_bit_identical():
    first = _life(intervals=512)
    for _ in range(5):
        assert _life(intervals=512) == first


def test_life_converges_under_interval_refinement():
    reference = _life(intervals=20000).predicted_cycles
    errors = [
        abs(_life(intervals=n).predicted_cycles - reference) / reference
        for n in (50, 100, 200, 400)
    ]
    assert all(later < earlier for earlier, later in zip(errors, errors[1:]))
    assert errors[0] < 1.0e-8  # already very accurate at 50 intervals


def test_measured_convergence_order_is_fourth():
    reference = _life(intervals=100000).predicted_cycles

    def err(n):
        return abs(_life(intervals=n).predicted_cycles - reference) / reference

    for coarse in (50, 100, 200):
        assert 13.0 < err(coarse) / err(2 * coarse) < 19.0, coarse


def test_rejects_bad_types():
    with pytest.raises(TypeError):
        cycles_to_finite_width_fracture(
            1.0e-3, CANONICAL_CYCLE, CANONICAL_GEOMETRY, LAW, K30
        )
    with pytest.raises(TypeError):
        cycles_to_finite_width_fracture(1.0e-3, CANONICAL_CYCLE, PANEL, "law", K30)
    with pytest.raises(ValueError, match="even"):
        _life(intervals=101)


def test_solver_options_are_passed_through():
    result = _life(tolerance=1.0e-6)
    assert result.critical.tolerance == 1.0e-6


def test_delta_k_at_critical_ratio_holds_for_several_toughnesses():
    for k_ic_mpa in (20.0, 30.0, 45.0):
        toughness = FractureToughness(
            name="t",
            k_ic=mpa_sqrt_m_to_pa_sqrt_m(k_ic_mpa),
            source_note="ILLUSTRATIVE (test)",
        )
        result = _life(toughness=toughness)
        assert result.delta_k_at_critical / toughness.k_ic == pytest.approx(
            5.0 / 6.0, rel=1e-9
        )
        assert result.delta_k_at_critical == pytest.approx(
            delta_k_for_geometry(
                result.critical_crack_length, CANONICAL_CYCLE, PANEL
            ),
            rel=1e-12,
        )
