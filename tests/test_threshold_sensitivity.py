"""AC-AL. Threshold-aware sensitivity sweeps and their transitions."""

from __future__ import annotations

import math

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GROWTH_THRESHOLD,
    CANONICAL_PARIS_LAW,
    CrackGrowthThreshold,
    GrowthState,
    cycles_to_finite_width_fracture,
    initial_crack_sensitivity_with_threshold,
    max_stress_sensitivity_with_threshold,
    mpa_sqrt_m_to_pa_sqrt_m,
    stress_range_sensitivity_with_threshold,
    threshold_sensitivity,
    toughness_sensitivity_with_threshold,
    width_sensitivity_with_threshold,
)

PANEL = CANONICAL_FINITE_WIDTH_GEOMETRY
LAW = CANONICAL_PARIS_LAW
TOU = CANONICAL_FRACTURE_TOUGHNESS
T4 = CANONICAL_GROWTH_THRESHOLD
A0 = 1.0e-3

M3_LIFE = cycles_to_finite_width_fracture(
    A0, CANONICAL_CYCLE, PANEL, LAW, TOU
).predicted_cycles


def _threshold(mpa: float) -> CrackGrowthThreshold:
    return CrackGrowthThreshold(
        name="t",
        delta_k_threshold=mpa_sqrt_m_to_pa_sqrt_m(mpa),
        source_note="ILLUSTRATIVE (test fixture)",
    )


def _decreasing(values):
    return all(later < earlier for earlier, later in zip(values, values[1:]))


# --------------------------------------------------------------------------
# Threshold sweep: a step, not a trend
# --------------------------------------------------------------------------
THRESHOLD_SWEEP = (None, 2.0, 3.0, 4.0, 5.0, 5.5, 5.6, 5.7, 6.0, 8.0, 10.0)


def _threshold_points():
    return threshold_sensitivity(
        tuple(None if v is None else _threshold(v) for v in THRESHOLD_SWEEP),
        A0, CANONICAL_CYCLE, PANEL, LAW, TOU,
    )


def test_threshold_sweep_crosses_between_five_point_six_and_five_point_seven():
    """delta_K(a0) = 5.6064 MPa*sqrt(m), so the crossing must sit there."""
    points = {p.parameter_value: p for p in _threshold_points()}
    active = [v for v, p in points.items() if p.state is GrowthState.ACTIVE_GROWTH]
    arrested = [v for v, p in points.items() if p.state.is_arrested]
    assert max(active) == pytest.approx(5.6e6)
    assert min(arrested) == pytest.approx(5.7e6)


def test_threshold_sweep_life_takes_only_two_values():
    lives = [p.predicted_cycles for p in _threshold_points()]
    finite = [v for v in lives if math.isfinite(v)]
    assert all(v == M3_LIFE for v in finite)
    assert set(lives) == {M3_LIFE, math.inf}


def test_threshold_sweep_threshold_crack_size_grows_with_threshold():
    sizes = [
        p.threshold_crack_length
        for p in _threshold_points()
        if p.threshold_crack_length is not None
    ]
    assert all(later > earlier for earlier, later in zip(sizes, sizes[1:]))


def test_threshold_sweep_ratio_falls_with_threshold():
    ratios = [
        p.threshold_ratio for p in _threshold_points() if p.threshold_ratio is not None
    ]
    assert _decreasing(ratios)


def test_disabled_row_reports_no_threshold_fields():
    disabled = _threshold_points()[0]
    assert disabled.parameter_value == 0.0
    assert disabled.delta_k_threshold is None
    assert disabled.threshold_ratio is None
    assert disabled.threshold_crack_length is None
    assert disabled.predicted_cycles == M3_LIFE


# --------------------------------------------------------------------------
# AI/AJ. initial-crack transition
# --------------------------------------------------------------------------
A0_SWEEP = (0.1e-3, 0.25e-3, 0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3, 8.0e-3)


def _initial_points():
    return initial_crack_sensitivity_with_threshold(
        A0_SWEEP, CANONICAL_CYCLE, PANEL, LAW, TOU, T4
    )


def test_initial_crack_sweep_shows_a_real_arrest_to_growth_transition():
    points = _initial_points()
    states = [p.state for p in points]
    assert states[:3] == [GrowthState.ARRESTED_BELOW_THRESHOLD] * 3
    assert all(s is GrowthState.ACTIVE_GROWTH for s in states[3:])
    # the transition happens between 0.5 mm and 1.0 mm
    assert points[2].predicted_cycles == math.inf
    assert math.isfinite(points[3].predicted_cycles)


def test_initial_crack_sweep_ratio_rises_monotonically():
    ratios = [p.threshold_ratio for p in _initial_points()]
    assert all(later > earlier for earlier, later in zip(ratios, ratios[1:]))


def test_active_life_decreases_with_larger_initial_crack():
    lives = [
        p.predicted_cycles for p in _initial_points() if p.state.is_growing
    ]
    assert _decreasing(lives)
    assert all(math.isfinite(v) for v in lives)


def test_arrested_rows_share_one_fracture_boundary():
    points = _initial_points()
    for point in points:
        assert point.critical_crack_length == pytest.approx(
            points[0].critical_crack_length, rel=1e-9
        )


# --------------------------------------------------------------------------
# AG/AH. stress-range transition
# --------------------------------------------------------------------------
SIGMA_MAX_SWEEP = tuple(v * 1.0e6 for v in (40, 60, 80, 100, 120, 140, 160, 180))


def _stress_points():
    return stress_range_sensitivity_with_threshold(
        SIGMA_MAX_SWEEP, 20.0e6, A0, PANEL, LAW, TOU, T4
    )


def test_stress_range_sweep_transitions_from_arrest_to_growth():
    points = _stress_points()
    states = [p.state for p in points]
    assert states[:3] == [GrowthState.ARRESTED_BELOW_THRESHOLD] * 3
    assert all(s is GrowthState.ACTIVE_GROWTH for s in states[3:])
    # transition between sigma_max = 80 and 100 MPa
    assert points[2].predicted_cycles == math.inf
    assert math.isfinite(points[3].predicted_cycles)


def test_stress_range_sweep_life_is_finite_after_the_transition():
    active = [p for p in _stress_points() if p.state.is_growing]
    assert all(math.isfinite(p.predicted_cycles) for p in active)
    assert _decreasing([p.predicted_cycles for p in active])


def test_stress_range_sweep_moves_both_driving_force_and_boundary():
    points = _stress_points()
    assert all(
        later > earlier
        for earlier, later in zip(
            [p.initial_delta_k for p in points], [p.initial_delta_k for p in points][1:]
        )
    )
    assert _decreasing([p.critical_crack_length for p in points])


# --------------------------------------------------------------------------
# AD/AE/AF. fixed delta_sigma: threshold state frozen, only a_c moves
# --------------------------------------------------------------------------
def _fixed_range_points():
    return max_stress_sensitivity_with_threshold(
        tuple(v * 1.0e6 for v in (100, 120, 140, 160, 180)),
        100.0e6, A0, PANEL, LAW, TOU, T4,
    )


def test_fixed_range_sweep_leaves_delta_k_and_threshold_state_unchanged():
    points = _fixed_range_points()
    for point in points[1:]:
        assert point.initial_delta_k == pytest.approx(
            points[0].initial_delta_k, rel=1e-12
        )
        assert point.threshold_ratio == pytest.approx(
            points[0].threshold_ratio, rel=1e-12
        )
        assert point.threshold_crack_length == pytest.approx(
            points[0].threshold_crack_length, rel=1e-12
        )
    assert {p.state for p in points} == {GrowthState.ACTIVE_GROWTH}


def test_fixed_range_sweep_life_changes_only_through_the_fracture_boundary():
    points = _fixed_range_points()
    assert _decreasing([p.critical_crack_length for p in points])
    assert _decreasing([p.predicted_cycles for p in points])
    # ... while the driving range is identical throughout
    assert len({round(p.initial_delta_k, 6) for p in points}) == 1


# --------------------------------------------------------------------------
# AK. width affects the threshold state only through Y(a0)
# --------------------------------------------------------------------------
WIDTHS = (0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.5)


def _width_points():
    return width_sensitivity_with_threshold(
        WIDTHS, A0, CANONICAL_CYCLE, LAW, TOU, T4
    )


def test_width_raises_delta_k_at_the_initial_crack_through_geometry():
    points = _width_points()
    assert _decreasing([p.initial_geometry_factor for p in points])
    assert _decreasing([p.initial_delta_k for p in points])
    assert _decreasing([p.threshold_ratio for p in points])


def test_width_effect_on_threshold_state_is_small_here():
    """Measured, not assumed: at a0/W <= 0.025 the ratio moves only ~0.15 %,
    so every width in this sweep stays on the same side of the threshold."""
    points = _width_points()
    ratios = [p.threshold_ratio for p in points]
    assert {p.state for p in points} == {GrowthState.ACTIVE_GROWTH}
    assert (max(ratios) - min(ratios)) / min(ratios) < 2.0e-3


def test_width_sweep_life_and_boundary_still_move_with_width():
    points = _width_points()
    assert all(
        later > earlier
        for earlier, later in zip(
            [p.critical_crack_length for p in points],
            [p.critical_crack_length for p in points][1:],
        )
    )
    assert all(
        later > earlier
        for earlier, later in zip(
            [p.predicted_cycles for p in points],
            [p.predicted_cycles for p in points][1:],
        )
    )


def test_width_sweep_rejects_a_panel_too_narrow_for_the_flaw():
    with pytest.raises(ValueError, match="not admissible"):
        width_sensitivity_with_threshold(
            (1.0e-3,), A0, CANONICAL_CYCLE, LAW, TOU, T4
        )


# --------------------------------------------------------------------------
# AC/AL. toughness moves a_c but never the threshold state
# --------------------------------------------------------------------------
K_IC_SWEEP = tuple(mpa_sqrt_m_to_pa_sqrt_m(v) for v in (20, 25, 30, 35, 40, 50, 60))


def _toughness_points():
    return toughness_sensitivity_with_threshold(
        K_IC_SWEEP, A0, CANONICAL_CYCLE, PANEL, LAW, TOU, T4
    )


def test_threshold_classification_is_independent_of_toughness():
    """K_IC does not enter delta_K, so it cannot change the threshold verdict."""
    points = _toughness_points()
    assert {p.state for p in points} == {GrowthState.ACTIVE_GROWTH}
    for point in points[1:]:
        assert point.initial_delta_k == pytest.approx(
            points[0].initial_delta_k, rel=1e-12
        )
        assert point.threshold_ratio == pytest.approx(
            points[0].threshold_ratio, rel=1e-12
        )
        assert point.threshold_crack_length == pytest.approx(
            points[0].threshold_crack_length, rel=1e-12
        )


def test_toughness_still_moves_the_fracture_boundary_and_the_life():
    points = _toughness_points()
    assert all(
        later > earlier
        for earlier, later in zip(
            [p.critical_crack_length for p in points],
            [p.critical_crack_length for p in points][1:],
        )
    )
    assert all(
        later > earlier
        for earlier, later in zip(
            [p.predicted_cycles for p in points],
            [p.predicted_cycles for p in points][1:],
        )
    )


def test_sweep_points_expose_active_and_arrested_helpers():
    active = _initial_points()[-1]
    arrested = _initial_points()[0]
    assert active.active_growth and not active.arrested
    assert arrested.arrested and not arrested.active_growth
