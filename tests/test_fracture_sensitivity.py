"""AC-AH. Fracture sensitivity sweeps and the residual-strength table."""

from __future__ import annotations

import math

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GEOMETRY,
    CANONICAL_PARIS_LAW,
    ILLUSTRATIVE_TOUGHNESS_DISCLAIMER,
    ThroughCrackGeometry,
    critical_crack_length,
    geometry_factor_sensitivity,
    initial_crack_sensitivity_to_fracture,
    max_stress_sensitivity_at_fixed_range,
    mpa_sqrt_m_to_pa_sqrt_m,
    residual_strength_table,
    stress_range_sensitivity_at_fixed_min,
    toughness_sensitivity,
)

UNIT = CANONICAL_GEOMETRY
LAW = CANONICAL_PARIS_LAW
TOUGH = CANONICAL_FRACTURE_TOUGHNESS
A0 = 1.0e-3

K_IC_SWEEP = tuple(mpa_sqrt_m_to_pa_sqrt_m(v) for v in (20, 25, 30, 35, 40, 50, 60))
SIGMA_MAX_SWEEP = (100.0e6, 120.0e6, 140.0e6, 160.0e6, 180.0e6)
Y_SWEEP = (0.8, 0.9, 1.0, 1.1, 1.2)
A0_SWEEP = (0.25e-3, 0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3)


def _monotonic_decreasing(values):
    return all(later < earlier for earlier, later in zip(values, values[1:]))


def _monotonic_increasing(values):
    return all(later > earlier for earlier, later in zip(values, values[1:]))


# --------------------------------------------------------------------------
# AC. toughness sensitivity
# --------------------------------------------------------------------------
def _toughness_points(intervals=1000):
    return toughness_sensitivity(
        K_IC_SWEEP, A0, CANONICAL_CYCLE, UNIT, LAW, TOUGH, intervals=intervals
    )


def test_toughness_sweep_increases_critical_size_and_life():
    points = _toughness_points()
    assert _monotonic_increasing([p.critical_crack_length for p in points])
    assert _monotonic_increasing([p.numerical_cycles for p in points])


def test_toughness_sweep_critical_size_is_exactly_quadratic():
    points = _toughness_points()
    reference = points[0]
    for point in points:
        ratio = point.k_ic / reference.k_ic
        assert point.critical_crack_length == pytest.approx(
            reference.critical_crack_length * ratio**2, rel=1e-12
        )


def test_toughness_sweep_leaves_growth_driving_force_untouched():
    """K_IC moves the endpoint only; delta_K and da/dN at a0 are unchanged."""
    points = _toughness_points()
    for point in points[1:]:
        assert point.initial_delta_k == pytest.approx(points[0].initial_delta_k)
        assert point.initial_growth_rate == pytest.approx(
            points[0].initial_growth_rate
        )


def test_toughness_sweep_initial_margin_rises_with_toughness():
    assert _monotonic_increasing([p.initial_margin for p in _toughness_points()])
    assert _monotonic_decreasing(
        [p.initial_utilization for p in _toughness_points()]
    )


def test_toughness_sweep_preserves_illustrative_provenance():
    """A sweep of illustrative values must stay flagged as illustrative."""
    assert ILLUSTRATIVE_TOUGHNESS_DISCLAIMER in TOUGH.source_note
    points = _toughness_points()
    assert len(points) == len(K_IC_SWEEP)


def test_toughness_sweep_agrees_with_the_analytical_reference():
    # The K_IC = 60 MPa*sqrt(m) case puts a_c at 79.6 mm (a_c/a0 ~ 80), so the
    # grid is refined for this check; see the Milestone 2 span note.
    for point in _toughness_points(intervals=8000):
        assert point.numerical_cycles == pytest.approx(
            point.analytical_cycles, rel=1e-8
        )


# --------------------------------------------------------------------------
# AD/AE. sigma_max at fixed delta_sigma -- the K_max vs delta_K distinction
# --------------------------------------------------------------------------
def _fixed_range_points():
    return max_stress_sensitivity_at_fixed_range(
        SIGMA_MAX_SWEEP, 100.0e6, A0, UNIT, LAW, TOUGH
    )


def test_fixed_range_sweep_holds_the_stress_range_constant():
    for point in _fixed_range_points():
        assert point.stress_range == pytest.approx(100.0e6)


def test_fixed_range_sweep_leaves_growth_rate_unchanged():
    points = _fixed_range_points()
    for point in points[1:]:
        assert point.initial_delta_k == pytest.approx(
            points[0].initial_delta_k, rel=1e-12
        )
        assert point.initial_growth_rate == pytest.approx(
            points[0].initial_growth_rate, rel=1e-12
        )


def test_fixed_range_sweep_shrinks_the_critical_size_and_the_life():
    points = _fixed_range_points()
    assert _monotonic_decreasing([p.critical_crack_length for p in points])
    assert _monotonic_decreasing([p.numerical_cycles for p in points])


def test_fixed_range_sweep_critical_size_is_exactly_inverse_square():
    points = _fixed_range_points()
    reference = points[0]
    for point in points:
        ratio = point.sigma_max / reference.sigma_max
        assert point.critical_crack_length == pytest.approx(
            reference.critical_crack_length * ratio**-2, rel=1e-12
        )


def test_fixed_range_sweep_life_falls_only_because_the_endpoint_moves():
    """Identical growth rates everywhere, so the entire life reduction comes
    from the fracture boundary moving inward."""
    points = _fixed_range_points()
    assert points[-1].numerical_cycles < points[0].numerical_cycles
    assert points[-1].initial_growth_rate == pytest.approx(
        points[0].initial_growth_rate, rel=1e-12
    )


def test_fixed_range_sweep_rejects_non_positive_range():
    with pytest.raises(ValueError, match="strictly positive"):
        max_stress_sensitivity_at_fixed_range(
            SIGMA_MAX_SWEEP, 0.0, A0, UNIT, LAW, TOUGH
        )


# --------------------------------------------------------------------------
# Stress range at fixed sigma_min -- deliberately distinct from the above
# --------------------------------------------------------------------------
def _fixed_min_points():
    return stress_range_sensitivity_at_fixed_min(
        (50.0e6, 75.0e6, 100.0e6, 125.0e6, 150.0e6), 20.0e6, A0, UNIT, LAW, TOUGH
    )


def test_fixed_min_sweep_changes_both_driving_forces():
    points = _fixed_min_points()
    # growth driving force rises ...
    assert _monotonic_increasing([p.initial_delta_k for p in points])
    assert _monotonic_increasing([p.initial_growth_rate for p in points])
    # ... and the fracture boundary moves inward at the same time
    assert _monotonic_decreasing([p.critical_crack_length for p in points])
    assert _monotonic_decreasing([p.numerical_cycles for p in points])


def test_fixed_min_sweep_differs_from_the_fixed_range_sweep():
    """The two sweeps must not be conflated: one holds delta_sigma, the other
    holds sigma_min, and only the latter changes the growth rate."""
    fixed_range = _fixed_range_points()
    fixed_min = _fixed_min_points()
    range_rates = {round(p.initial_growth_rate, 24) for p in fixed_range}
    min_rates = {round(p.initial_growth_rate, 24) for p in fixed_min}
    assert len(range_rates) == 1
    assert len(min_rates) == len(fixed_min)


def test_fixed_min_sweep_sigma_max_tracks_the_range():
    for point in _fixed_min_points():
        assert point.sigma_max == pytest.approx(20.0e6 + point.parameter_value)
        assert point.stress_range == pytest.approx(point.parameter_value)


# --------------------------------------------------------------------------
# AF. geometry-factor sensitivity
# --------------------------------------------------------------------------
def _geometry_points():
    return geometry_factor_sensitivity(Y_SWEEP, A0, CANONICAL_CYCLE, LAW, TOUGH)


def test_geometry_sweep_shrinks_critical_size_and_life():
    points = _geometry_points()
    assert _monotonic_decreasing([p.critical_crack_length for p in points])
    assert _monotonic_decreasing([p.numerical_cycles for p in points])


def test_geometry_sweep_critical_size_is_exactly_inverse_square():
    points = _geometry_points()
    reference = next(p for p in points if p.parameter_value == 1.0)
    for point in points:
        assert point.critical_crack_length == pytest.approx(
            reference.critical_crack_length * point.parameter_value**-2, rel=1e-12
        )


def test_geometry_sweep_raises_the_growth_rate_too():
    """Y hurts twice: faster growth AND a nearer fracture boundary."""
    points = _geometry_points()
    assert _monotonic_increasing([p.initial_growth_rate for p in points])
    reference = next(p for p in points if p.parameter_value == 1.0)
    for point in points:
        assert point.initial_growth_rate == pytest.approx(
            reference.initial_growth_rate * point.parameter_value**LAW.m, rel=1e-12
        )


def test_geometry_sweep_life_falls_faster_than_either_effect_alone():
    points = _geometry_points()
    low = next(p for p in points if p.parameter_value == 0.8)
    high = next(p for p in points if p.parameter_value == 1.2)
    assert high.numerical_cycles < low.numerical_cycles
    # the combined effect exceeds the growth-rate effect (Y**-m) alone
    growth_only = low.numerical_cycles * (1.2 / 0.8) ** -LAW.m
    assert high.numerical_cycles < growth_only


# --------------------------------------------------------------------------
# AG. initial-crack sensitivity
# --------------------------------------------------------------------------
def _initial_crack_points(intervals=1000):
    return initial_crack_sensitivity_to_fracture(
        A0_SWEEP, CANONICAL_CYCLE, UNIT, LAW, TOUGH, intervals=intervals
    )


def test_initial_crack_sweep_life_decreases_monotonically():
    assert _monotonic_decreasing(
        [p.numerical_cycles for p in _initial_crack_points()]
    )


def test_initial_crack_sweep_shares_one_fracture_boundary():
    points = _initial_crack_points()
    a_c = critical_crack_length(
        CANONICAL_CYCLE, UNIT, TOUGH
    ).critical_crack_length
    for point in points:
        assert point.critical_crack_length == pytest.approx(a_c, rel=1e-12)


def test_initial_crack_sweep_fraction_and_utilization_rise():
    points = _initial_crack_points()
    assert _monotonic_increasing([p.initial_crack_fraction for p in points])
    assert _monotonic_increasing([p.initial_utilization for p in points])
    for point in points:
        # utilization = sqrt(a0 / a_c) for constant Y
        assert point.initial_utilization == pytest.approx(
            math.sqrt(point.initial_crack_fraction), rel=1e-12
        )
        assert point.initial_crack_fraction < 1.0


def test_initial_crack_sweep_agrees_with_the_analytical_reference():
    # a0 = 0.25 mm against a_c = 19.894 mm is a span of ~80, so refine the grid.
    for point in _initial_crack_points(intervals=8000):
        assert point.numerical_cycles == pytest.approx(
            point.analytical_cycles, rel=1e-8
        )


# --------------------------------------------------------------------------
# AH. residual-strength table
# --------------------------------------------------------------------------
TABLE_LENGTHS = (0.5e-3, 1.0e-3, 2.0e-3, 5.0e-3, 10.0e-3, 20.0e-3)


def _table():
    return residual_strength_table(TABLE_LENGTHS, CANONICAL_CYCLE, UNIT, TOUGH)


def test_residual_strength_table_is_monotonic():
    rows = _table()
    assert _monotonic_increasing([r.k_max for r in rows])
    assert _monotonic_decreasing([r.residual_strength for r in rows])
    assert _monotonic_increasing([r.utilization for r in rows])
    assert _monotonic_decreasing([r.margin for r in rows])


def test_residual_strength_table_passes_below_critical_and_fails_above():
    rows = _table()
    a_c = critical_crack_length(
        CANONICAL_CYCLE, UNIT, TOUGH
    ).critical_crack_length
    for row in rows:
        assert row.passes == (row.crack_length <= a_c)
    # the canonical boundary is 19.894 mm, so only the 20 mm row fails
    assert [r.passes for r in rows] == [True, True, True, True, True, False]


def test_residual_strength_table_rows_are_internally_consistent():
    for row in _table():
        assert row.utilization == pytest.approx(row.k_max / row.k_ic, rel=1e-12)
        assert row.margin == pytest.approx(row.k_ic / row.k_max - 1.0, rel=1e-12)
        assert row.utilization == pytest.approx(
            row.applied_sigma_max / row.residual_strength, rel=1e-12
        )
        assert row.applied_sigma_max == 120.0e6


def test_residual_strength_table_hand_values():
    rows = _table()
    by_length = {round(r.crack_length, 9): r for r in rows}
    # a = 1 mm: sigma_res = 30e6 / sqrt(pi*0.001) = 5.3523723484583133e8 Pa
    assert by_length[0.001].residual_strength == pytest.approx(
        5.3523723484583133e8, rel=1e-12
    )
    # a = 20 mm > a_c = 19.894 mm, so the screen fails
    # a = 20 mm: sigma_res = 30e6 / sqrt(pi*0.02) = 1.1968268412042981e8 Pa
    assert by_length[0.02].residual_strength == pytest.approx(
        1.1968268412042981e8, rel=1e-12
    )
    assert by_length[0.02].residual_strength < 120.0e6
    assert not by_length[0.02].passes


def test_sweeps_reject_wrong_geometry_type():
    with pytest.raises(TypeError):
        toughness_sensitivity(K_IC_SWEEP, A0, CANONICAL_CYCLE, "geo", LAW, TOUGH)


def test_geometry_sweep_accepts_plain_factors():
    points = geometry_factor_sensitivity((1.0,), A0, CANONICAL_CYCLE, LAW, TOUGH)
    assert points[0].geometry_factor == 1.0
    assert points[0].critical_crack_length == pytest.approx(
        critical_crack_length(
            CANONICAL_CYCLE, ThroughCrackGeometry(1.0), TOUGH
        ).critical_crack_length
    )
