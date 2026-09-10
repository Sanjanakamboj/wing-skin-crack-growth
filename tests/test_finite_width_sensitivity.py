"""AG-AM. Finite-width sensitivity sweeps."""

from __future__ import annotations

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_PARIS_LAW,
    FiniteWidthCenterCrack,
    ThroughCrackGeometry,
    finite_width_max_stress_sensitivity,
    finite_width_stress_range_sensitivity,
    geometry_amplification_table,
    initial_crack_sensitivity_finite_width,
    mpa_sqrt_m_to_pa_sqrt_m,
    stress_intensity,
    toughness_sensitivity_finite_width,
    width_sensitivity,
)

PANEL = CANONICAL_FINITE_WIDTH_GEOMETRY
LAW = CANONICAL_PARIS_LAW
K30 = CANONICAL_FRACTURE_TOUGHNESS
A0 = 1.0e-3

WIDTHS = (0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.5)
A0_SWEEP = (0.25e-3, 0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3, 8.0e-3)
SIGMA_MAX_SWEEP = (100.0e6, 120.0e6, 140.0e6, 160.0e6, 180.0e6)
K_IC_SWEEP = tuple(mpa_sqrt_m_to_pa_sqrt_m(v) for v in (20, 25, 30, 35, 40, 50, 60))


def _increasing(values):
    return all(later > earlier for earlier, later in zip(values, values[1:]))


def _decreasing(values):
    return all(later < earlier for earlier, later in zip(values, values[1:]))


# --------------------------------------------------------------------------
# AG. width sensitivity
# --------------------------------------------------------------------------
def _width_points():
    return width_sensitivity(WIDTHS, A0, CANONICAL_CYCLE, LAW, K30)


def test_wider_panel_gives_lower_y_larger_critical_size_and_longer_life():
    points = _width_points()
    assert _decreasing([p.initial_geometry_factor for p in points])
    assert _increasing([p.critical_crack_length for p in points])
    assert _increasing([p.predicted_cycles for p in points])


def test_width_sweep_ligament_fraction_rises_with_width():
    assert _increasing([p.ligament_fraction_at_critical for p in _width_points()])
    for point in _width_points():
        assert 0.0 < point.ligament_fraction_at_critical < 1.0


def test_width_sweep_converges_toward_the_infinite_plate():
    points = _width_points()
    widest = points[-1]
    assert widest.initial_geometry_factor == pytest.approx(1.0, abs=1e-4)
    assert widest.geometry_factor_at_critical == pytest.approx(1.0, abs=1e-2)


def test_width_sweep_rejects_a_width_too_narrow_for_the_flaw():
    with pytest.raises(ValueError, match="not admissible"):
        width_sensitivity((1.0e-3,), A0, CANONICAL_CYCLE, LAW, K30)


def test_width_sweep_records_the_ratio_of_flaw_to_width():
    for point in _width_points():
        assert point.initial_crack_width_ratio == pytest.approx(
            A0 / point.plate_width
        )


# --------------------------------------------------------------------------
# AH. initial-crack sensitivity
# --------------------------------------------------------------------------
def _initial_points():
    return initial_crack_sensitivity_finite_width(
        A0_SWEEP, CANONICAL_CYCLE, PANEL, LAW, K30
    )


def test_life_decreases_as_the_initial_flaw_grows():
    assert _decreasing([p.predicted_cycles for p in _initial_points()])


def test_initial_flaw_sweep_utilization_and_y_rise():
    points = _initial_points()
    assert _increasing([p.initial_utilization for p in points])
    assert _increasing([p.initial_geometry_factor for p in points])
    for point in points:
        assert point.initial_utilization < 1.0
        assert point.critical_crack_length == pytest.approx(
            points[0].critical_crack_length, rel=1e-9
        )


# --------------------------------------------------------------------------
# AI/AJ. sigma_max at fixed delta_sigma
# --------------------------------------------------------------------------
def _fixed_range_points():
    return finite_width_max_stress_sensitivity(
        SIGMA_MAX_SWEEP, 100.0e6, A0, PANEL, LAW, K30
    )


def test_fixed_range_sweep_leaves_delta_k_and_growth_rate_unchanged():
    points = _fixed_range_points()
    for point in points[1:]:
        assert point.stress_range == pytest.approx(100.0e6)
        assert point.initial_delta_k == pytest.approx(
            points[0].initial_delta_k, rel=1e-12
        )
        assert point.initial_growth_rate == pytest.approx(
            points[0].initial_growth_rate, rel=1e-12
        )


def test_fixed_range_sweep_shrinks_critical_size_and_life():
    points = _fixed_range_points()
    assert _decreasing([p.critical_crack_length for p in points])
    assert _decreasing([p.predicted_cycles for p in points])


def test_fixed_range_sweep_life_falls_only_because_the_endpoint_moves():
    points = _fixed_range_points()
    assert points[-1].predicted_cycles < points[0].predicted_cycles
    assert points[-1].initial_growth_rate == pytest.approx(
        points[0].initial_growth_rate, rel=1e-12
    )


# --------------------------------------------------------------------------
# AK. stress range at fixed sigma_min -- deliberately distinct
# --------------------------------------------------------------------------
def _fixed_min_points():
    return finite_width_stress_range_sensitivity(
        (80.0e6, 100.0e6, 120.0e6, 140.0e6, 160.0e6, 180.0e6),
        20.0e6,
        A0,
        PANEL,
        LAW,
        K30,
    )


def test_fixed_min_sweep_moves_both_driving_forces():
    points = _fixed_min_points()
    assert _increasing([p.stress_range for p in points])
    assert _increasing([p.initial_delta_k for p in points])
    assert _increasing([p.initial_growth_rate for p in points])
    assert _decreasing([p.critical_crack_length for p in points])
    assert _decreasing([p.predicted_cycles for p in points])


def test_the_two_stress_sweeps_are_not_the_same_experiment():
    fixed_range_rates = {
        round(p.initial_growth_rate, 24) for p in _fixed_range_points()
    }
    fixed_min_rates = {round(p.initial_growth_rate, 24) for p in _fixed_min_points()}
    assert len(fixed_range_rates) == 1
    assert len(fixed_min_rates) == len(_fixed_min_points())


# --------------------------------------------------------------------------
# AL/AM. toughness sensitivity, and the broken K_IC**2 law
# --------------------------------------------------------------------------
def _toughness_points():
    return toughness_sensitivity_finite_width(
        K_IC_SWEEP, A0, CANONICAL_CYCLE, PANEL, LAW, K30
    )


def test_life_and_critical_size_increase_with_toughness():
    points = _toughness_points()
    assert _increasing([p.critical_crack_length for p in points])
    assert _increasing([p.predicted_cycles for p in points])


def test_finite_width_critical_size_is_not_proportional_to_toughness_squared():
    """A key Milestone 3 finding: the exact Milestone 2 scaling a_c ~ K_IC**2
    does NOT survive finite width, because Y(a_c) itself grows with a_c."""
    points = _toughness_points()
    ratios = [p.critical_crack_length / p.k_ic**2 for p in points]
    # Normalise before comparing: the raw ratios are ~1e-17, where the default
    # absolute tolerance of pytest.approx would swamp any relative comparison.
    normalised = [r / ratios[0] for r in ratios]
    assert not all(
        n == pytest.approx(1.0, rel=1e-3) for n in normalised
    ), normalised
    # the ratio falls as toughness rises: bigger cracks are penalised more
    assert _decreasing(ratios)
    assert normalised[-1] < 0.75
    # for reference, Milestone 2 gives exactly 1.0 for every entry
    assert normalised[-1] == pytest.approx(0.4606, rel=1e-3)


def test_infinite_plate_critical_size_is_still_proportional_to_toughness_squared():
    """Contrast: the same normalised check on the Milestone 2 model is flat."""
    from crackgrowth import CANONICAL_GEOMETRY, critical_crack_length, FractureToughness

    ratios = []
    for k_ic in K_IC_SWEEP:
        toughness = FractureToughness(
            name="t", k_ic=k_ic, source_note="ILLUSTRATIVE (test)"
        )
        a_c = critical_crack_length(
            CANONICAL_CYCLE, CANONICAL_GEOMETRY, toughness
        ).critical_crack_length
        ratios.append(a_c / k_ic**2)
    normalised = [r / ratios[0] for r in ratios]
    assert all(n == pytest.approx(1.0, rel=1e-12) for n in normalised)


def test_toughness_gain_is_smaller_than_the_infinite_plate_would_predict():
    points = _toughness_points()
    low, high = points[0], points[-1]  # 20 -> 60 MPa*sqrt(m)
    infinite_plate_expectation = (high.k_ic / low.k_ic) ** 2
    actual = high.critical_crack_length / low.critical_crack_length
    assert actual < infinite_plate_expectation
    assert infinite_plate_expectation == pytest.approx(9.0, rel=1e-9)


def test_toughness_sweep_leaves_the_growth_driving_force_untouched():
    points = _toughness_points()
    for point in points[1:]:
        assert point.initial_delta_k == pytest.approx(points[0].initial_delta_k)
        assert point.initial_growth_rate == pytest.approx(
            points[0].initial_growth_rate
        )


def test_toughness_sweep_ligament_fraction_falls_as_toughness_rises():
    assert _decreasing(
        [p.ligament_fraction_at_critical for p in _toughness_points()]
    )


# --------------------------------------------------------------------------
# Geometry amplification table
# --------------------------------------------------------------------------
TABLE_LENGTHS = (0.5e-3, 1.0e-3, 2.0e-3, 5.0e-3, 10.0e-3, 15.0e-3, 20.0e-3, 30.0e-3, 40.0e-3)


def test_amplification_table_is_monotonic_and_consistent():
    rows = geometry_amplification_table(TABLE_LENGTHS, CANONICAL_CYCLE, PANEL)
    assert _increasing([r.geometry_factor for r in rows])
    assert _increasing([r.crack_width_ratio for r in rows])
    assert _increasing([r.k_max for r in rows])
    assert _decreasing([r.ligament_fraction for r in rows])
    for row in rows:
        assert row.amplification == pytest.approx(row.geometry_factor, rel=1e-12)
        assert row.total_crack_length == pytest.approx(2.0 * row.crack_length)
        assert row.k_max > row.k_max_infinite_plate
        assert row.delta_k == pytest.approx(row.k_max * 100.0 / 120.0, rel=1e-12)


def test_amplification_table_matches_the_infinite_plate_helper():
    unit = ThroughCrackGeometry(geometry_factor=1.0)
    for row in geometry_amplification_table(TABLE_LENGTHS, CANONICAL_CYCLE, PANEL):
        assert row.k_max_infinite_plate == pytest.approx(
            stress_intensity(120.0e6, row.crack_length, unit), rel=1e-12
        )


def test_amplification_table_rejects_inadmissible_crack_lengths():
    with pytest.raises(ValueError, match="not admissible"):
        geometry_amplification_table((10.0e-3, 60.0e-3), CANONICAL_CYCLE, PANEL)


def test_amplification_is_modest_for_small_cracks_and_large_for_big_ones():
    rows = {
        round(r.crack_length, 9): r
        for r in geometry_amplification_table(TABLE_LENGTHS, CANONICAL_CYCLE, PANEL)
    }
    assert rows[0.001].amplification == pytest.approx(1.00025, rel=1e-4)
    assert rows[0.04].amplification == pytest.approx(1.79891, rel=1e-4)


def test_narrower_panel_amplifies_more_at_the_same_crack_length():
    a = (10.0e-3,)
    narrow = geometry_amplification_table(
        a, CANONICAL_CYCLE, FiniteWidthCenterCrack(0.05)
    )[0]
    wide = geometry_amplification_table(
        a, CANONICAL_CYCLE, FiniteWidthCenterCrack(0.5)
    )[0]
    assert narrow.amplification > wide.amplification
    assert narrow.k_max > wide.k_max
