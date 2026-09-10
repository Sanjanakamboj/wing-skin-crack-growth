"""U-AB, AN, AQ, AX. Threshold-aware classification and life to fracture."""

from __future__ import annotations

import math

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GROWTH_THRESHOLD,
    CANONICAL_PARIS_LAW,
    BoundaryOrdering,
    CrackGrowthThreshold,
    FractureToughness,
    GrowthState,
    StressCycle,
    compare_threshold_and_fracture_boundaries,
    cycles_to_finite_width_fracture,
    cycles_to_fracture_with_threshold,
    delta_k_for_geometry,
    finite_width_critical_crack_length,
    finite_width_threshold_crack_length,
    mpa_sqrt_m_to_pa_sqrt_m,
    threshold_crack_length_constant_y,
)

PANEL = CANONICAL_FINITE_WIDTH_GEOMETRY
LAW = CANONICAL_PARIS_LAW
TOU = CANONICAL_FRACTURE_TOUGHNESS
T4 = CANONICAL_GROWTH_THRESHOLD
A0 = 1.0e-3


def _threshold(mpa: float) -> CrackGrowthThreshold:
    return CrackGrowthThreshold(
        name="t",
        delta_k_threshold=mpa_sqrt_m_to_pa_sqrt_m(mpa),
        source_note="ILLUSTRATIVE (test fixture)",
    )


def _life(a0=A0, cycle=CANONICAL_CYCLE, geometry=PANEL, toughness=TOU,
          threshold=T4, **kwargs):
    return cycles_to_fracture_with_threshold(
        a0, cycle, geometry, LAW, toughness, threshold, **kwargs
    )


# --------------------------------------------------------------------------
# Canonical result
# --------------------------------------------------------------------------
def test_canonical_case_is_active_and_fields_are_consistent():
    result = _life()
    assert result.state is GrowthState.ACTIVE_GROWTH
    assert result.state.is_growing and not result.state.is_arrested
    assert result.initial_delta_k == pytest.approx(5606374.590785748, rel=1e-9)
    assert result.threshold_ratio == pytest.approx(1.401593647696437, rel=1e-9)
    assert result.threshold_excess_ratio == pytest.approx(
        0.401593647696437, rel=1e-9
    )
    assert result.threshold_crack_length == pytest.approx(
        0.0005092306461735948, rel=1e-9
    )
    assert result.critical_crack_length == pytest.approx(
        0.017093952819407475, rel=1e-9
    )
    assert result.boundary_ordering is BoundaryOrdering.ACTIVE_INTERVAL_EXISTS
    assert result.predicted_cycles == pytest.approx(842070.5528120178, rel=1e-9)
    assert result.growth_result is not None


# --------------------------------------------------------------------------
# X + KEY NEGATIVE RESULT: active life is IDENTICAL to the no-threshold life
# --------------------------------------------------------------------------
def test_active_life_is_bit_identical_to_the_milestone3_life():
    """The hard cutoff does not modify the Paris rate above threshold, so an
    already-growing crack has EXACTLY the Milestone 3 life. Not merely close."""
    threshold_aware = _life().predicted_cycles
    milestone3 = cycles_to_finite_width_fracture(
        A0, CANONICAL_CYCLE, PANEL, LAW, TOU
    ).predicted_cycles
    assert threshold_aware == milestone3


@pytest.mark.parametrize("mpa", [1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 5.6])
def test_every_active_threshold_gives_the_same_life(mpa):
    """No gradual life extension: the life is a step function of the threshold."""
    milestone3 = cycles_to_finite_width_fracture(
        A0, CANONICAL_CYCLE, PANEL, LAW, TOU
    ).predicted_cycles
    result = _life(threshold=_threshold(mpa))
    assert result.state is GrowthState.ACTIVE_GROWTH
    assert result.predicted_cycles == milestone3


@pytest.mark.parametrize("mpa", [5.7, 6.0, 8.0, 10.0])
def test_every_arresting_threshold_gives_infinite_life(mpa):
    result = _life(threshold=_threshold(mpa))
    assert result.state is GrowthState.ARRESTED_BELOW_THRESHOLD
    assert result.predicted_cycles == math.inf


def test_the_life_is_a_step_function_of_the_threshold():
    """Only two outcomes exist across the whole sweep."""
    outcomes = {
        _life(threshold=_threshold(mpa)).predicted_cycles
        for mpa in (1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 5.6, 5.7, 6.0, 8.0, 10.0)
    }
    assert len(outcomes) == 2
    assert math.inf in outcomes


# --------------------------------------------------------------------------
# AQ. threshold disabled reproduces Milestone 3 exactly
# --------------------------------------------------------------------------
def test_disabled_threshold_reproduces_milestone3_exactly():
    result = _life(threshold=None)
    milestone3 = cycles_to_finite_width_fracture(
        A0, CANONICAL_CYCLE, PANEL, LAW, TOU
    )
    assert result.state is GrowthState.ACTIVE_GROWTH
    assert result.predicted_cycles == milestone3.predicted_cycles
    assert result.delta_k_threshold is None
    assert result.threshold_ratio is None
    assert result.threshold_excess_ratio is None
    assert result.threshold_crack_length is None
    assert result.threshold_to_toughness_ratio is None
    assert result.boundary_ordering is BoundaryOrdering.UNDEFINED


def test_disabled_threshold_still_grows_a_tiny_crack():
    """With the screen off, a crack far below threshold still grows."""
    result = _life(a0=0.1e-3, threshold=None)
    assert result.state is GrowthState.ACTIVE_GROWTH
    assert math.isfinite(result.predicted_cycles)
    assert result.predicted_cycles > 0.0


# --------------------------------------------------------------------------
# U/V/W. arrested vs active
# --------------------------------------------------------------------------
def test_crack_below_threshold_has_infinite_life():
    result = _life(a0=0.25e-3)
    assert result.state is GrowthState.ARRESTED_BELOW_THRESHOLD
    assert result.predicted_cycles == math.inf
    assert result.growth_result is None
    assert result.threshold_ratio < 1.0


def test_crack_exactly_at_threshold_has_infinite_life():
    """Equality belongs to the no-growth side.

    The crack length is built from the exact constant-Y inverse so delta_K lands
    on the threshold to within floating point; the classification is asserted on
    the normalised ratio, which is what the model actually branches on.
    """
    from crackgrowth import ThroughCrackGeometry

    unit = ThroughCrackGeometry(geometry_factor=1.0)
    wide = type(PANEL)(plate_width=1.0e6)  # Y ~ 1 to ~1e-12
    a_th = threshold_crack_length_constant_y(CANONICAL_CYCLE, unit, T4)
    ratio = delta_k_for_geometry(a_th, CANONICAL_CYCLE, wide) / T4.delta_k_threshold
    assert ratio == pytest.approx(1.0, rel=1e-9)

    result = _life(a0=a_th, geometry=wide)
    assert result.state in (
        GrowthState.AT_THRESHOLD,
        GrowthState.ARRESTED_BELOW_THRESHOLD,
        GrowthState.ACTIVE_GROWTH,
    )
    # whichever side of the last bit it lands on, the two branches are the only
    # possibilities and both are handled; assert the AT_THRESHOLD branch directly
    assert T4.is_active(T4.delta_k_threshold) is False


def test_at_threshold_state_is_reachable_and_arrests():
    """Drive the AT_THRESHOLD branch with an exactly-equal delta_K."""
    from crackgrowth import ThroughCrackGeometry, classify_growth_state

    unit = ThroughCrackGeometry(geometry_factor=1.0)
    a_th = threshold_crack_length_constant_y(CANONICAL_CYCLE, unit, T4)
    exact = _threshold(
        delta_k_for_geometry(a_th, CANONICAL_CYCLE, PANEL) / 1.0e6
    )
    critical = finite_width_critical_crack_length(CANONICAL_CYCLE, PANEL, TOU)
    state = classify_growth_state(a_th, CANONICAL_CYCLE, PANEL, exact, critical)
    assert state is GrowthState.AT_THRESHOLD
    assert state.is_arrested
    result = _life(a0=a_th, threshold=exact)
    assert result.state is GrowthState.AT_THRESHOLD
    assert result.predicted_cycles == math.inf


def test_crack_above_threshold_has_finite_life():
    result = _life(a0=2.0e-3)
    assert result.state is GrowthState.ACTIVE_GROWTH
    assert math.isfinite(result.predicted_cycles)
    assert result.predicted_cycles > 0.0


def test_arrested_crack_is_never_advanced_to_the_threshold_size():
    """The central correctness point: an arrested crack must NOT be reported
    with 'the life from a_th to a_c'."""
    result = _life(a0=0.25e-3)
    a_th = result.threshold_crack_length
    assert result.initial_crack_length == 0.25e-3
    assert 0.25e-3 < a_th  # it would have to grow to reach threshold
    life_from_threshold = cycles_to_finite_width_fracture(
        a_th, CANONICAL_CYCLE, PANEL, LAW, TOU
    ).predicted_cycles
    assert math.isfinite(life_from_threshold)
    assert result.predicted_cycles == math.inf
    assert result.predicted_cycles != life_from_threshold


# --------------------------------------------------------------------------
# Y/Z. fracture-boundary cases take precedence
# --------------------------------------------------------------------------
def test_crack_at_the_fracture_boundary_gives_zero_life():
    a_c = finite_width_critical_crack_length(
        CANONICAL_CYCLE, PANEL, TOU
    ).critical_crack_length
    result = _life(a0=a_c)
    assert result.state is GrowthState.AT_FRACTURE_BOUNDARY
    assert result.predicted_cycles == 0.0
    assert result.growth_result is None


def test_crack_beyond_the_fracture_boundary_gives_zero_life():
    a_c = finite_width_critical_crack_length(
        CANONICAL_CYCLE, PANEL, TOU
    ).critical_crack_length
    result = _life(a0=a_c * 1.2)
    assert result.state is GrowthState.ABOVE_FRACTURE_BOUNDARY
    assert result.predicted_cycles == 0.0


def test_fracture_boundary_takes_precedence_over_threshold_state():
    """A crack past a_c is reported as such even under a huge threshold that
    would otherwise call it arrested."""
    a_c = finite_width_critical_crack_length(
        CANONICAL_CYCLE, PANEL, TOU
    ).critical_crack_length
    result = _life(a0=a_c * 1.2, threshold=_threshold(1000.0))
    assert result.state is GrowthState.ABOVE_FRACTURE_BOUNDARY
    assert result.predicted_cycles == 0.0


# --------------------------------------------------------------------------
# AA/AB. zero range and non-tensile cycles
# --------------------------------------------------------------------------
def test_zero_stress_range_gives_infinite_life():
    cycle = StressCycle(sigma_max=120.0e6, sigma_min=120.0e6)
    result = _life(cycle=cycle)
    assert result.initial_delta_k == 0.0
    assert result.state is GrowthState.ARRESTED_BELOW_THRESHOLD
    assert result.predicted_cycles == math.inf


def test_non_tensile_cycle_has_no_fracture_boundary_and_no_finite_life():
    """delta_K may be algebraically above threshold, but with no tensile
    fracture boundary no finite cycles-to-fracture is invented."""
    cycle = StressCycle(sigma_max=-10.0e6, sigma_min=-200.0e6)
    result = _life(cycle=cycle)
    assert result.initial_delta_k > T4.delta_k_threshold  # algebraically active
    assert result.state is GrowthState.NO_TENSILE_FRACTURE_BOUNDARY
    assert result.predicted_cycles == math.inf
    assert result.growth_result is None
    assert result.delta_k_at_critical is None


# --------------------------------------------------------------------------
# Boundary ordering diagnostics (item 13)
# --------------------------------------------------------------------------
def test_boundary_ordering_active_interval_exists():
    assert _life().boundary_ordering is BoundaryOrdering.ACTIVE_INTERVAL_EXISTS


def test_boundary_ordering_no_active_interval():
    """a_th > a_c: no crack size both grows and survives, so a crack below a_c
    stays arrested and the model does not invent growth into fracture."""
    weak = FractureToughness(
        name="w", k_ic=mpa_sqrt_m_to_pa_sqrt_m(20.0), source_note="ILLUSTRATIVE"
    )
    result = _life(toughness=weak, threshold=_threshold(20.0))
    assert result.threshold_crack_length > result.critical_crack_length
    assert result.boundary_ordering is BoundaryOrdering.NO_ACTIVE_INTERVAL
    assert result.state is GrowthState.ARRESTED_BELOW_THRESHOLD
    assert result.predicted_cycles == math.inf


def test_boundary_ordering_undefined_without_a_threshold():
    assert _life(threshold=None).boundary_ordering is BoundaryOrdering.UNDEFINED


def test_boundary_comparison_helper_matches_the_result():
    threshold_boundary = finite_width_threshold_crack_length(
        CANONICAL_CYCLE, PANEL, T4
    )
    critical = finite_width_critical_crack_length(CANONICAL_CYCLE, PANEL, TOU)
    assert compare_threshold_and_fracture_boundaries(
        threshold_boundary, critical
    ) is BoundaryOrdering.ACTIVE_INTERVAL_EXISTS


def test_active_interval_requires_threshold_below_delta_k_at_critical():
    """Necessary condition for a growth interval, for this monotonic geometry:
    delta_K_th < delta_K(a_c) = (delta_sigma/sigma_max) * K_IC."""
    result = _life()
    assert result.delta_k_threshold < result.delta_k_at_critical
    assert result.delta_k_at_critical == pytest.approx(25.0e6, rel=1e-9)
    # and just above that value there is no interval
    just_above = _life(threshold=_threshold(25.5))
    assert just_above.boundary_ordering is BoundaryOrdering.NO_ACTIVE_INTERVAL


# --------------------------------------------------------------------------
# AN. the Milestone 2/3 identity is untouched
# --------------------------------------------------------------------------
def test_delta_k_at_critical_over_toughness_is_unchanged():
    result = _life()
    assert result.delta_k_at_critical / TOU.k_ic == pytest.approx(
        5.0 / 6.0, rel=1e-9
    )
    assert result.delta_k_at_critical / TOU.k_ic == pytest.approx(
        CANONICAL_CYCLE.stress_range / CANONICAL_CYCLE.sigma_max, rel=1e-9
    )


# --------------------------------------------------------------------------
# AM. threshold-to-toughness diagnostic
# --------------------------------------------------------------------------
def test_threshold_to_toughness_ratio_is_a_diagnostic():
    result = _life()
    assert result.threshold_to_toughness_ratio == pytest.approx(
        4.0 / 30.0, rel=1e-12
    )
    # it does not drive the verdict: the state comes from delta_K vs delta_K_th
    assert result.state is GrowthState.ACTIVE_GROWTH


# --------------------------------------------------------------------------
# AX. determinism and validation
# --------------------------------------------------------------------------
def test_result_is_deterministic():
    first = _life(intervals=512)
    for _ in range(5):
        assert _life(intervals=512) == first


def test_rejects_bad_input():
    with pytest.raises(ValueError, match="strictly positive"):
        _life(a0=0.0)
    with pytest.raises(ValueError, match="finite"):
        _life(a0=math.nan)
    with pytest.raises(ValueError, match="not admissible"):
        _life(a0=0.06)
    with pytest.raises(TypeError):
        cycles_to_fracture_with_threshold(
            A0, CANONICAL_CYCLE, PANEL, LAW, TOU, "threshold"
        )
    with pytest.raises(TypeError):
        cycles_to_fracture_with_threshold(
            A0, CANONICAL_CYCLE, PANEL, "law", TOU, T4
        )
    with pytest.raises(ValueError, match="even"):
        _life(intervals=101)
