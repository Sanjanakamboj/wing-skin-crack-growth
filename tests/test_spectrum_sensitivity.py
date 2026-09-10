"""AJ, AK, AM-AS. Spectrum sensitivity sweeps and the sequence-order study."""

from __future__ import annotations

import functools
import math

import pytest

from crackgrowth import (
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GROWTH_THRESHOLD,
    CANONICAL_PARIS_LAW,
    CANONICAL_SPECTRUM,
    CrackGrowthThreshold,
    SpectrumStatus,
    active_block_count,
    block_count_sensitivity,
    initial_flaw_sensitivity_under_spectrum,
    mpa_sqrt_m_to_pa_sqrt_m,
    scaled_spectrum,
    simulate_repeated_spectrum,
    spectrum_scale_sensitivity,
    spectrum_with_block_count,
    threshold_sensitivity_under_spectrum,
    toughness_sensitivity_under_spectrum,
    width_sensitivity_under_spectrum,
)

G = CANONICAL_FINITE_WIDTH_GEOMETRY
LAW = CANONICAL_PARIS_LAW
TOU = CANONICAL_FRACTURE_TOUGHNESS
TH = CANONICAL_GROWTH_THRESHOLD
A0 = 1.0e-3
S = CANONICAL_SPECTRUM

#: These sweeps assert TRENDS and STATUSES, not exact cycle counts, so they use
#: a looser block-solve tolerance than the default to keep the suite quick. The
#: exact canonical numbers are locked in test_spectrum_life.py at the default
#: tolerance.
SWEEP_TOL = 1.0e-9


def _threshold(mpa):
    return CrackGrowthThreshold(
        name=f"{mpa}",
        delta_k_threshold=mpa_sqrt_m_to_pa_sqrt_m(mpa),
        source_note="ILLUSTRATIVE (test fixture)",
    )


def _decreasing(values):
    return all(b < a for a, b in zip(values, values[1:]))


def _increasing(values):
    return all(b > a for a, b in zip(values, values[1:]))


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def test_scaled_spectrum_scales_stresses_but_not_counts_or_order():
    scaled = scaled_spectrum(S, 1.5)
    assert scaled.block_names == S.block_names
    assert [b.cycle_count for b in scaled] == [b.cycle_count for b in S]
    for original, new in zip(S, scaled):
        assert new.sigma_max == pytest.approx(original.sigma_max * 1.5)
        assert new.sigma_min == pytest.approx(original.sigma_min * 1.5)
        assert new.stress_ratio == pytest.approx(original.stress_ratio)


def test_spectrum_with_block_count_replaces_or_removes():
    fewer = spectrum_with_block_count(S, "C severe gust", 5)
    assert [b.cycle_count for b in fewer] == [10_000, 1_000, 5]
    removed = spectrum_with_block_count(S, "C severe gust", 0)
    assert removed.block_names == ("A low-amplitude", "B manoeuvre")


def test_spectrum_with_block_count_rejects_bad_input():
    with pytest.raises(ValueError, match="no block named"):
        spectrum_with_block_count(S, "nope", 5)
    with pytest.raises(ValueError, match="non-negative"):
        spectrum_with_block_count(S, "C severe gust", -1)
    with pytest.raises(TypeError):
        spectrum_with_block_count(S, "C severe gust", True)


def test_removing_the_last_block_is_rejected():
    single = spectrum_with_block_count(
        spectrum_with_block_count(S, "C severe gust", 0), "B manoeuvre", 0
    )
    with pytest.raises(ValueError, match="empty spectrum"):
        spectrum_with_block_count(single, "A low-amplitude", 0)


def test_active_block_count_at_the_initial_flaw():
    assert active_block_count(A0, S, G, TH) == 2  # low block is arrested
    assert active_block_count(A0, S, G, None) == 3
    assert active_block_count(4.0e-3, S, G, TH) == 3
    assert active_block_count(0.25e-3, S, G, TH) == 0


# --------------------------------------------------------------------------
# AM. stress scale
# --------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def _scale_points():
    return spectrum_scale_sensitivity(
        (0.75, 0.9, 1.0, 1.1, 1.25), A0, S, G, LAW, TOU, TH,
        tolerance=SWEEP_TOL,
    )


def test_higher_stress_scale_shortens_life_monotonically():
    points = _scale_points()
    assert all(p.status is SpectrumStatus.FRACTURE_REACHED for p in points)
    assert _decreasing([p.total_cycles for p in points])
    assert _decreasing([p.completed_full_spectra for p in points])


def test_low_stress_scale_leaves_more_blocks_arrested():
    points = _scale_points()
    assert points[0].parameter_value == pytest.approx(0.75)
    assert points[0].active_blocks_at_start == 1
    assert points[-1].active_blocks_at_start == 2
    assert points[0].blocks_activating_later == 2


def test_severe_block_governs_fracture_across_the_scale_sweep():
    assert {p.fracture_block_name for p in _scale_points()} == {"C severe gust"}


# --------------------------------------------------------------------------
# AN. severe-block count
# --------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def _severe_points():
    return block_count_sensitivity(
        "C severe gust", (0, 10, 50, 100, 200), A0, S, G, LAW, TOU, TH,
        tolerance=SWEEP_TOL,
    )


def test_more_severe_cycles_shorten_life():
    points = _severe_points()[1:]  # skip the removed-block row
    assert _decreasing([p.total_cycles for p in points])
    assert _decreasing([p.completed_full_spectra for p in points])


def test_removing_the_severe_block_moves_fracture_to_the_next_block():
    removed = _severe_points()[0]
    assert removed.parameter_value == 0.0
    assert removed.total_blocks == 2
    assert removed.fracture_block_name == "B manoeuvre"
    assert removed.minimum_block_critical_crack_length > (
        _severe_points()[1].minimum_block_critical_crack_length
    )


def test_severe_block_governs_fracture_whenever_present():
    assert {p.fracture_block_name for p in _severe_points()[1:]} == {"C severe gust"}


# --------------------------------------------------------------------------
# AO. low-block count -- the counter-intuitive case
# --------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def _low_points():
    return block_count_sensitivity(
        "A low-amplitude", (0, 1_000, 10_000, 50_000), A0, S, G, LAW, TOU, TH,
        tolerance=SWEEP_TOL,
    )


def test_more_arrested_low_cycles_raise_the_total_cycle_count():
    """Counter-intuitive but correct: while the low block is arrested its
    cycles add to the total without adding growth, so total cycles to fracture
    RISE. The spectra needed still fall, because once it activates it does
    contribute."""
    points = _low_points()
    assert _increasing([p.total_cycles for p in points])
    assert _decreasing([p.completed_full_spectra for p in points])


def test_removing_the_low_block_gives_the_shortest_cycle_count():
    points = _low_points()
    assert points[0].parameter_value == 0.0
    assert points[0].total_cycles < points[-1].total_cycles / 10.0
    assert points[0].active_blocks_at_start == 2
    assert points[0].total_blocks == 2


# --------------------------------------------------------------------------
# AP. threshold
# --------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def _threshold_points():
    return threshold_sensitivity_under_spectrum(
        (None, _threshold(2.0), _threshold(4.0), _threshold(6.0), _threshold(8.0)),
        A0, S, G, LAW, TOU, tolerance=SWEEP_TOL,
    )


def test_a_high_enough_threshold_arrests_the_whole_spectrum():
    last = _threshold_points()[-1]
    assert last.status is SpectrumStatus.SPECTRUM_ARRESTED
    assert last.total_cycles == math.inf
    assert last.active_blocks_at_start == 0


def test_a_low_threshold_matches_the_disabled_screen():
    """delta_K_th = 2 is below every block's delta_K(a0), so nothing is
    arrested and the answer equals the no-threshold run exactly."""
    disabled, low = _threshold_points()[0], _threshold_points()[1]
    assert disabled.active_blocks_at_start == 3
    assert low.active_blocks_at_start == 3
    assert low.total_cycles == pytest.approx(disabled.total_cycles, rel=1e-12)


def test_raising_the_threshold_lengthens_the_predicted_life():
    """Non-decreasing, not strictly increasing: any threshold below every
    block's delta_K(a0) arrests nothing and gives the identical no-threshold
    life, so the sweep starts with a flat step before it rises."""
    finite = [p for p in _threshold_points() if math.isfinite(p.total_cycles)]
    cycles = [p.total_cycles for p in finite]
    assert all(b >= a for a, b in zip(cycles, cycles[1:]))
    assert cycles[-1] > cycles[0]
    assert cycles[0] == pytest.approx(cycles[1], rel=1e-12)  # the flat step


def test_disabled_threshold_reports_no_later_activation():
    assert _threshold_points()[0].blocks_activating_later == 0


# --------------------------------------------------------------------------
# AQ. toughness
# --------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def _toughness_points():
    return toughness_sensitivity_under_spectrum(
        tuple(mpa_sqrt_m_to_pa_sqrt_m(v) for v in (20, 30, 40, 50)),
        A0, S, G, LAW, TOU, TH, tolerance=SWEEP_TOL,
    )


def test_toughness_does_not_change_threshold_activity():
    """K_IC does not enter delta_K, so it cannot change which blocks are
    active at a given crack length."""
    points = _toughness_points()
    assert {p.active_blocks_at_start for p in points} == {2}
    assert {p.blocks_activating_later for p in points} == {1}


def test_higher_toughness_lengthens_life_with_diminishing_returns():
    points = _toughness_points()
    assert _increasing([p.total_cycles for p in points])
    assert _increasing([p.minimum_block_critical_crack_length for p in points])
    first_gain = points[1].total_cycles - points[0].total_cycles
    last_gain = points[-1].total_cycles - points[-2].total_cycles
    assert last_gain < first_gain


# --------------------------------------------------------------------------
# AR. width
# --------------------------------------------------------------------------
def test_wider_panel_lengthens_life_and_converges():
    points = width_sensitivity_under_spectrum(
        (0.05, 0.075, 0.1, 0.2, 0.5), A0, S, LAW, TOU, TH, tolerance=SWEEP_TOL
    )
    assert _increasing([p.total_cycles for p in points])
    assert _increasing([p.minimum_block_critical_crack_length for p in points])
    # the widest panels differ only slightly -- convergence to the infinite plate
    assert points[-1].total_cycles / points[-2].total_cycles < 1.01


def test_width_sweep_rejects_a_panel_too_narrow_for_the_flaw():
    with pytest.raises(ValueError, match="not admissible"):
        width_sensitivity_under_spectrum((1.0e-3,), A0, S, LAW, TOU, TH)


# --------------------------------------------------------------------------
# AS. initial flaw
# --------------------------------------------------------------------------
@functools.lru_cache(maxsize=None)
def _flaw_points():
    return initial_flaw_sensitivity_under_spectrum(
        (0.25e-3, 0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3), S, G, LAW, TOU, TH,
        tolerance=SWEEP_TOL,
    )


def test_a_small_enough_flaw_arrests_the_whole_spectrum():
    smallest = _flaw_points()[0]
    assert smallest.status is SpectrumStatus.SPECTRUM_ARRESTED
    assert smallest.active_blocks_at_start == 0
    assert smallest.total_cycles == math.inf


def test_larger_flaws_activate_more_blocks_and_shorten_life():
    points = _flaw_points()
    active = [p for p in points if math.isfinite(p.total_cycles)]
    assert _decreasing([p.total_cycles for p in active])
    assert [p.active_blocks_at_start for p in points] == [0, 1, 2, 2, 3]


def test_largest_flaw_leaves_nothing_to_activate_later():
    assert _flaw_points()[-1].active_blocks_at_start == 3
    assert _flaw_points()[-1].blocks_activating_later == 0


# --------------------------------------------------------------------------
# AJ/AK. sequence-order study
# --------------------------------------------------------------------------
ORDERS = ((0, 1, 2), (2, 1, 0), (1, 0, 2), (2, 0, 1))


@functools.lru_cache(maxsize=None)
def _fixed_spectra_crack(order, repeats):
    """Crack length after a FIXED number of spectra, which removes the
    block-granularity quantisation present in a cycles-to-fracture comparison."""
    # The DEFAULT (tight) tolerance is used here deliberately: the sequence
    # effect being measured is finer than the bisection error at SWEEP_TOL,
    # which would swamp it.
    return simulate_repeated_spectrum(
        A0, S.reordered(order), G, LAW, TOU, TH, max_spectrum_repeats=repeats
    ).final_crack_length


def test_order_is_irrelevant_while_the_low_block_is_still_arrested():
    """Before the low block activates, only the two always-active blocks grow
    the crack, and their order does not matter: the growth per pass is
    order-independent to numerical precision."""
    cracks = [_fixed_spectra_crack(order, 200) for order in ORDERS]
    for value in cracks[1:]:
        assert value == pytest.approx(cracks[0], rel=1e-7)


def test_sequence_order_matters_once_the_low_block_activates():
    """A real, if small, sequence effect, and it appears ONLY through threshold
    activation: applying the low block FIRST each pass catches it at the
    smallest crack of that pass, where it is most often still arrested, so it
    contributes less growth.

    Nothing here is overload retardation or closure memory -- the only state
    carried between blocks is the crack length.
    """
    cracks = {order: _fixed_spectra_crack(order, 400) for order in ORDERS}
    low_first = cracks[(0, 1, 2)]
    others = [cracks[o] for o in ORDERS[1:]]
    assert all(v > low_first for v in others)
    relative = max(abs(v - low_first) / low_first for v in others)
    assert 1.0e-4 < relative < 1.0e-2
    assert relative == pytest.approx(3.58e-3, rel=5e-2)


def test_orders_that_do_not_start_with_the_low_block_nearly_coincide():
    values = [_fixed_spectra_crack(o, 400) for o in ORDERS[1:]]
    for value in values[1:]:
        assert value == pytest.approx(values[0], rel=1e-6)


def test_sequence_order_result_is_deterministic():
    for order in ORDERS:
        first = simulate_repeated_spectrum(
            A0, S.reordered(order), G, LAW, TOU, TH, max_spectrum_repeats=50
        )
        second = simulate_repeated_spectrum(
            A0, S.reordered(order), G, LAW, TOU, TH, max_spectrum_repeats=50
        )
        assert first.final_crack_length == second.final_crack_length


def test_reordering_does_not_change_the_cycles_per_spectrum():
    for order in ORDERS:
        assert S.reordered(order).cycles_per_spectrum == S.cycles_per_spectrum
