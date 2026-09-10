"""S-AL, AT-AU, BE. Repeated-spectrum simulation."""

from __future__ import annotations

import functools
import math

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GROWTH_THRESHOLD,
    CANONICAL_PARIS_LAW,
    CANONICAL_SPECTRUM,
    CrackGrowthThreshold,
    LoadSpectrum,
    SpectrumBlock,
    SpectrumStatus,
    StressCycle,
    block_boundary_table,
    cycles_to_finite_width_fracture,
    finite_width_critical_crack_length,
    k_max_for_geometry,
    minimum_block_critical_crack_length,
    mpa_sqrt_m_to_pa_sqrt_m,
    simulate_repeated_spectrum,
)

G = CANONICAL_FINITE_WIDTH_GEOMETRY
LAW = CANONICAL_PARIS_LAW
TOU = CANONICAL_FRACTURE_TOUGHNESS
TH = CANONICAL_GROWTH_THRESHOLD


def _sim(a0=1.0e-3, spectrum=CANONICAL_SPECTRUM, threshold=TH, toughness=TOU, **kw):
    return simulate_repeated_spectrum(
        a0, spectrum, G, LAW, toughness, threshold, **kw
    )


@functools.lru_cache(maxsize=None)
def _canonical():
    """The full canonical run, computed once and shared across tests.

    It takes a couple of seconds; the result is immutable, so caching it is
    safe and keeps the suite quick.
    """
    return _sim()


def _one_block(cycle, n, name="CA"):
    return LoadSpectrum(
        blocks=(SpectrumBlock(name=name, stress_cycle=cycle, cycle_count=n),),
        name="single-block",
        source_note="ILLUSTRATIVE (test fixture)",
    )


# --------------------------------------------------------------------------
# AG. one-block spectrum reproduces the Milestone 4 constant-amplitude life
# --------------------------------------------------------------------------
M4_LIFE = cycles_to_finite_width_fracture(
    1.0e-3, CANONICAL_CYCLE, G, LAW, TOU
).predicted_cycles


@pytest.mark.parametrize("n", [10_000, 50_000])
def test_single_block_spectrum_reproduces_constant_amplitude_life(n):
    """The strongest Milestone 5 regression: a spectrum of one repeated block
    must reproduce the verified constant-amplitude answer."""
    result = _sim(spectrum=_one_block(CANONICAL_CYCLE, n))
    assert result.status is SpectrumStatus.FRACTURE_REACHED
    assert result.total_cycles == pytest.approx(M4_LIFE, rel=1e-7)
    assert result.total_cycles == pytest.approx(842070.5528, rel=1e-7)
    assert result.final_fracture_utilization == pytest.approx(1.0, rel=1e-9)


def test_single_block_regression_improves_with_larger_blocks():
    coarse = _sim(spectrum=_one_block(CANONICAL_CYCLE, 10_000)).total_cycles
    fine = _sim(spectrum=_one_block(CANONICAL_CYCLE, 50_000)).total_cycles
    assert abs(fine - M4_LIFE) <= abs(coarse - M4_LIFE)


# --------------------------------------------------------------------------
# V/W/X. sequencing and cycle accounting
# --------------------------------------------------------------------------
def test_blocks_execute_in_spectrum_order():
    result = _sim(max_spectrum_repeats=3)
    names = [h.advance.block_name for h in result.history]
    assert names == list(CANONICAL_SPECTRUM.block_names) * 3
    repeats = [h.spectrum_repeat for h in result.history]
    assert repeats == [1, 1, 1, 2, 2, 2, 3, 3, 3]
    assert [h.block_index for h in result.history] == [0, 1, 2] * 3


def test_total_cycles_equal_the_sum_of_completed_block_cycles():
    result = _sim(max_spectrum_repeats=5)
    assert result.total_cycles == pytest.approx(
        sum(h.advance.completed_cycles for h in result.history)
    )
    assert result.total_cycles == pytest.approx(
        5 * CANONICAL_SPECTRUM.cycles_per_spectrum
    )


def test_cycles_before_each_block_are_consistent():
    result = _sim(max_spectrum_repeats=4)
    running = 0.0
    for execution in result.history:
        assert execution.cycles_before == pytest.approx(running)
        running += execution.advance.completed_cycles


def test_crack_length_is_carried_between_blocks():
    result = _sim(max_spectrum_repeats=3)
    for previous, current in zip(result.history, result.history[1:]):
        assert current.advance.start_crack_length == (
            previous.advance.end_crack_length
        )


# --------------------------------------------------------------------------
# Canonical result
# --------------------------------------------------------------------------
def test_canonical_spectrum_result():
    result = _canonical()
    assert result.status is SpectrumStatus.FRACTURE_REACHED
    assert result.completed_full_spectra == 565
    assert result.total_cycles == pytest.approx(6_282_500.0, rel=1e-9)
    assert result.fracture_block_name == "C severe gust"
    assert result.fracture_block_index == 2
    assert result.fracture_spectrum_repeat == 566
    assert result.final_crack_length == pytest.approx(12.0031e-3, rel=1e-4)
    assert result.cycles_per_spectrum == 11_100
    assert result.partial_spectrum_cycles == pytest.approx(11_000.0)


def test_canonical_total_cycles_decompose_exactly():
    result = _canonical()
    assert result.total_cycles == pytest.approx(
        result.completed_full_spectra * result.cycles_per_spectrum
        + result.partial_spectrum_cycles
    )


def test_canonical_fracture_is_detected_at_the_start_of_the_severe_block():
    """The crack passes the severe block's boundary during the milder blocks,
    so the severe block fractures on its very first cycle."""
    result = _canonical()
    assert result.fracture_cycle_within_block == 0.0
    assert result.final_crack_length > result.minimum_block_critical_crack_length
    assert result.final_fracture_utilization > 1.0
    assert result.final_k_max > TOU.k_ic


def test_canonical_is_deterministic():
    first = _canonical()
    second = _canonical()
    assert first.total_cycles == second.total_cycles
    assert first.final_crack_length == second.final_crack_length
    assert first.status == second.status


# --------------------------------------------------------------------------
# AB. a block arrested at a0 activates later
# --------------------------------------------------------------------------
def test_low_block_starts_arrested_and_activates_later():
    result = _canonical()
    low = next(c for c in result.contributions if c.block_name == "A low-amplitude")
    assert result.history[0].advance.block_name == "A low-amplitude"
    assert result.history[0].advance.crack_extension == 0.0
    assert not result.history[0].advance.active_at_start
    assert low.first_active_repeat is not None
    assert low.first_active_repeat > 1
    assert low.crack_extension > 0.0


# --------------------------------------------------------------------------
# AD/AE/AF. growth contributions
# --------------------------------------------------------------------------
def test_contributions_sum_to_the_total_crack_extension():
    result = _canonical()
    assert sum(c.crack_extension for c in result.contributions) == pytest.approx(
        result.total_crack_extension, rel=1e-9
    )


def test_contribution_fractions_sum_to_one():
    result = _canonical()
    assert sum(c.extension_fraction for c in result.contributions) == pytest.approx(
        1.0, rel=1e-9
    )


def test_contribution_cycles_sum_to_the_total():
    result = _canonical()
    assert sum(c.executed_cycles for c in result.contributions) == pytest.approx(
        result.total_cycles, rel=1e-9
    )


def test_only_the_fracture_block_is_tagged():
    result = _canonical()
    tagged = [c.block_name for c in result.contributions if c.triggered_fracture]
    assert tagged == ["C severe gust"]


def test_contributions_cover_every_block_exactly_once():
    result = _canonical()
    assert sorted(c.block_name for c in result.contributions) == sorted(
        CANONICAL_SPECTRUM.block_names
    )


# --------------------------------------------------------------------------
# Y/Z. arrested spectrum
# --------------------------------------------------------------------------
def test_all_blocks_below_threshold_gives_an_arrested_spectrum():
    high = CrackGrowthThreshold(
        name="high",
        delta_k_threshold=mpa_sqrt_m_to_pa_sqrt_m(20.0),
        source_note="ILLUSTRATIVE (test)",
    )
    result = _sim(threshold=high)
    assert result.status is SpectrumStatus.SPECTRUM_ARRESTED
    assert result.final_crack_length == result.initial_crack_length
    assert result.total_crack_extension == 0.0
    assert result.cycles_to_fracture == math.inf


def test_arrest_is_detected_after_a_single_unchanged_repeat():
    """No looping to the repeat guard: one unchanged pass proves it."""
    high = CrackGrowthThreshold(
        name="high",
        delta_k_threshold=mpa_sqrt_m_to_pa_sqrt_m(20.0),
        source_note="ILLUSTRATIVE (test)",
    )
    result = _sim(threshold=high, max_spectrum_repeats=100_000)
    assert result.completed_full_spectra == 1
    assert len(result.history) == len(CANONICAL_SPECTRUM)
    assert result.total_cycles == pytest.approx(
        float(CANONICAL_SPECTRUM.cycles_per_spectrum)
    )


def test_arrested_cycles_still_count_as_applied_cycles():
    high = CrackGrowthThreshold(
        name="high",
        delta_k_threshold=mpa_sqrt_m_to_pa_sqrt_m(20.0),
        source_note="ILLUSTRATIVE (test)",
    )
    result = _sim(threshold=high)
    assert result.total_cycles > 0.0
    assert result.cycles_to_fracture == math.inf


def test_small_initial_flaw_arrests_the_whole_spectrum():
    result = _sim(a0=0.25e-3)
    assert result.status is SpectrumStatus.SPECTRUM_ARRESTED
    assert result.cycles_to_fracture == math.inf


# --------------------------------------------------------------------------
# AA. mixed active/arrested
# --------------------------------------------------------------------------
def test_mixed_spectrum_low_contributes_nothing_until_it_activates():
    result = _sim(max_spectrum_repeats=5)
    first_pass = result.history[:3]
    assert first_pass[0].advance.crack_extension == 0.0  # low, arrested
    assert first_pass[1].advance.crack_extension > 0.0  # manoeuvre
    assert first_pass[2].advance.crack_extension > 0.0  # severe
    assert result.final_crack_length > result.initial_crack_length


# --------------------------------------------------------------------------
# AH. threshold disabled
# --------------------------------------------------------------------------
def test_disabling_the_threshold_shortens_the_predicted_life():
    with_threshold = _canonical().total_cycles
    without = _sim(threshold=None).total_cycles
    assert without < with_threshold
    assert without == pytest.approx(4_018_100.0, rel=1e-9)
    assert 1.0 - without / with_threshold == pytest.approx(0.360, rel=1e-2)


def test_disabled_threshold_makes_every_block_active_from_the_first_pass():
    result = _sim(threshold=None, max_spectrum_repeats=1)
    assert all(h.advance.crack_extension > 0.0 for h in result.history)
    assert all(h.advance.active_at_start for h in result.history)


# --------------------------------------------------------------------------
# AI. repeat guard
# --------------------------------------------------------------------------
def test_max_repeat_guard_stops_without_extrapolating():
    result = _sim(max_spectrum_repeats=10)
    assert result.status is SpectrumStatus.MAX_REPEATS_REACHED
    assert result.completed_full_spectra == 10
    assert result.total_cycles == pytest.approx(
        10.0 * CANONICAL_SPECTRUM.cycles_per_spectrum
    )
    assert result.cycles_to_fracture == math.inf
    assert result.fracture_block_name is None


# --------------------------------------------------------------------------
# T. later blocks are not executed after fracture
# --------------------------------------------------------------------------
def test_no_blocks_execute_after_fracture():
    spectrum = LoadSpectrum(
        blocks=(
            SpectrumBlock(
                name="severe",
                stress_cycle=StressCycle(sigma_max=150.0e6, sigma_min=20.0e6),
                cycle_count=10_000_000,
            ),
            SpectrumBlock(
                name="never-runs",
                stress_cycle=StressCycle(sigma_max=60.0e6, sigma_min=20.0e6),
                cycle_count=100,
            ),
        ),
        name="fracture-first",
        source_note="ILLUSTRATIVE (test fixture)",
    )
    result = _sim(spectrum=spectrum)
    assert result.status is SpectrumStatus.FRACTURE_REACHED
    assert result.fracture_block_name == "severe"
    assert [h.advance.block_name for h in result.history] == ["severe"]
    never = next(c for c in result.contributions if c.block_name == "never-runs")
    assert never.executions == 0
    assert never.executed_cycles == 0.0


# --------------------------------------------------------------------------
# INITIAL_FLAW_ALREADY_CRITICAL
# --------------------------------------------------------------------------
def test_initial_flaw_already_critical():
    a_c = finite_width_critical_crack_length(
        CANONICAL_SPECTRUM[0].stress_cycle, G, TOU
    ).critical_crack_length
    spectrum = LoadSpectrum(
        blocks=(CANONICAL_SPECTRUM[0],),
        name="low only",
        source_note="ILLUSTRATIVE (test fixture)",
    )
    result = _sim(a0=a_c * 1.05, spectrum=spectrum)
    assert result.status is SpectrumStatus.INITIAL_FLAW_ALREADY_CRITICAL
    assert result.total_cycles == 0.0
    assert result.cycles_to_fracture == 0.0


# --------------------------------------------------------------------------
# AT. non-tensile spectrum
# --------------------------------------------------------------------------
def test_non_tensile_spectrum_has_no_fracture_boundary():
    spectrum = LoadSpectrum(
        blocks=(
            SpectrumBlock(
                name="compressive",
                stress_cycle=StressCycle(sigma_max=-10.0e6, sigma_min=-200.0e6),
                cycle_count=1_000,
            ),
        ),
        name="compression only",
        source_note="ILLUSTRATIVE (test fixture)",
    )
    assert spectrum[0].stress_range > 0.0  # algebraically above threshold
    result = _sim(spectrum=spectrum)
    assert result.status is SpectrumStatus.NO_TENSILE_FRACTURE_BOUNDARY
    assert result.cycles_to_fracture == math.inf
    assert result.fracture_block_name is None
    assert result.minimum_block_critical_crack_length == math.inf


# --------------------------------------------------------------------------
# AL. no Miner's rule anywhere
# --------------------------------------------------------------------------
def test_no_miner_damage_sum_is_exposed_as_a_life():
    """Guard against a damage-sum interpretation creeping in. The reported life
    must equal the summed executed cycles, not any n/N ratio."""
    result = _sim(max_spectrum_repeats=4)
    assert result.total_cycles == pytest.approx(
        sum(h.advance.completed_cycles for h in result.history)
    )
    # The contribution fractions are crack extension, not damage: they are
    # driven by growth, and an arrested block contributes exactly zero.
    first_pass_low = result.history[0].advance
    assert first_pass_low.completed_cycles > 0.0
    assert first_pass_low.crack_extension == 0.0


def test_a_miner_style_ratio_would_give_a_different_answer():
    """Sanity contrast: cycles applied are NOT proportional to extension."""
    result = _canonical()
    low = next(c for c in result.contributions if c.block_name == "A low-amplitude")
    severe = next(c for c in result.contributions if c.block_name == "C severe gust")
    cycle_ratio = low.executed_cycles / severe.executed_cycles
    extension_ratio = low.crack_extension / severe.crack_extension
    # The low block runs 100x more cycles (slightly more than 100x, because the
    # severe block's final execution fractures at its very first cycle and so
    # completes none of its 100).
    assert cycle_ratio == pytest.approx(100.177, rel=1e-3)
    # It does NOT produce 100x the crack extension. A cycle-count-weighted
    # damage sum would badly misrank the blocks.
    assert extension_ratio == pytest.approx(4.606, rel=1e-3)
    assert extension_ratio < cycle_ratio / 20.0


# --------------------------------------------------------------------------
# diagnostics
# --------------------------------------------------------------------------
def test_minimum_block_critical_crack_length_is_the_severe_block():
    envelope = minimum_block_critical_crack_length(CANONICAL_SPECTRUM, G, TOU)
    severe = finite_width_critical_crack_length(
        CANONICAL_SPECTRUM[2].stress_cycle, G, TOU
    ).critical_crack_length
    assert envelope == pytest.approx(severe, rel=1e-12)
    assert envelope == pytest.approx(11.8589e-3, rel=1e-4)


def test_crack_may_legitimately_exceed_the_envelope_between_severe_blocks():
    """The envelope is a diagnostic only: the simulation must not terminate
    merely because the crack passed another block's boundary."""
    result = _canonical()
    envelope = result.minimum_block_critical_crack_length
    exceeded = [
        h for h in result.history if h.advance.end_crack_length > envelope
    ]
    assert exceeded, "expected the crack to pass the envelope before fracture"
    # ... and those executions were milder blocks that did not fracture
    assert any(not h.advance.fracture_occurred for h in exceeded)


def test_block_boundary_table_reports_per_block_values():
    rows = block_boundary_table(CANONICAL_SPECTRUM, G, TOU, TH)
    assert [r["name"] for r in rows] == list(CANONICAL_SPECTRUM.block_names)
    assert [r["critical_crack_length"] for r in rows] == pytest.approx(
        [0.03173584, 0.01940922, 0.01185894], rel=1e-5
    )
    assert [r["threshold_crack_length"] for r in rows] == pytest.approx(
        [2.0330e-3, 0.6286e-3, 0.3013e-3], rel=1e-3
    )


def test_block_boundary_table_without_a_threshold():
    rows = block_boundary_table(CANONICAL_SPECTRUM, G, TOU, None)
    assert all(r["threshold_crack_length"] is None for r in rows)


# --------------------------------------------------------------------------
# BF. validation
# --------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [0.0, -1.0e-3, math.nan, math.inf])
def test_invalid_initial_crack_rejected(bad):
    with pytest.raises(ValueError):
        _sim(a0=bad)


def test_wrong_types_and_guards_rejected():
    with pytest.raises(TypeError):
        simulate_repeated_spectrum(1.0e-3, "spectrum", G, LAW, TOU, TH)
    with pytest.raises(TypeError):
        simulate_repeated_spectrum(
            1.0e-3, CANONICAL_SPECTRUM, G, LAW, TOU, TH, max_spectrum_repeats=True
        )
    with pytest.raises(ValueError, match="at least 1"):
        _sim(max_spectrum_repeats=0)
    with pytest.raises(ValueError, match="at least 1"):
        _sim(max_history_records=0)


def test_history_cap_does_not_change_the_answer():
    full = _sim(max_spectrum_repeats=20)
    capped = _sim(max_spectrum_repeats=20, max_history_records=5)
    assert capped.total_cycles == full.total_cycles
    assert capped.final_crack_length == full.final_crack_length
    assert capped.history_truncated
    assert len(capped.history) == 5
    for a, b in zip(full.contributions, capped.contributions):
        assert a.crack_extension == pytest.approx(b.crack_extension)
        assert a.executed_cycles == pytest.approx(b.executed_cycles)


def test_final_k_max_matches_the_fracture_block():
    result = _canonical()
    expected = k_max_for_geometry(
        result.final_crack_length,
        CANONICAL_SPECTRUM[result.fracture_block_index].stress_cycle,
        G,
    )
    assert result.final_k_max == pytest.approx(expected, rel=1e-12)
    assert result.final_fracture_utilization == pytest.approx(
        expected / TOU.k_ic, rel=1e-12
    )
