"""G-R, U, AC, AW, AX, AY. Single-block advancement."""

from __future__ import annotations

import math

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GROWTH_THRESHOLD,
    CANONICAL_PARIS_LAW,
    BlockAdvanceStatus,
    CrackGrowthThreshold,
    FractureToughness,
    SpectrumBlock,
    StressCycle,
    advance_crack_through_block,
    delta_k_for_geometry,
    finite_width_critical_crack_length,
    finite_width_threshold_crack_length,
    integrate_crack_growth_log_grid,
    k_max_for_geometry,
    mpa_sqrt_m_to_pa_sqrt_m,
)

G = CANONICAL_FINITE_WIDTH_GEOMETRY
LAW = CANONICAL_PARIS_LAW
TOU = CANONICAL_FRACTURE_TOUGHNESS
TH = CANONICAL_GROWTH_THRESHOLD

LOW = SpectrumBlock(
    name="A", stress_cycle=StressCycle(sigma_max=70.0e6, sigma_min=20.0e6),
    cycle_count=10_000,
)
MED = SpectrumBlock(
    name="B", stress_cycle=StressCycle(sigma_max=110.0e6, sigma_min=20.0e6),
    cycle_count=1_000,
)
SEVERE = SpectrumBlock(
    name="C", stress_cycle=StressCycle(sigma_max=150.0e6, sigma_min=20.0e6),
    cycle_count=100,
)


def _advance(a0, block, threshold=TH, **kw):
    return advance_crack_through_block(a0, block, G, LAW, TOU, threshold, **kw)


# --------------------------------------------------------------------------
# G/H/I. per-block boundaries
# --------------------------------------------------------------------------
def test_each_block_has_its_own_threshold_verdict_at_the_same_crack():
    a0 = 1.0e-3
    assert not _advance(a0, LOW).active_at_start
    assert _advance(a0, MED).active_at_start
    assert _advance(a0, SEVERE).active_at_start


def test_each_block_has_its_own_critical_crack_size():
    sizes = [
        _advance(1.0e-3, b).block_critical_crack_length for b in (LOW, MED, SEVERE)
    ]
    assert sizes == pytest.approx([0.03173584, 0.01940922, 0.01185894], rel=1e-5)
    # a higher sigma_max gives a SMALLER fracture boundary
    assert sizes[0] > sizes[1] > sizes[2]


def test_block_boundaries_are_not_a_spectrum_average():
    """Each block keeps its OWN boundary; no single spectrum-wide value is used.

    The comparison uses a CYCLE-WEIGHTED mean sigma_max, which is the averaging
    someone would plausibly reach for. (A plain arithmetic mean of 70/110/150
    happens to equal block B's own 110 MPa, so it would not discriminate.)
    """
    sizes = [
        _advance(1.0e-3, b).block_critical_crack_length for b in (LOW, MED, SEVERE)
    ]
    assert len(set(sizes)) == 3  # three distinct boundaries

    blocks = (LOW, MED, SEVERE)
    total = sum(b.cycle_count for b in blocks)
    weighted_sigma = sum(b.sigma_max * b.cycle_count for b in blocks) / total
    assert weighted_sigma == pytest.approx(74.324e6, rel=1e-3)
    averaged = finite_width_critical_crack_length(
        StressCycle(sigma_max=weighted_sigma, sigma_min=20.0e6), G, TOU
    ).critical_crack_length
    for size in sizes:
        assert size != pytest.approx(averaged, rel=1e-6)
    # and the averaged boundary is far more permissive than the governing one
    assert averaged > min(sizes) * 2.0


# --------------------------------------------------------------------------
# J/K/L. arrested vs active
# --------------------------------------------------------------------------
def test_arrested_block_produces_zero_extension_but_consumes_its_cycles():
    result = _advance(1.0e-3, LOW)
    assert result.status is BlockAdvanceStatus.ARRESTED_BELOW_THRESHOLD
    assert result.crack_extension == 0.0
    assert result.end_crack_length == result.start_crack_length
    assert result.completed_cycles == float(LOW.cycle_count)
    assert result.solver_iterations == 0
    assert not result.fracture_occurred


def test_active_block_advances_the_crack():
    result = _advance(1.0e-3, MED)
    assert result.status is BlockAdvanceStatus.ADVANCED
    assert result.end_crack_length > result.start_crack_length
    assert result.crack_extension > 0.0
    assert result.completed_cycles == float(MED.cycle_count)
    assert result.active_at_start and result.active_at_end


def test_disabled_threshold_makes_the_low_block_active():
    arrested = _advance(1.0e-3, LOW, threshold=TH)
    active = _advance(1.0e-3, LOW, threshold=None)
    assert arrested.crack_extension == 0.0
    assert active.crack_extension > 0.0
    assert active.status is BlockAdvanceStatus.ADVANCED


# --------------------------------------------------------------------------
# M. the block root round-trips through the verified life integral
# --------------------------------------------------------------------------
@pytest.mark.parametrize("block", [MED, SEVERE])
@pytest.mark.parametrize("a0", [1.0e-3, 3.0e-3, 8.0e-3])
def test_block_end_reproduces_the_requested_cycle_count(block, a0):
    result = _advance(a0, block)
    assert result.status is BlockAdvanceStatus.ADVANCED
    cycles = integrate_crack_growth_log_grid(
        result.start_crack_length,
        result.end_crack_length,
        block.stress_cycle,
        G,
        LAW,
        intervals=200,
    ).predicted_cycles
    assert cycles == pytest.approx(float(block.cycle_count), rel=1e-6)


def test_block_growth_exceeds_the_explicit_euler_estimate():
    """da/dN increases with a, so an explicit Euler step under-predicts growth.
    This is exactly the property the bracket seeding relies on."""
    a0 = 1.0e-3
    result = _advance(a0, MED)
    euler = MED.cycle_count * LAW.growth_rate(
        delta_k_for_geometry(a0, MED.stress_cycle, G)
    )
    assert result.crack_extension > euler
    assert result.crack_extension == pytest.approx(euler, rel=1e-2)


# --------------------------------------------------------------------------
# N/O/P. solver determinism and bounds
# --------------------------------------------------------------------------
def test_block_solver_is_deterministic():
    first = _advance(1.0e-3, MED)
    for _ in range(4):
        assert _advance(1.0e-3, MED) == first


def test_block_end_stays_within_the_bracket():
    result = _advance(1.0e-3, SEVERE)
    assert result.start_crack_length < result.end_crack_length
    assert result.end_crack_length < result.block_critical_crack_length
    assert result.end_crack_length < G.max_admissible_crack_length


def test_block_never_evaluates_at_half_width():
    """The panel rejects a >= W/2, so a probe there would raise, not return."""
    result = _advance(1.0e-3, SEVERE, tolerance=1.0e-15, max_iterations=400)
    assert math.isfinite(result.end_geometry_factor)
    assert result.end_crack_length < G.max_admissible_crack_length


def test_tighter_tolerance_refines_the_end_crack_length():
    coarse = _advance(1.0e-3, MED, tolerance=1.0e-9)
    fine = _advance(1.0e-3, MED, tolerance=1.0e-15)
    assert fine.end_crack_length == pytest.approx(
        coarse.end_crack_length, abs=1.0e-9
    )


# --------------------------------------------------------------------------
# AX/AY. numerical convergence
# --------------------------------------------------------------------------
@pytest.mark.parametrize("intervals", [50, 100, 200, 500, 1000])
def test_end_crack_length_is_insensitive_to_interval_count(intervals):
    """The per-block span is short, so the log-grid rule is already exact
    there; the solve tolerance, not the interval count, sets the accuracy."""
    reference = _advance(1.0e-3, MED, intervals=2000).end_crack_length
    assert _advance(1.0e-3, MED, intervals=intervals).end_crack_length == (
        pytest.approx(reference, rel=1e-9)
    )


@pytest.mark.parametrize("tolerance", [1.0e-6, 1.0e-8, 1.0e-10, 1.0e-12])
def test_end_crack_length_converges_with_tolerance(tolerance):
    reference = _advance(1.0e-3, MED, tolerance=1.0e-15).end_crack_length
    assert _advance(1.0e-3, MED, tolerance=tolerance).end_crack_length == (
        pytest.approx(reference, abs=max(tolerance, 1e-15))
    )


# --------------------------------------------------------------------------
# Q/R. fracture inside a block
# --------------------------------------------------------------------------
def test_fracture_inside_a_block_reports_a_partial_cycle_count():
    """Start just below the severe block's boundary with far too many cycles."""
    a_c = finite_width_critical_crack_length(
        SEVERE.stress_cycle, G, TOU
    ).critical_crack_length
    block = SpectrumBlock(
        name="C-long", stress_cycle=SEVERE.stress_cycle, cycle_count=10_000_000
    )
    result = _advance(a_c * 0.9, block)
    assert result.status is BlockAdvanceStatus.FRACTURED_WITHIN_BLOCK
    assert result.fracture_occurred
    assert 0.0 < result.completed_cycles < float(block.cycle_count)
    assert result.fracture_cycle_within_block == pytest.approx(
        result.completed_cycles
    )
    assert result.end_crack_length == pytest.approx(a_c, rel=1e-9)


def test_fracture_end_state_sits_on_the_boundary():
    a_c = finite_width_critical_crack_length(
        SEVERE.stress_cycle, G, TOU
    ).critical_crack_length
    block = SpectrumBlock(
        name="C-long", stress_cycle=SEVERE.stress_cycle, cycle_count=10_000_000
    )
    result = _advance(a_c * 0.9, block)
    k_max = k_max_for_geometry(result.end_crack_length, SEVERE.stress_cycle, G)
    assert k_max / TOU.k_ic == pytest.approx(1.0, rel=1e-8)
    assert result.end_fracture_utilization == pytest.approx(1.0, rel=1e-8)


def test_completed_cycles_are_not_rounded_to_the_block_end():
    a_c = finite_width_critical_crack_length(
        SEVERE.stress_cycle, G, TOU
    ).critical_crack_length
    block = SpectrumBlock(
        name="C-long", stress_cycle=SEVERE.stress_cycle, cycle_count=10_000_000
    )
    result = _advance(a_c * 0.99, block)
    assert result.completed_cycles != float(block.cycle_count)
    assert isinstance(result.completed_cycles, float)


# --------------------------------------------------------------------------
# U. fracture detected at the START of a newly encountered severe block
# --------------------------------------------------------------------------
def test_fracture_detected_immediately_when_block_begins_over_its_boundary():
    """A crack grown by earlier blocks may already exceed THIS block's
    boundary; no growth within the block is required to detect it."""
    a_c = finite_width_critical_crack_length(
        SEVERE.stress_cycle, G, TOU
    ).critical_crack_length
    result = _advance(a_c * 1.01, SEVERE)
    assert result.status is BlockAdvanceStatus.ALREADY_AT_OR_BEYOND_FRACTURE
    assert result.fracture_occurred
    assert result.completed_cycles == 0.0
    assert result.fracture_cycle_within_block == 0.0
    assert result.crack_extension == 0.0
    assert result.start_fracture_utilization > 1.0


def test_a_milder_block_survives_a_crack_that_fractures_a_severe_block():
    """The same crack is safe under the low block and critical under the
    severe one -- which is why fracture must be checked per block."""
    a_c_severe = finite_width_critical_crack_length(
        SEVERE.stress_cycle, G, TOU
    ).critical_crack_length
    a = a_c_severe * 1.01
    assert _advance(a, SEVERE).fracture_occurred
    low = _advance(a, LOW)
    assert not low.fracture_occurred
    assert low.start_fracture_utilization < 1.0


# --------------------------------------------------------------------------
# AC. threshold state is recomputed, never cached from a0
# --------------------------------------------------------------------------
def test_arrested_block_becomes_active_at_a_larger_crack():
    a_th = finite_width_threshold_crack_length(
        LOW.stress_cycle, G, TH
    ).threshold_crack_length
    assert a_th == pytest.approx(2.033e-3, rel=1e-3)
    assert not _advance(a_th * 0.99, LOW).active_at_start
    assert _advance(a_th * 1.01, LOW).active_at_start
    assert _advance(a_th * 0.99, LOW).crack_extension == 0.0
    assert _advance(a_th * 1.01, LOW).crack_extension > 0.0


# --------------------------------------------------------------------------
# AV. a non-tensile block cannot fracture
# --------------------------------------------------------------------------
def test_non_tensile_block_grows_but_never_fractures():
    """No closure is modelled, so a compression-only cycle still has a positive
    algebraic delta_K -- but it must not create a tensile fracture."""
    block = SpectrumBlock(
        name="compressive",
        stress_cycle=StressCycle(sigma_max=-10.0e6, sigma_min=-200.0e6),
        cycle_count=1_000,
    )
    result = _advance(1.0e-3, block)
    assert block.stress_range > 0.0
    assert result.block_critical_crack_length == math.inf
    assert not result.fracture_occurred
    assert result.start_fracture_utilization == 0.0
    assert result.status is BlockAdvanceStatus.ADVANCED
    assert result.crack_extension > 0.0


def test_zero_range_block_is_arrested():
    block = SpectrumBlock(
        name="steady",
        stress_cycle=StressCycle(sigma_max=100.0e6, sigma_min=100.0e6),
        cycle_count=500,
    )
    result = _advance(1.0e-3, block)
    assert result.start_delta_k == 0.0
    assert result.crack_extension == 0.0
    assert result.status is BlockAdvanceStatus.ARRESTED_BELOW_THRESHOLD
    assert result.completed_cycles == 500.0


# --------------------------------------------------------------------------
# AW/BF. validation
# --------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [0.0, -1.0e-3, math.nan, math.inf])
def test_invalid_crack_length_rejected(bad):
    with pytest.raises(ValueError):
        _advance(bad, MED)


def test_crack_beyond_half_width_rejected():
    with pytest.raises(ValueError, match="not admissible"):
        _advance(0.06, MED)


def test_wrong_types_rejected():
    with pytest.raises(TypeError):
        advance_crack_through_block(1.0e-3, "block", G, LAW, TOU, TH)
    with pytest.raises(TypeError):
        advance_crack_through_block(1.0e-3, MED, G, "law", TOU, TH)
    with pytest.raises(TypeError):
        advance_crack_through_block(1.0e-3, MED, G, LAW, TOU, "threshold")
    with pytest.raises(ValueError, match="strictly positive"):
        _advance(1.0e-3, MED, tolerance=0.0)
    with pytest.raises(ValueError, match="at least 1"):
        _advance(1.0e-3, MED, max_iterations=0)


def test_high_threshold_arrests_every_block():
    high = CrackGrowthThreshold(
        name="high",
        delta_k_threshold=mpa_sqrt_m_to_pa_sqrt_m(30.0),
        source_note="ILLUSTRATIVE (test)",
    )
    for block in (LOW, MED, SEVERE):
        result = _advance(1.0e-3, block, threshold=high)
        assert result.crack_extension == 0.0


def test_low_toughness_makes_the_start_already_critical():
    weak = FractureToughness(
        name="weak",
        k_ic=mpa_sqrt_m_to_pa_sqrt_m(3.0),
        source_note="ILLUSTRATIVE (test)",
    )
    result = advance_crack_through_block(1.0e-3, SEVERE, G, LAW, weak, TH)
    assert result.fracture_occurred
    assert result.completed_cycles == 0.0
