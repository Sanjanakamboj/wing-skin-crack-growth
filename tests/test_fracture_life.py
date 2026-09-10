"""U-Z, AC-AG, AN. Life from the initial flaw to the fracture boundary."""

from __future__ import annotations

import math

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GEOMETRY,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    FlawAdmissibility,
    FractureToughness,
    InitialFlawBeyondCriticalError,
    StressCycle,
    ThroughCrackGeometry,
    analytical_cycles_to_crack_length,
    assess_initial_flaw,
    crack_growth_rate,
    critical_crack_length,
    cycles_to_critical_crack,
    mpa_sqrt_m_to_pa_sqrt_m,
)

UNIT = ThroughCrackGeometry(geometry_factor=1.0)
K30 = FractureToughness(
    name="test", k_ic=30.0e6, source_note="ILLUSTRATIVE (test fixture)"
)


def _life(
    a0=1.0e-3,
    cycle=CANONICAL_CYCLE,
    geometry=UNIT,
    toughness=K30,
    law=CANONICAL_PARIS_LAW,
    **kwargs,
):
    return cycles_to_critical_crack(
        a0, cycle, geometry, law, toughness, **kwargs
    )


# --------------------------------------------------------------------------
# U/V/W. initial-flaw admissibility
# --------------------------------------------------------------------------
def test_initial_flaw_below_critical_gives_finite_life():
    result = _life(a0=1.0e-3)
    assert result.admissibility is FlawAdmissibility.BELOW_CRITICAL
    assert result.predicted_cycles > 0.0
    assert math.isfinite(result.predicted_cycles)
    assert result.has_finite_life
    assert result.growth_result is not None


def test_initial_flaw_exactly_at_critical_gives_zero_life():
    a_c = critical_crack_length(
        CANONICAL_CYCLE, UNIT, K30
    ).critical_crack_length
    result = _life(a0=a_c)
    assert result.admissibility is FlawAdmissibility.AT_CRITICAL
    assert result.predicted_cycles == 0.0
    assert result.analytical_cycles == 0.0
    assert result.growth_result is None


def test_initial_flaw_above_critical_is_reported_honestly():
    a_c = critical_crack_length(
        CANONICAL_CYCLE, UNIT, K30
    ).critical_crack_length
    result = _life(a0=a_c * 1.5)
    assert result.admissibility is FlawAdmissibility.ABOVE_CRITICAL
    assert result.predicted_cycles == 0.0
    assert result.growth_result is None
    # the flaw is genuinely outside the boundary, not merely at it
    assert not result.initial_assessment.passes
    assert result.initial_assessment.utilization > 1.0
    assert result.initial_assessment.margin < 0.0


def test_initial_flaw_above_critical_raises_in_strict_mode():
    a_c = critical_crack_length(
        CANONICAL_CYCLE, UNIT, K30
    ).critical_crack_length
    with pytest.raises(InitialFlawBeyondCriticalError, match="already exceeds"):
        _life(a0=a_c * 1.5, strict=True)


def test_strict_mode_does_not_affect_admissible_flaws():
    assert _life(a0=1.0e-3, strict=True).predicted_cycles == pytest.approx(
        _life(a0=1.0e-3).predicted_cycles
    )


def test_integration_never_runs_backwards():
    """For a0 > a_c the integrator must not be called at all."""
    a_c = critical_crack_length(
        CANONICAL_CYCLE, UNIT, K30
    ).critical_crack_length
    result = _life(a0=a_c * 2.0)
    assert result.growth_result is None
    assert result.predicted_cycles >= 0.0


def test_assess_initial_flaw_classifies_directly():
    critical = critical_crack_length(CANONICAL_CYCLE, UNIT, K30)
    a_c = critical.critical_crack_length
    assert assess_initial_flaw(a_c * 0.5, critical) is FlawAdmissibility.BELOW_CRITICAL
    assert assess_initial_flaw(a_c, critical) is FlawAdmissibility.AT_CRITICAL
    assert assess_initial_flaw(a_c * 2.0, critical) is FlawAdmissibility.ABOVE_CRITICAL


def test_no_tensile_boundary_gives_infinite_life():
    cycle = StressCycle(sigma_max=-10.0e6, sigma_min=-100.0e6)
    result = _life(cycle=cycle)
    assert result.admissibility is FlawAdmissibility.NO_TENSILE_BOUNDARY
    assert result.critical_crack_length == math.inf
    assert result.predicted_cycles == math.inf
    assert result.analytical_cycles == math.inf
    assert result.growth_result is None
    assert result.critical_assessment is None
    assert result.delta_k_at_critical is None


def test_zero_stress_range_gives_infinite_life_despite_finite_boundary():
    """A finite fracture size exists, but the crack never grows to it."""
    cycle = StressCycle(sigma_max=120.0e6, sigma_min=120.0e6)
    result = _life(cycle=cycle)
    assert result.admissibility is FlawAdmissibility.BELOW_CRITICAL
    assert math.isfinite(result.critical_crack_length)
    assert result.predicted_cycles == math.inf
    assert result.analytical_cycles == math.inf


# --------------------------------------------------------------------------
# Numerical edge cases (item 25)
# --------------------------------------------------------------------------
def test_a_zero_marginally_below_critical_gives_small_positive_life():
    a_c = critical_crack_length(
        CANONICAL_CYCLE, UNIT, K30
    ).critical_crack_length
    for fraction in (1.0 - 1.0e-6, 1.0 - 1.0e-9, 1.0 - 1.0e-12):
        result = _life(a0=a_c * fraction)
        assert result.admissibility is FlawAdmissibility.BELOW_CRITICAL
        assert result.predicted_cycles > 0.0, fraction
        assert math.isfinite(result.predicted_cycles)


def test_marginal_life_shrinks_toward_zero_as_a_zero_approaches_critical():
    a_c = critical_crack_length(
        CANONICAL_CYCLE, UNIT, K30
    ).critical_crack_length
    lives = [
        _life(a0=a_c * f).predicted_cycles
        for f in (1.0 - 1.0e-3, 1.0 - 1.0e-6, 1.0 - 1.0e-9)
    ]
    assert all(later < earlier for earlier, later in zip(lives, lives[1:]))
    assert all(v > 0.0 for v in lives)


# --------------------------------------------------------------------------
# X/Y/Z. numerical life, agreement, determinism
# --------------------------------------------------------------------------
@pytest.mark.parametrize("k_ic_mpa", [20.0, 25.0, 30.0, 40.0, 60.0])
@pytest.mark.parametrize("a0", [0.25e-3, 1.0e-3, 4.0e-3])
def test_numerical_life_to_critical_is_positive(k_ic_mpa, a0):
    toughness = FractureToughness(
        name="s",
        k_ic=mpa_sqrt_m_to_pa_sqrt_m(k_ic_mpa),
        source_note="ILLUSTRATIVE (test)",
    )
    result = _life(a0=a0, toughness=toughness)
    assert result.predicted_cycles > 0.0
    assert math.isfinite(result.predicted_cycles)


@pytest.mark.parametrize("k_ic_mpa", [20.0, 30.0, 45.0, 60.0])
@pytest.mark.parametrize("sigma_max", [80.0e6, 120.0e6, 200.0e6])
@pytest.mark.parametrize("y", [0.8, 1.0, 1.2])
def test_numerical_matches_analytical_life_to_fracture(k_ic_mpa, sigma_max, y):
    toughness = FractureToughness(
        name="s",
        k_ic=mpa_sqrt_m_to_pa_sqrt_m(k_ic_mpa),
        source_note="ILLUSTRATIVE (test)",
    )
    cycle = StressCycle(sigma_max=sigma_max, sigma_min=0.2 * sigma_max)
    result = cycles_to_critical_crack(
        0.5e-3,
        cycle,
        ThroughCrackGeometry(y),
        CANONICAL_PARIS_LAW,
        toughness,
        intervals=20000,
    )
    # This sweep deliberately includes extreme spans (up to a_c/a0 ~ 560, i.e.
    # a 280 mm critical crack -- far outside any sensible wing skin, but a
    # worthwhile stress test of the quadrature). The fixed uniform grid needs
    # more intervals as that span grows, so 20000 intervals and a 1e-7 tolerance
    # are used here; see test_required_intervals_grow_with_crack_size_span.
    assert result.predicted_cycles == pytest.approx(
        result.analytical_cycles, rel=1e-7
    )


def test_required_intervals_grow_with_crack_size_span():
    """The fixed uniform Simpson grid resolves 1/(da/dN) less well as a_c/a0
    grows, because the integrand is steepest just above a0.

    This is a documented property of the Milestone 1 quadrature, not an error:
    the convergence order is unaffected. A wide-span case simply needs more
    intervals to reach a given tolerance, so callers integrating to a distant
    fracture boundary should raise ``intervals`` accordingly.
    """
    toughness = FractureToughness(
        name="wide",
        k_ic=mpa_sqrt_m_to_pa_sqrt_m(60.0),
        source_note="ILLUSTRATIVE (test)",
    )
    cycle = StressCycle(sigma_max=80.0e6, sigma_min=16.0e6)
    geometry = ThroughCrackGeometry(1.2)
    a0 = 0.5e-3
    a_c = critical_crack_length(cycle, geometry, toughness).critical_crack_length
    assert a_c / a0 > 200.0

    analytical = analytical_cycles_to_crack_length(
        a0, a_c, cycle, geometry, CANONICAL_PARIS_LAW
    )

    def relative_error(intervals: int) -> float:
        n = cycles_to_critical_crack(
            a0, cycle, geometry, CANONICAL_PARIS_LAW, toughness, intervals=intervals
        ).predicted_cycles
        return abs(n - analytical) / analytical

    coarse, fine = relative_error(1000), relative_error(20000)
    # even the coarse grid is accurate to better than 0.1 % ...
    assert coarse < 1.0e-3
    # ... and refinement still converges at the expected fourth order.
    assert fine < 1.0e-8
    assert relative_error(2000) / relative_error(4000) > 10.0


def test_narrow_span_is_accurate_at_the_default_interval_count():
    """The canonical case (a_c/a0 ~ 20) needs no special interval count."""
    result = _life(a0=1.0e-3)
    assert result.critical_crack_length / result.initial_crack_length < 25.0
    assert result.predicted_cycles == pytest.approx(
        result.analytical_cycles, rel=1e-8
    )


def test_analytical_reference_is_the_unchanged_m1_closed_form():
    """The analytical life must equal the Milestone 1 closed form evaluated at
    a_final = a_critical -- no separate fracture-growth equation exists."""
    a_c = critical_crack_length(
        CANONICAL_CYCLE, UNIT, K30
    ).critical_crack_length
    expected = analytical_cycles_to_crack_length(
        1.0e-3, a_c, CANONICAL_CYCLE, UNIT, CANONICAL_PARIS_LAW
    )
    assert _life(a0=1.0e-3).analytical_cycles == expected


def test_repeated_life_calculation_is_bit_for_bit_identical():
    first = _life(a0=1.0e-3, intervals=512)
    for _ in range(5):
        assert _life(a0=1.0e-3, intervals=512) == first


# --------------------------------------------------------------------------
# AC. life increases with toughness
# --------------------------------------------------------------------------
def test_life_and_critical_size_increase_with_toughness():
    lives, sizes = [], []
    for k_ic_mpa in (20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0):
        result = _life(
            toughness=FractureToughness(
                name="s",
                k_ic=mpa_sqrt_m_to_pa_sqrt_m(k_ic_mpa),
                source_note="ILLUSTRATIVE (test)",
            )
        )
        lives.append(result.predicted_cycles)
        sizes.append(result.critical_crack_length)
    assert all(later > earlier for earlier, later in zip(sizes, sizes[1:]))
    assert all(later > earlier for earlier, later in zip(lives, lives[1:]))


def test_life_gain_from_toughness_has_diminishing_returns():
    """Because most life is spent at small crack sizes, doubling K_IC does not
    come close to doubling the life."""
    low = _life(
        toughness=FractureToughness(
            name="l", k_ic=mpa_sqrt_m_to_pa_sqrt_m(25.0), source_note="ILLUSTRATIVE"
        )
    ).predicted_cycles
    high = _life(
        toughness=FractureToughness(
            name="h", k_ic=mpa_sqrt_m_to_pa_sqrt_m(50.0), source_note="ILLUSTRATIVE"
        )
    ).predicted_cycles
    assert high > low
    assert high < 1.3 * low


# --------------------------------------------------------------------------
# AD/AE. sigma_max at FIXED delta_sigma
# --------------------------------------------------------------------------
def _fixed_range_cycle(sigma_max: float) -> StressCycle:
    return StressCycle(sigma_max=sigma_max, sigma_min=sigma_max - 100.0e6)


SIGMA_MAX_SWEEP = (100.0e6, 120.0e6, 140.0e6, 160.0e6, 180.0e6)


def test_growth_rate_unchanged_when_sigma_max_varies_at_fixed_range():
    """delta_K -- and hence Paris growth -- depends only on delta_sigma."""
    a = 2.0e-3
    rates = [
        crack_growth_rate(a, _fixed_range_cycle(s), UNIT, CANONICAL_PARIS_LAW)
        for s in SIGMA_MAX_SWEEP
    ]
    for rate in rates[1:]:
        assert rate == pytest.approx(rates[0], rel=1e-12)


def test_critical_size_shrinks_as_sigma_max_rises_at_fixed_range():
    sizes = [
        critical_crack_length(
            _fixed_range_cycle(s), UNIT, K30
        ).critical_crack_length
        for s in SIGMA_MAX_SWEEP
    ]
    assert all(later < earlier for earlier, later in zip(sizes, sizes[1:]))


def test_life_decreases_with_sigma_max_at_fixed_range():
    """The life falls purely because the endpoint moves inward, since the
    growth rate at every crack length is identical across this sweep."""
    lives = [
        _life(cycle=_fixed_range_cycle(s)).predicted_cycles for s in SIGMA_MAX_SWEEP
    ]
    assert all(later < earlier for earlier, later in zip(lives, lives[1:]))


def test_stress_ratio_changes_across_the_fixed_range_sweep():
    ratios = [_fixed_range_cycle(s).stress_ratio for s in SIGMA_MAX_SWEEP]
    assert all(later > earlier for earlier, later in zip(ratios, ratios[1:]))


# --------------------------------------------------------------------------
# AF. life decreases with Y
# --------------------------------------------------------------------------
def test_life_and_critical_size_decrease_with_geometry_factor():
    lives, sizes = [], []
    for y in (0.8, 0.9, 1.0, 1.1, 1.2):
        result = _life(geometry=ThroughCrackGeometry(y))
        lives.append(result.predicted_cycles)
        sizes.append(result.critical_crack_length)
    assert all(later < earlier for earlier, later in zip(sizes, sizes[1:]))
    assert all(later < earlier for earlier, later in zip(lives, lives[1:]))


def test_critical_size_follows_inverse_square_y_within_the_sweep():
    base = _life(geometry=ThroughCrackGeometry(1.0)).critical_crack_length
    for y in (0.8, 0.9, 1.1, 1.2):
        assert _life(
            geometry=ThroughCrackGeometry(y)
        ).critical_crack_length == pytest.approx(base * y**-2, rel=1e-12)


# --------------------------------------------------------------------------
# AG. life decreases with larger initial flaw
# --------------------------------------------------------------------------
def test_life_decreases_with_larger_initial_flaw():
    lives = [
        _life(a0=a0).predicted_cycles
        for a0 in (0.25e-3, 0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3)
    ]
    assert all(later < earlier for earlier, later in zip(lives, lives[1:]))


def test_initial_flaw_fraction_and_utilization_rise_together():
    a_c = critical_crack_length(
        CANONICAL_CYCLE, UNIT, K30
    ).critical_crack_length
    previous = 0.0
    for a0 in (0.25e-3, 0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3):
        result = _life(a0=a0)
        utilization = result.initial_assessment.utilization
        assert utilization > previous
        # utilization = sqrt(a0 / a_c) for constant Y
        assert utilization == pytest.approx(math.sqrt(a0 / a_c), rel=1e-12)
        previous = utilization


# --------------------------------------------------------------------------
# AN. result fields independently verified
# --------------------------------------------------------------------------
def test_canonical_life_to_fracture_result_fields():
    result = cycles_to_critical_crack(
        CANONICAL_INITIAL_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
        CANONICAL_FRACTURE_TOUGHNESS,
    )
    assert result.initial_crack_length == 1.0e-3
    assert result.critical_crack_length == pytest.approx(
        0.019894367886486918, rel=1e-12
    )
    assert result.admissibility is FlawAdmissibility.BELOW_CRITICAL

    # analytical life, hand-checkable via the M1 closed form at a_f = a_c
    assert result.analytical_cycles == pytest.approx(881160.779753657, rel=1e-9)
    assert result.predicted_cycles == pytest.approx(881160.7850256478, rel=1e-9)
    assert result.predicted_cycles == pytest.approx(
        result.analytical_cycles, rel=1e-8
    )

    # endpoints carried through from the unchanged M1 integrator
    growth = result.growth_result
    assert growth is not None
    assert growth.initial_crack_length == 1.0e-3
    assert growth.final_crack_length == pytest.approx(result.critical_crack_length)
    assert growth.initial_delta_k == pytest.approx(5.604991216397929e6, rel=1e-12)
    assert growth.initial_growth_rate == pytest.approx(
        1.7608599228871055e-09, rel=1e-10
    )
    assert growth.final_delta_k == pytest.approx(25.0e6, rel=1e-9)
    assert growth.final_growth_rate == pytest.approx(1.5625e-07, rel=1e-9)

    # fracture endpoint
    assert result.critical_assessment is not None
    assert result.critical_assessment.k_max == pytest.approx(30.0e6, rel=1e-12)
    assert result.critical_assessment.utilization == pytest.approx(1.0, rel=1e-12)
    assert result.delta_k_at_critical == pytest.approx(25.0e6, rel=1e-12)

    # initial screen
    assert result.initial_assessment.k_max == pytest.approx(
        6.725989459677515e6, rel=1e-12
    )
    assert result.initial_assessment.margin == pytest.approx(
        3.460310290381928, rel=1e-12
    )
    assert result.initial_assessment.passes


def test_life_to_fracture_rejects_bad_input():
    with pytest.raises(ValueError, match="strictly positive"):
        _life(a0=0.0)
    with pytest.raises(ValueError, match="finite"):
        _life(a0=math.nan)
    with pytest.raises(TypeError):
        cycles_to_critical_crack(
            1.0e-3, CANONICAL_CYCLE, UNIT, "law", K30
        )
    with pytest.raises(TypeError):
        cycles_to_critical_crack(
            1.0e-3, CANONICAL_CYCLE, UNIT, CANONICAL_PARIS_LAW, "toughness"
        )
    with pytest.raises(ValueError, match="even"):
        _life(intervals=101)
