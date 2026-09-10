"""T-AJ. Analytical reference, numerical quadrature, and their agreement.

Independence: the tests of the analytical closed form use hand-derived literal
expected values, never the numerical integrator; the convergence and
determinism tests exercise the numerical integrator alone. Only the explicit
cross-check tests compare the two, which is the point of that check.
"""

from __future__ import annotations

import math

import pytest

from crackgrowth import (
    CrackGrowthResult,
    ParisLaw,
    StressCycle,
    ThroughCrackGeometry,
    analytical_cycles_to_crack_length,
    crack_growth_rate,
    cycles_to_crack_length,
    delta_stress_intensity,
)


def _law(C: float, m: float) -> ParisLaw:
    return ParisLaw(C=C, m=m, name="test", source_note="ILLUSTRATIVE (test)")


UNIT = ThroughCrackGeometry(geometry_factor=1.0)
CYCLE_100 = StressCycle(sigma_max=120.0e6, sigma_min=20.0e6)


# --------------------------------------------------------------------------
# T. analytical life, hand calculation, m != 2
# --------------------------------------------------------------------------
def test_analytical_life_hand_calculation_m_equals_three():
    """m = 3, C = 1e-29, Y = 1, delta_sigma = 100 MPa, a: 1 mm -> 10 mm.

    p = 1 - 3/2 = -0.5
    af**p - a0**p = 0.01**-0.5 - 0.001**-0.5 = 10 - 31.6227766017 = -21.6227766017
    (Y*ds*sqrt(pi))**3 = (1.772453850905516e8)**3 = 5.568327996831708e24
    denominator = -0.5 * 1e-29 * 5.568327996831708e24 = -2.784163998415854e-5
    N = -21.6227766017 / -2.784163998415854e-5 = 776634.4444503565 cycles
    """
    n = analytical_cycles_to_crack_length(
        1.0e-3, 10.0e-3, CYCLE_100, UNIT, _law(1.0e-29, 3.0)
    )
    assert n == pytest.approx(776634.4444503565, rel=1e-12)


def test_analytical_life_hand_calculation_with_geometry_factor():
    """Y = 1.2, m = 3, C = 1e-29, a: 2 mm -> 8 mm."""
    n = analytical_cycles_to_crack_length(
        2.0e-3, 8.0e-3, CYCLE_100, ThroughCrackGeometry(1.2), _law(1.0e-29, 3.0)
    )
    assert n == pytest.approx(232389.47508994734, rel=1e-12)


def test_analytical_life_hand_calculation_m_below_two():
    """m = 1.5 -> p = +0.25, so numerator and denominator are both positive."""
    C, m = 1.0e-25, 1.5
    a0, af = 1.0e-3, 10.0e-3
    p = 0.25
    expected = (af**p - a0**p) / (
        p * C * (1.0 * 100.0e6 * math.sqrt(math.pi)) ** m
    )
    n = analytical_cycles_to_crack_length(a0, af, CYCLE_100, UNIT, _law(C, m))
    assert n == pytest.approx(expected, rel=1e-12)
    assert n == pytest.approx(2.346027657952109e12, rel=1e-9)
    assert n > 0.0


# --------------------------------------------------------------------------
# U. analytical life, hand calculation, m == 2 (logarithmic branch)
# --------------------------------------------------------------------------
def test_analytical_life_hand_calculation_m_equals_two():
    """m = 2, C = 1e-22, Y = 1, delta_sigma = 100 MPa, a: 1 mm -> 10 mm.

    (Y*ds*sqrt(pi))**2 = (1.772453850905516e8)**2 = 3.141592653589793e16
    C * that = 1e-22 * 3.141592653589793e16 = 3.141592653589793e-6
    ln(af/a0) = ln(10) = 2.302585092994046
    N = 2.302585092994046 / 3.141592653589793e-6 = 732935.5988794279 cycles
    """
    n = analytical_cycles_to_crack_length(
        1.0e-3, 10.0e-3, CYCLE_100, UNIT, _law(1.0e-22, 2.0)
    )
    assert n == pytest.approx(732935.5988794279, rel=1e-12)


def test_m_equals_two_branch_is_the_limit_of_the_power_branch():
    """The log branch must agree with the power branch as m -> 2."""
    a0, af, C = 1.0e-3, 10.0e-3, 1.0e-22
    at_two = analytical_cycles_to_crack_length(a0, af, CYCLE_100, UNIT, _law(C, 2.0))
    for m in (2.0 + 1.0e-6, 2.0 - 1.0e-6):
        # C is re-scaled so both cases share the same rate at a = 1 m.
        near = analytical_cycles_to_crack_length(
            a0, af, CYCLE_100, UNIT, _law(C, m)
        )
        assert near == pytest.approx(at_two, rel=1e-4)


# --------------------------------------------------------------------------
# V. analytical life is positive for m > 2 (the double-negative sign case)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("m", [2.5, 3.0, 3.5, 4.0, 5.0])
def test_analytical_life_positive_for_m_above_two(m):
    """For m > 2, p < 0 so numerator and denominator are both negative."""
    a0, af = 1.0e-3, 10.0e-3
    p = 1.0 - 0.5 * m
    assert p < 0.0
    assert af**p - a0**p < 0.0  # numerator negative
    n = analytical_cycles_to_crack_length(a0, af, CYCLE_100, UNIT, _law(1.0e-30, m))
    assert n > 0.0
    assert math.isfinite(n)


def test_analytical_life_hand_calculation_m_equals_four():
    n = analytical_cycles_to_crack_length(
        1.0e-3, 10.0e-3, CYCLE_100, UNIT, _law(1.0e-36, 4.0)
    )
    assert n == pytest.approx(911890.6527810405, rel=1e-12)


# --------------------------------------------------------------------------
# W. numerical life positive
# --------------------------------------------------------------------------
@pytest.mark.parametrize("m", [1.5, 2.0, 2.5, 3.0, 4.0])
def test_numerical_life_is_positive_and_finite(m):
    result = cycles_to_crack_length(
        1.0e-3, 10.0e-3, CYCLE_100, UNIT, _law(1.0e-30, m)
    )
    assert result.predicted_cycles > 0.0
    assert math.isfinite(result.predicted_cycles)


# --------------------------------------------------------------------------
# X. numerical vs analytical agreement across the parameter space
# --------------------------------------------------------------------------
@pytest.mark.parametrize("m", [1.0, 1.5, 1.9, 2.0, 2.1, 2.5, 3.0, 3.5, 4.0, 5.0])
@pytest.mark.parametrize("a_range", [(1.0e-3, 10.0e-3), (0.5e-3, 3.0e-3), (2.0e-3, 50.0e-3)])
@pytest.mark.parametrize("sigma_max", [40.0e6, 120.0e6, 300.0e6])
@pytest.mark.parametrize("geometry_factor", [0.8, 1.0, 1.12, 2.0])
def test_numerical_matches_analytical(m, a_range, sigma_max, geometry_factor):
    a0, af = a_range
    cycle = StressCycle(sigma_max=sigma_max, sigma_min=0.2 * sigma_max)
    geometry = ThroughCrackGeometry(geometry_factor)
    law = _law(1.0e-30, m)
    numerical = cycles_to_crack_length(
        a0, af, cycle, geometry, law, intervals=2000
    ).predicted_cycles
    analytical = analytical_cycles_to_crack_length(a0, af, cycle, geometry, law)
    assert numerical == pytest.approx(analytical, rel=1e-8)


# --------------------------------------------------------------------------
# Y. convergence under interval refinement
# --------------------------------------------------------------------------
def test_numerical_result_converges_under_refinement():
    a0, af, law = 1.0e-3, 10.0e-3, _law(1.0e-29, 3.0)
    analytical = analytical_cycles_to_crack_length(a0, af, CYCLE_100, UNIT, law)
    errors = []
    for intervals in (4, 8, 16, 32, 64):
        n = cycles_to_crack_length(
            a0, af, CYCLE_100, UNIT, law, intervals=intervals
        ).predicted_cycles
        errors.append(abs(n - analytical) / analytical)
    assert all(
        finer < coarser for coarser, finer in zip(errors, errors[1:])
    ), errors
    assert errors[-1] < errors[0] / 100.0


def test_simpson_error_falls_at_the_expected_order():
    """Composite Simpson is O(h**4): halving h must cut the error ~16-fold.

    The integrand ``a ** (-m/2)`` is steep near ``a0``, so the asymptotic order
    is only reached once the grid resolves that region; the check is therefore
    made on a refined grid rather than on a handful of intervals.
    """
    a0, af, law = 1.0e-3, 10.0e-3, _law(1.0e-29, 3.0)
    analytical = analytical_cycles_to_crack_length(a0, af, CYCLE_100, UNIT, law)

    def err(n_intervals: int) -> float:
        n = cycles_to_crack_length(
            a0, af, CYCLE_100, UNIT, law, intervals=n_intervals
        ).predicted_cycles
        return abs(n - analytical)

    for coarse in (256, 512):
        ratio = err(coarse) / err(2 * coarse)
        assert 14.0 < ratio < 18.0, (coarse, ratio)


# --------------------------------------------------------------------------
# Z. determinism
# --------------------------------------------------------------------------
def test_repeated_integration_is_bit_for_bit_identical():
    args = (1.0e-3, 10.0e-3, CYCLE_100, UNIT, _law(1.0e-29, 3.0))
    first = cycles_to_crack_length(*args, intervals=512)
    for _ in range(5):
        assert cycles_to_crack_length(*args, intervals=512) == first


# --------------------------------------------------------------------------
# AA-AE. monotonic trends
# --------------------------------------------------------------------------
def _life(a0=1.0e-3, af=10.0e-3, cycle=CYCLE_100, geometry=UNIT, law=None) -> float:
    law = law or _law(1.0e-29, 3.0)
    return cycles_to_crack_length(a0, af, cycle, geometry, law).predicted_cycles


# AA. life decreases when stress range increases
def test_life_decreases_with_increasing_stress_range():
    lives = [
        _life(cycle=StressCycle(sigma_max=s, sigma_min=0.0))
        for s in (50.0e6, 100.0e6, 150.0e6, 200.0e6)
    ]
    assert all(later < earlier for earlier, later in zip(lives, lives[1:]))


def test_life_scales_as_stress_range_to_the_minus_m():
    m = 3.0
    law = _law(1.0e-29, m)
    base = _life(cycle=StressCycle(sigma_max=100.0e6, sigma_min=0.0), law=law)
    doubled = _life(cycle=StressCycle(sigma_max=200.0e6, sigma_min=0.0), law=law)
    assert doubled == pytest.approx(base * 2.0**-m, rel=1e-9)


# AB. life decreases when Y increases
def test_life_decreases_with_increasing_geometry_factor():
    lives = [
        _life(geometry=ThroughCrackGeometry(y)) for y in (0.8, 1.0, 1.2, 1.5, 2.0)
    ]
    assert all(later < earlier for earlier, later in zip(lives, lives[1:]))


def test_life_scales_as_geometry_factor_to_the_minus_m():
    m = 3.0
    law = _law(1.0e-29, m)
    base = _life(geometry=ThroughCrackGeometry(1.0), law=law)
    doubled = _life(geometry=ThroughCrackGeometry(2.0), law=law)
    assert doubled == pytest.approx(base * 2.0**-m, rel=1e-9)


# AC. life decreases when C increases
def test_life_decreases_with_increasing_paris_coefficient():
    lives = [_life(law=_law(c, 3.0)) for c in (1.0e-30, 1.0e-29, 1.0e-28)]
    assert all(later < earlier for earlier, later in zip(lives, lives[1:]))


def test_life_is_inversely_proportional_to_paris_coefficient():
    assert _life(law=_law(1.0e-28, 3.0)) == pytest.approx(
        _life(law=_law(1.0e-29, 3.0)) / 10.0, rel=1e-9
    )


# AD. life decreases as the initial crack gets larger, for a fixed target
def test_life_decreases_with_larger_initial_crack():
    lives = [_life(a0=a0) for a0 in (0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3)]
    assert all(later < earlier for earlier, later in zip(lives, lives[1:]))


# AE. life increases as the target gets farther away
def test_life_increases_with_farther_target():
    lives = [_life(af=af) for af in (2.0e-3, 5.0e-3, 10.0e-3, 25.0e-3)]
    assert all(later > earlier for earlier, later in zip(lives, lives[1:]))


def test_life_is_additive_over_a_split_crack_interval():
    law = _law(1.0e-29, 3.0)
    whole = _life(a0=1.0e-3, af=10.0e-3, law=law)
    first = _life(a0=1.0e-3, af=4.0e-3, law=law)
    second = _life(a0=4.0e-3, af=10.0e-3, law=law)
    assert first + second == pytest.approx(whole, rel=1e-6)


# --------------------------------------------------------------------------
# AF. zero stress range
# --------------------------------------------------------------------------
def test_zero_stress_range_gives_infinite_life_numerically():
    cycle = StressCycle(sigma_max=80.0e6, sigma_min=80.0e6)
    result = cycles_to_crack_length(
        1.0e-3, 10.0e-3, cycle, UNIT, _law(1.0e-29, 3.0)
    )
    assert result.predicted_cycles == math.inf
    assert result.initial_delta_k == 0.0
    assert result.final_delta_k == 0.0
    assert result.initial_growth_rate == 0.0
    assert result.final_growth_rate == 0.0


def test_zero_stress_range_gives_infinite_life_analytically():
    cycle = StressCycle(sigma_max=0.0, sigma_min=0.0)
    assert (
        analytical_cycles_to_crack_length(
            1.0e-3, 10.0e-3, cycle, UNIT, _law(1.0e-29, 3.0)
        )
        == math.inf
    )


# --------------------------------------------------------------------------
# AG. invalid crack-length ordering
# --------------------------------------------------------------------------
@pytest.mark.parametrize("a0, af", [(10.0e-3, 1.0e-3), (5.0e-3, 5.0e-3)])
@pytest.mark.parametrize(
    "func", [cycles_to_crack_length, analytical_cycles_to_crack_length]
)
def test_final_must_exceed_initial(func, a0, af):
    with pytest.raises(ValueError, match="a_final must be strictly greater"):
        func(a0, af, CYCLE_100, UNIT, _law(1.0e-29, 3.0))


@pytest.mark.parametrize("bad", [0.0, -1.0e-3])
@pytest.mark.parametrize(
    "func", [cycles_to_crack_length, analytical_cycles_to_crack_length]
)
def test_non_positive_crack_lengths_rejected(func, bad):
    with pytest.raises(ValueError, match="strictly positive"):
        func(bad, 10.0e-3, CYCLE_100, UNIT, _law(1.0e-29, 3.0))


# --------------------------------------------------------------------------
# AH. non-finite input rejection
# --------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
@pytest.mark.parametrize(
    "func", [cycles_to_crack_length, analytical_cycles_to_crack_length]
)
def test_non_finite_crack_lengths_rejected(func, bad):
    with pytest.raises(ValueError, match="finite"):
        func(bad, 10.0e-3, CYCLE_100, UNIT, _law(1.0e-29, 3.0))
    with pytest.raises(ValueError, match="finite"):
        func(1.0e-3, bad, CYCLE_100, UNIT, _law(1.0e-29, 3.0))


@pytest.mark.parametrize(
    "func", [cycles_to_crack_length, analytical_cycles_to_crack_length]
)
def test_wrong_argument_types_rejected(func):
    law = _law(1.0e-29, 3.0)
    with pytest.raises(TypeError):
        func(1.0e-3, 10.0e-3, "cycle", UNIT, law)
    with pytest.raises(TypeError):
        func(1.0e-3, 10.0e-3, CYCLE_100, "geometry", law)
    with pytest.raises(TypeError):
        func(1.0e-3, 10.0e-3, CYCLE_100, UNIT, "law")


# --------------------------------------------------------------------------
# AI. even-interval requirement (composite Simpson)
# --------------------------------------------------------------------------
@pytest.mark.parametrize("odd", [3, 5, 101, 999])
def test_odd_interval_count_rejected(odd):
    with pytest.raises(ValueError, match="even"):
        cycles_to_crack_length(
            1.0e-3, 10.0e-3, CYCLE_100, UNIT, _law(1.0e-29, 3.0), intervals=odd
        )


@pytest.mark.parametrize("bad", [0, 1, -2])
def test_too_few_intervals_rejected(bad):
    with pytest.raises(ValueError, match="at least 2"):
        cycles_to_crack_length(
            1.0e-3, 10.0e-3, CYCLE_100, UNIT, _law(1.0e-29, 3.0), intervals=bad
        )


@pytest.mark.parametrize("bad", [2.0, "8", None, True])
def test_non_integer_interval_count_rejected(bad):
    with pytest.raises(TypeError, match="intervals must be an int"):
        cycles_to_crack_length(
            1.0e-3, 10.0e-3, CYCLE_100, UNIT, _law(1.0e-29, 3.0), intervals=bad
        )


def test_minimum_even_interval_count_is_accepted():
    result = cycles_to_crack_length(
        1.0e-3, 10.0e-3, CYCLE_100, UNIT, _law(1.0e-29, 3.0), intervals=2
    )
    assert result.integration_intervals == 2
    assert result.predicted_cycles > 0.0


# --------------------------------------------------------------------------
# AJ. result endpoint fields verified independently
# --------------------------------------------------------------------------
def test_result_endpoint_fields_are_independently_correct():
    a0, af = 1.0e-3, 10.0e-3
    law = _law(1.0e-29, 3.0)
    result = cycles_to_crack_length(a0, af, CYCLE_100, UNIT, law, intervals=250)

    assert isinstance(result, CrackGrowthResult)
    assert result.initial_crack_length == a0
    assert result.final_crack_length == af
    assert result.integration_intervals == 250

    # delta_K endpoints against literal hand values, in Pa*sqrt(m)
    assert result.initial_delta_k == pytest.approx(5.604991216397929e6, rel=1e-12)
    assert result.final_delta_k == pytest.approx(1.7724538509055160e7, rel=1e-12)
    # and against the stress-intensity module
    assert result.initial_delta_k == pytest.approx(
        delta_stress_intensity(a0, CYCLE_100, UNIT)
    )
    assert result.final_delta_k == pytest.approx(
        delta_stress_intensity(af, CYCLE_100, UNIT)
    )

    # growth-rate endpoints against literal hand values, in m/cycle
    assert result.initial_growth_rate == pytest.approx(1.7608599228871055e-09, rel=1e-10)
    assert result.final_growth_rate == pytest.approx(5.568327996831707e-08, rel=1e-10)
    assert result.initial_growth_rate == pytest.approx(
        crack_growth_rate(a0, CYCLE_100, UNIT, law)
    )
    assert result.final_growth_rate == pytest.approx(
        crack_growth_rate(af, CYCLE_100, UNIT, law)
    )

    # the crack grows faster at the longer crack length
    assert result.final_growth_rate > result.initial_growth_rate
    assert result.final_delta_k > result.initial_delta_k


def test_result_is_immutable():
    result = cycles_to_crack_length(
        1.0e-3, 10.0e-3, CYCLE_100, UNIT, _law(1.0e-29, 3.0)
    )
    with pytest.raises(Exception):
        result.predicted_cycles = 1.0
