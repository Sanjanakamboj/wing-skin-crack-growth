"""Z-AE, AD. Log-grid integrator: validity, convergence, geometry awareness."""

from __future__ import annotations

import math

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_PARIS_LAW,
    FiniteWidthCenterCrack,
    ParisLaw,
    StressCycle,
    ThroughCrackGeometry,
    analytical_cycles_to_crack_length,
    cycles_to_crack_length,
    integrate_crack_growth_log_grid,
)

UNIT = ThroughCrackGeometry(geometry_factor=1.0)
PANEL = CANONICAL_FINITE_WIDTH_GEOMETRY
LAW = CANONICAL_PARIS_LAW


def _law(C=1.0e-29, m=3.0):
    return ParisLaw(C=C, m=m, name="t", source_note="ILLUSTRATIVE (test)")


def _log(a0, af, cycle=CANONICAL_CYCLE, geometry=UNIT, law=LAW, **kwargs):
    return integrate_crack_growth_log_grid(a0, af, cycle, geometry, law, **kwargs)


# --------------------------------------------------------------------------
# AB. constant-Y agreement with the analytical closed form
# --------------------------------------------------------------------------
@pytest.mark.parametrize("m", [1.5, 2.0, 2.5, 3.0, 4.0])
@pytest.mark.parametrize("a_range", [(1.0e-3, 10.0e-3), (0.5e-3, 3.0e-3)])
def test_log_grid_matches_analytical_for_constant_geometry(m, a_range):
    a0, af = a_range
    law = _law(1.0e-30, m)
    analytical = analytical_cycles_to_crack_length(a0, af, CANONICAL_CYCLE, UNIT, law)
    assert _log(a0, af, law=law, intervals=200).predicted_cycles == pytest.approx(
        analytical, rel=1e-9
    )


def test_log_grid_canonical_constant_y_life():
    analytical = analytical_cycles_to_crack_length(
        1.0e-3, 10.0e-3, CANONICAL_CYCLE, UNIT, LAW
    )
    assert _log(1.0e-3, 10.0e-3).predicted_cycles == pytest.approx(
        analytical, rel=1e-12
    )
    assert analytical == pytest.approx(776634.4444503565, rel=1e-12)


# --------------------------------------------------------------------------
# AC. wide-span convergence, and comparison with the linear grid
# --------------------------------------------------------------------------
@pytest.mark.parametrize("span", [10.0, 100.0, 1000.0])
def test_log_grid_converges_on_wide_spans(span):
    a0 = 1.0e-5
    af = a0 * span
    analytical = analytical_cycles_to_crack_length(a0, af, CANONICAL_CYCLE, UNIT, LAW)
    assert _log(a0, af, intervals=200).predicted_cycles == pytest.approx(
        analytical, rel=1e-9
    )


def test_log_grid_error_falls_at_fourth_order():
    a0, af = 1.0e-5, 1.0e-2
    analytical = analytical_cycles_to_crack_length(a0, af, CANONICAL_CYCLE, UNIT, LAW)

    def err(n):
        return abs(_log(a0, af, intervals=n).predicted_cycles - analytical) / analytical

    for coarse in (20, 40, 80):
        assert 14.0 < err(coarse) / err(2 * coarse) < 18.0, coarse


def test_log_grid_beats_the_linear_grid_on_a_wide_span():
    """Measured, not assumed: on af/a0 = 1000 the log grid is far more accurate
    at equal interval count. The linear grid remains correct, just inefficient."""
    a0, af = 1.0e-5, 1.0e-2
    analytical = analytical_cycles_to_crack_length(a0, af, CANONICAL_CYCLE, UNIT, LAW)
    for n in (100, 500, 1000):
        linear = abs(
            cycles_to_crack_length(
                a0, af, CANONICAL_CYCLE, UNIT, LAW, intervals=n
            ).predicted_cycles
            - analytical
        )
        log_grid = abs(_log(a0, af, intervals=n).predicted_cycles - analytical)
        assert log_grid < linear, n


# --------------------------------------------------------------------------
# AA. even-interval requirement and validation
# --------------------------------------------------------------------------
@pytest.mark.parametrize("odd", [3, 5, 101, 999])
def test_odd_interval_count_rejected(odd):
    with pytest.raises(ValueError, match="even"):
        _log(1.0e-3, 10.0e-3, intervals=odd)


@pytest.mark.parametrize("bad", [0, 1, -2])
def test_too_few_intervals_rejected(bad):
    with pytest.raises(ValueError, match="at least 2"):
        _log(1.0e-3, 10.0e-3, intervals=bad)


@pytest.mark.parametrize("bad", [2.0, "8", None, True])
def test_non_integer_interval_count_rejected(bad):
    with pytest.raises(TypeError, match="intervals must be an int"):
        _log(1.0e-3, 10.0e-3, intervals=bad)


@pytest.mark.parametrize("a0, af", [(10.0e-3, 1.0e-3), (5.0e-3, 5.0e-3)])
def test_final_must_exceed_initial(a0, af):
    with pytest.raises(ValueError, match="a_final must be strictly greater"):
        _log(a0, af)


@pytest.mark.parametrize("bad", [0.0, -1.0e-3, math.nan, math.inf])
def test_invalid_crack_lengths_rejected(bad):
    with pytest.raises(ValueError):
        _log(bad, 10.0e-3)
    with pytest.raises(ValueError):
        _log(1.0e-3, bad)


def test_wrong_types_rejected():
    with pytest.raises(TypeError):
        integrate_crack_growth_log_grid(1.0e-3, 10.0e-3, "cycle", UNIT, LAW)
    with pytest.raises(TypeError):
        integrate_crack_growth_log_grid(1.0e-3, 10.0e-3, CANONICAL_CYCLE, UNIT, "law")


# --------------------------------------------------------------------------
# AN. zero stress range
# --------------------------------------------------------------------------
def test_zero_stress_range_gives_infinite_life():
    cycle = StressCycle(sigma_max=80.0e6, sigma_min=80.0e6)
    result = _log(1.0e-3, 10.0e-3, cycle=cycle)
    assert result.predicted_cycles == math.inf
    assert result.initial_delta_k == 0.0
    assert result.initial_growth_rate == 0.0


def test_zero_stress_range_with_finite_width_also_infinite():
    cycle = StressCycle(sigma_max=80.0e6, sigma_min=80.0e6)
    assert _log(1.0e-3, 10.0e-3, cycle=cycle, geometry=PANEL).predicted_cycles == (
        math.inf
    )


# --------------------------------------------------------------------------
# Determinism and bounds
# --------------------------------------------------------------------------
def test_repeated_integration_is_bit_for_bit_identical():
    first = _log(1.0e-3, 10.0e-3, intervals=512)
    for _ in range(5):
        assert _log(1.0e-3, 10.0e-3, intervals=512) == first


def test_integration_never_leaves_the_interval():
    """Every abscissa exp(u0 + i*h) lies within [a0, af] by construction."""
    a0, af, n = 1.0e-3, 10.0e-3, 100
    u0, h = math.log(a0), (math.log(af) - math.log(a0)) / n
    for i in range(n + 1):
        a = a0 if i == 0 else (af if i == n else math.exp(u0 + i * h))
        assert a0 <= a <= af


# --------------------------------------------------------------------------
# AE/AF. geometry is re-evaluated throughout the integration
# --------------------------------------------------------------------------
def test_finite_width_life_differs_from_infinite_plate():
    a0, af = 1.0e-3, 15.0e-3
    infinite = _log(a0, af, geometry=UNIT).predicted_cycles
    finite = _log(a0, af, geometry=PANEL).predicted_cycles
    assert finite < infinite


def test_freezing_y_at_a0_changes_the_answer():
    """If Y were hoisted out of the integral (frozen at a0) the result would
    differ measurably -- proof that Y(a) is re-evaluated at every abscissa."""
    a0, af = 1.0e-3, 15.0e-3
    true_life = _log(a0, af, geometry=PANEL).predicted_cycles
    frozen = ThroughCrackGeometry(
        geometry_factor=PANEL.geometry_factor_at(a0)
    )
    frozen_life = _log(a0, af, geometry=frozen).predicted_cycles
    assert frozen_life != pytest.approx(true_life, rel=1e-4)
    # freezing Y understates the driving force, so it OVERestimates the life
    assert frozen_life > true_life


def test_endpoint_geometry_factors_are_reported():
    a0, af = 1.0e-3, 15.0e-3
    result = _log(a0, af, geometry=PANEL)
    assert result.initial_geometry_factor == pytest.approx(
        PANEL.geometry_factor_at(a0), rel=1e-12
    )
    assert result.final_geometry_factor == pytest.approx(
        PANEL.geometry_factor_at(af), rel=1e-12
    )
    assert result.final_geometry_factor > result.initial_geometry_factor
    assert result.integration_method == "composite Simpson on a log(a) grid"
    assert result.final_growth_rate > result.initial_growth_rate


def test_wider_panel_converges_to_the_constant_y_result():
    a0, af = 1.0e-3, 10.0e-3
    infinite = _log(a0, af, geometry=UNIT).predicted_cycles
    for width, tolerance in ((1.0, 1e-3), (10.0, 1e-5), (1000.0, 1e-9)):
        panel = FiniteWidthCenterCrack(plate_width=width)
        assert _log(a0, af, geometry=panel).predicted_cycles == pytest.approx(
            infinite, rel=tolerance
        )


# --------------------------------------------------------------------------
# AD. the Milestone 1 linear-grid integrator is untouched
# --------------------------------------------------------------------------
def test_linear_grid_integrator_is_unchanged():
    result = cycles_to_crack_length(
        1.0e-3, 10.0e-3, CANONICAL_CYCLE, UNIT, LAW
    )
    assert result.predicted_cycles == pytest.approx(776634.4447, abs=1e-3)
    assert result.integration_intervals == 1000
    assert not hasattr(result, "integration_method")


def test_the_two_integrators_agree_for_the_canonical_case():
    linear = cycles_to_crack_length(
        1.0e-3, 10.0e-3, CANONICAL_CYCLE, UNIT, LAW
    ).predicted_cycles
    log_grid = _log(1.0e-3, 10.0e-3).predicted_cycles
    assert log_grid == pytest.approx(linear, rel=1e-8)
