"""A-T, AR, AS. Threshold representation, hard-cutoff rate, threshold boundary."""

from __future__ import annotations

import dataclasses
import math

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_GEOMETRY,
    CANONICAL_GROWTH_THRESHOLD,
    CANONICAL_PARIS_LAW,
    ILLUSTRATIVE_THRESHOLD_DISCLAIMER,
    CrackGrowthThreshold,
    FiniteWidthCenterCrack,
    StressCycle,
    ThresholdBoundaryStatus,
    ThroughCrackGeometry,
    crack_growth_rate,
    crack_growth_rate_for_geometry,
    delta_k_for_geometry,
    finite_width_threshold_crack_length,
    mpa_sqrt_m_to_pa_sqrt_m,
    threshold_crack_length_constant_y,
    thresholded_crack_growth_rate,
)

PANEL = CANONICAL_FINITE_WIDTH_GEOMETRY
UNIT = ThroughCrackGeometry(geometry_factor=1.0)
LAW = CANONICAL_PARIS_LAW
T4 = CANONICAL_GROWTH_THRESHOLD
A0 = 1.0e-3


def _threshold(mpa: float) -> CrackGrowthThreshold:
    return CrackGrowthThreshold(
        name="t",
        delta_k_threshold=mpa_sqrt_m_to_pa_sqrt_m(mpa),
        source_note="ILLUSTRATIVE (test fixture)",
    )


# --------------------------------------------------------------------------
# A. dataclass validation
# --------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [0.0, -1.0e6])
def test_non_positive_threshold_rejected(bad):
    with pytest.raises(ValueError, match="strictly positive"):
        CrackGrowthThreshold(
            name="x", delta_k_threshold=bad, source_note="ILLUSTRATIVE"
        )


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_non_finite_threshold_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        CrackGrowthThreshold(
            name="x", delta_k_threshold=bad, source_note="ILLUSTRATIVE"
        )


def test_threshold_is_immutable():
    with pytest.raises(dataclasses.FrozenInstanceError):
        T4.delta_k_threshold = 1.0


def test_empty_name_rejected():
    with pytest.raises(ValueError, match="name"):
        CrackGrowthThreshold(
            name="  ", delta_k_threshold=4.0e6, source_note="ILLUSTRATIVE"
        )


def test_notes_must_be_a_string():
    with pytest.raises(TypeError, match="notes"):
        CrackGrowthThreshold(
            name="x",
            delta_k_threshold=4.0e6,
            source_note="ILLUSTRATIVE",
            notes=object(),
        )


# --------------------------------------------------------------------------
# B. provenance
# --------------------------------------------------------------------------
def test_empty_provenance_rejected():
    with pytest.raises(ValueError, match="source_note"):
        CrackGrowthThreshold(name="x", delta_k_threshold=4.0e6, source_note="")


def test_canonical_threshold_is_flagged_illustrative():
    assert ILLUSTRATIVE_THRESHOLD_DISCLAIMER in T4.source_note


# --------------------------------------------------------------------------
# C. unit conversion (reuses the toughness helper, not a duplicate)
# --------------------------------------------------------------------------
def test_canonical_threshold_value_and_units():
    assert T4.delta_k_threshold == pytest.approx(4.0e6, rel=1e-12)
    assert T4.delta_k_threshold_in_mpa_sqrt_m == pytest.approx(4.0, rel=1e-12)


@pytest.mark.parametrize("mpa, pa", [(2.0, 2.0e6), (4.0, 4.0e6), (0.5, 5.0e5)])
def test_threshold_unit_conversion(mpa, pa):
    assert _threshold(mpa).delta_k_threshold == pytest.approx(pa, rel=1e-12)
    assert _threshold(mpa).delta_k_threshold_in_mpa_sqrt_m == pytest.approx(
        mpa, rel=1e-12
    )


# --------------------------------------------------------------------------
# D/E/F. the hard cutoff
# --------------------------------------------------------------------------
def test_rate_below_threshold_is_exactly_zero():
    point = thresholded_crack_growth_rate(0.25e-3, CANONICAL_CYCLE, PANEL, LAW, T4)
    assert point.delta_k < T4.delta_k_threshold
    assert point.active_growth is False
    assert point.growth_rate == 0.0
    assert point.threshold_ratio < 1.0
    assert point.threshold_excess_ratio < 0.0


def test_rate_exactly_at_threshold_is_zero():
    """Equality belongs to the NO-GROWTH side.

    The crack length is constructed algebraically so delta_K lands on the
    threshold as closely as floating point allows; the assertion is then made
    on the normalised ratio rather than on a raw difference near zero.
    """
    cycle = StressCycle(sigma_max=120.0e6, sigma_min=20.0e6)
    a_th = threshold_crack_length_constant_y(cycle, UNIT, T4)
    delta_k = delta_k_for_geometry(a_th, cycle, UNIT)
    assert delta_k / T4.delta_k_threshold == pytest.approx(1.0, rel=1e-12)

    # exact equality case, built directly from the threshold value
    assert T4.is_active(T4.delta_k_threshold) is False
    assert T4.threshold_ratio(T4.delta_k_threshold) == 1.0
    assert T4.threshold_excess_ratio(T4.delta_k_threshold) == 0.0


def test_rate_above_threshold_is_the_unmodified_paris_rate():
    """No (delta_K - delta_K_th) subtraction and no (delta_K**m - delta_K_th**m)."""
    for a in (1.0e-3, 2.0e-3, 5.0e-3, 10.0e-3):
        point = thresholded_crack_growth_rate(a, CANONICAL_CYCLE, PANEL, LAW, T4)
        assert point.active_growth is True
        unmodified = crack_growth_rate_for_geometry(
            a, CANONICAL_CYCLE, PANEL, LAW
        )
        assert point.growth_rate == unmodified


def test_rate_is_not_a_reduced_range_law():
    """Guard against the two common alternative near-threshold formulations."""
    a = 2.0e-3
    point = thresholded_crack_growth_rate(a, CANONICAL_CYCLE, PANEL, LAW, T4)
    reduced_range = LAW.C * (point.delta_k - T4.delta_k_threshold) ** LAW.m
    reduced_power = LAW.C * (
        point.delta_k**LAW.m - T4.delta_k_threshold**LAW.m
    )
    assert point.growth_rate != pytest.approx(reduced_range, rel=1e-6)
    assert point.growth_rate != pytest.approx(reduced_power, rel=1e-6)
    assert point.growth_rate > reduced_range
    assert point.growth_rate > reduced_power


# --------------------------------------------------------------------------
# G. the original no-threshold rate is untouched
# --------------------------------------------------------------------------
def test_original_paris_rate_applies_no_threshold():
    """A crack far below threshold still grows under the ORIGINAL API."""
    a = 0.1e-3
    delta_k = delta_k_for_geometry(a, CANONICAL_CYCLE, PANEL)
    assert delta_k < T4.delta_k_threshold
    assert crack_growth_rate_for_geometry(a, CANONICAL_CYCLE, PANEL, LAW) > 0.0
    assert crack_growth_rate(a, CANONICAL_CYCLE, UNIT, LAW) > 0.0
    # ... while the threshold-aware path reports zero
    assert (
        thresholded_crack_growth_rate(a, CANONICAL_CYCLE, PANEL, LAW, T4).growth_rate
        == 0.0
    )


# --------------------------------------------------------------------------
# H/I. threshold ratio and active flag, hand calculations
# --------------------------------------------------------------------------
def test_threshold_ratio_hand_calculation():
    # delta_K(a0) = 5606374.590785748, delta_K_th = 4e6
    # ratio = 1.401593647696437
    point = thresholded_crack_growth_rate(A0, CANONICAL_CYCLE, PANEL, LAW, T4)
    assert point.delta_k == pytest.approx(5606374.590785748, rel=1e-12)
    assert point.threshold_ratio == pytest.approx(1.401593647696437, rel=1e-12)
    assert point.threshold_excess_ratio == pytest.approx(
        0.401593647696437, rel=1e-12
    )
    assert point.active_growth is True


def test_growth_point_fields_are_consistent():
    point = thresholded_crack_growth_rate(2.0e-3, CANONICAL_CYCLE, PANEL, LAW, T4)
    assert point.crack_length == 2.0e-3
    assert point.geometry_factor == pytest.approx(
        PANEL.geometry_factor_at(2.0e-3), rel=1e-12
    )
    assert point.delta_k_threshold == T4.delta_k_threshold
    assert point.threshold_excess_ratio == pytest.approx(
        point.threshold_ratio - 1.0, rel=1e-12
    )
    assert point.active_growth == (point.threshold_ratio > 1.0)


def test_positive_excess_ratio_means_growing_not_safe():
    """Sign-convention lock: positive excess ratio is the ADVERSE case."""
    growing = thresholded_crack_growth_rate(
        4.0e-3, CANONICAL_CYCLE, PANEL, LAW, T4
    )
    arrested = thresholded_crack_growth_rate(
        0.25e-3, CANONICAL_CYCLE, PANEL, LAW, T4
    )
    assert growing.threshold_excess_ratio > 0.0 and growing.growth_rate > 0.0
    assert arrested.threshold_excess_ratio < 0.0 and arrested.growth_rate == 0.0


# --------------------------------------------------------------------------
# J/K/L/M. constant-Y threshold crack size
# --------------------------------------------------------------------------
def test_constant_y_threshold_crack_size_hand_calculation():
    # a_th = (1/pi) * (4e6 / (1.0 * 100e6))**2 = 0.0016/pi = 5.092958178940651e-4
    assert threshold_crack_length_constant_y(
        CANONICAL_CYCLE, UNIT, T4
    ) == pytest.approx(5.092958178940651e-4, rel=1e-12)


@pytest.mark.parametrize("factor", [0.5, 2.0, 3.0])
def test_threshold_size_scales_as_threshold_squared(factor):
    base = threshold_crack_length_constant_y(CANONICAL_CYCLE, UNIT, T4)
    scaled = threshold_crack_length_constant_y(
        CANONICAL_CYCLE, UNIT, _threshold(4.0 * factor)
    )
    assert scaled == pytest.approx(base * factor**2, rel=1e-12)


@pytest.mark.parametrize("factor", [0.5, 1.5, 2.0])
def test_threshold_size_scales_as_stress_range_to_the_minus_two(factor):
    base = threshold_crack_length_constant_y(
        StressCycle(sigma_max=100.0e6, sigma_min=0.0), UNIT, T4
    )
    scaled = threshold_crack_length_constant_y(
        StressCycle(sigma_max=100.0e6 * factor, sigma_min=0.0), UNIT, T4
    )
    assert scaled == pytest.approx(base * factor**-2, rel=1e-12)


@pytest.mark.parametrize("factor", [0.8, 1.2, 2.0])
def test_threshold_size_scales_as_y_to_the_minus_two(factor):
    base = threshold_crack_length_constant_y(CANONICAL_CYCLE, UNIT, T4)
    scaled = threshold_crack_length_constant_y(
        CANONICAL_CYCLE, ThroughCrackGeometry(factor), T4
    )
    assert scaled == pytest.approx(base * factor**-2, rel=1e-12)


def test_constant_y_threshold_size_infinite_for_zero_stress_range():
    cycle = StressCycle(sigma_max=80.0e6, sigma_min=80.0e6)
    assert threshold_crack_length_constant_y(cycle, UNIT, T4) == math.inf


# --------------------------------------------------------------------------
# N/O/P/Q. finite-width threshold root
# --------------------------------------------------------------------------
def test_finite_width_threshold_root_canonical_value():
    result = finite_width_threshold_crack_length(CANONICAL_CYCLE, PANEL, T4)
    assert result.status is ThresholdBoundaryStatus.FINITE_THRESHOLD_SIZE
    assert result.has_finite_boundary
    assert result.threshold_crack_length == pytest.approx(
        0.0005092306461735948, rel=1e-9
    )


def test_delta_k_at_threshold_root_equals_the_threshold():
    result = finite_width_threshold_crack_length(CANONICAL_CYCLE, PANEL, T4)
    # normalised comparison: the raw quantities are ~4e6, the difference ~1e-4
    assert result.delta_k_at_threshold / T4.delta_k_threshold == pytest.approx(
        1.0, rel=1e-9
    )
    assert delta_k_for_geometry(
        result.threshold_crack_length, CANONICAL_CYCLE, PANEL
    ) / T4.delta_k_threshold == pytest.approx(1.0, rel=1e-9)


def test_just_below_threshold_root_is_arrested_and_just_above_is_active():
    a_th = finite_width_threshold_crack_length(
        CANONICAL_CYCLE, PANEL, T4
    ).threshold_crack_length
    assert delta_k_for_geometry(a_th * 0.999, CANONICAL_CYCLE, PANEL) < (
        T4.delta_k_threshold
    )
    assert delta_k_for_geometry(a_th * 1.001, CANONICAL_CYCLE, PANEL) > (
        T4.delta_k_threshold
    )
    assert not thresholded_crack_growth_rate(
        a_th * 0.999, CANONICAL_CYCLE, PANEL, LAW, T4
    ).active_growth
    assert thresholded_crack_growth_rate(
        a_th * 1.001, CANONICAL_CYCLE, PANEL, LAW, T4
    ).active_growth


def test_threshold_result_diagnostics_are_consistent():
    result = finite_width_threshold_crack_length(CANONICAL_CYCLE, PANEL, T4)
    a_th = result.threshold_crack_length
    assert result.plate_width == PANEL.plate_width
    assert result.total_crack_length_at_threshold == pytest.approx(2.0 * a_th)
    assert result.residual_ligament == pytest.approx(
        PANEL.plate_width - 2.0 * a_th
    )
    assert result.ligament_fraction == pytest.approx(
        1.0 - 2.0 * a_th / PANEL.plate_width
    )
    assert 0.0 < result.ligament_fraction < 1.0
    assert result.geometry_factor_at_threshold == pytest.approx(
        PANEL.geometry_factor_at(a_th), rel=1e-12
    )
    assert result.stress_range == pytest.approx(100.0e6)
    assert result.iterations > 0


# --------------------------------------------------------------------------
# R/S/T. solver determinism and bounds
# --------------------------------------------------------------------------
def test_threshold_solver_is_deterministic():
    first = finite_width_threshold_crack_length(CANONICAL_CYCLE, PANEL, T4)
    for _ in range(5):
        assert finite_width_threshold_crack_length(CANONICAL_CYCLE, PANEL, T4) == first


def test_threshold_root_stays_below_half_width():
    for mpa in (2.0, 4.0, 10.0, 20.0, 40.0):
        result = finite_width_threshold_crack_length(
            CANONICAL_CYCLE, PANEL, _threshold(mpa)
        )
        if result.has_finite_boundary:
            assert 0.0 < result.threshold_crack_length < (
                PANEL.max_admissible_crack_length
            )
            assert math.isfinite(result.geometry_factor_at_threshold)


def test_threshold_root_never_evaluates_at_half_width():
    """Y diverges at W/2 and the panel rejects a >= W/2, so a probe there would
    raise rather than return."""
    result = finite_width_threshold_crack_length(
        CANONICAL_CYCLE, PANEL, T4, tolerance=1.0e-15, max_iterations=400
    )
    assert result.has_finite_boundary
    assert result.threshold_crack_length < PANEL.max_admissible_crack_length


def test_threshold_solver_respects_tolerance():
    coarse = finite_width_threshold_crack_length(
        CANONICAL_CYCLE, PANEL, T4, tolerance=1.0e-6
    )
    fine = finite_width_threshold_crack_length(
        CANONICAL_CYCLE, PANEL, T4, tolerance=1.0e-14
    )
    assert fine.iterations > coarse.iterations
    assert fine.threshold_crack_length == pytest.approx(
        coarse.threshold_crack_length, abs=1.0e-6
    )


def test_threshold_solver_reports_zero_stress_range():
    cycle = StressCycle(sigma_max=80.0e6, sigma_min=80.0e6)
    result = finite_width_threshold_crack_length(cycle, PANEL, T4)
    assert result.status is ThresholdBoundaryStatus.ZERO_STRESS_RANGE
    assert result.threshold_crack_length == math.inf
    assert result.iterations == 0


def test_threshold_solver_reports_already_above_at_lower_bound():
    tiny = _threshold(1.0e-9)
    result = finite_width_threshold_crack_length(
        CANONICAL_CYCLE, PANEL, tiny, lower_bound=1.0e-3
    )
    assert result.status is ThresholdBoundaryStatus.ALREADY_ABOVE_AT_LOWER_BOUND
    assert result.threshold_crack_length == math.inf


def test_threshold_solver_reports_no_root_for_an_unreachable_threshold():
    """delta_K stays below the threshold everywhere in the admissible interval.

    Because Y diverges as a -> W/2, delta_K at the upper bracket is enormous
    (~3.2e13 Pa*sqrt(m) here), so this branch is only reachable with a
    physically absurd threshold. It is implemented and tested honestly rather
    than assumed impossible -- the same treatment the Milestone 3 fracture
    solver gives its equivalent status.
    """
    result = finite_width_threshold_crack_length(
        CANONICAL_CYCLE, PANEL, _threshold(1.0e8)
    )
    assert result.status is ThresholdBoundaryStatus.NO_ROOT_IN_BRACKET
    assert result.threshold_crack_length == math.inf
    assert result.iterations == 0


def test_threshold_solver_rejects_a_lower_bound_above_the_bracket():
    with pytest.raises(ValueError, match="not below the admissible upper bracket"):
        finite_width_threshold_crack_length(
            CANONICAL_CYCLE, PANEL, T4, lower_bound=0.06
        )


# --------------------------------------------------------------------------
# AO/AP. constant-Y cross-check of the numerical root
# --------------------------------------------------------------------------
def test_numerical_root_matches_the_exact_constant_y_formula():
    """A very wide panel makes Y ~ 1, so the root must match the closed form.

    This independently verifies the bisection path against exact algebra."""
    exact = threshold_crack_length_constant_y(CANONICAL_CYCLE, UNIT, T4)
    wide = FiniteWidthCenterCrack(plate_width=1000.0)
    numerical = finite_width_threshold_crack_length(
        CANONICAL_CYCLE, wide, T4
    ).threshold_crack_length
    assert numerical == pytest.approx(exact, rel=1e-8)


def test_finite_width_threshold_size_differs_from_the_constant_y_value():
    """At W = 100 mm the correction is small but real: Y > 1 means the threshold
    is reached at a slightly SMALLER crack than the infinite plate predicts."""
    exact = threshold_crack_length_constant_y(CANONICAL_CYCLE, UNIT, T4)
    finite = finite_width_threshold_crack_length(
        CANONICAL_CYCLE, PANEL, T4
    ).threshold_crack_length
    assert finite < exact
    assert abs(finite - exact) / exact == pytest.approx(1.28e-4, rel=5e-2)


# --------------------------------------------------------------------------
# AR/AS. input rejection
# --------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [0.0, -1.0e-3])
def test_non_positive_crack_length_rejected(bad):
    with pytest.raises(ValueError, match="strictly positive"):
        thresholded_crack_growth_rate(bad, CANONICAL_CYCLE, PANEL, LAW, T4)


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_non_finite_crack_length_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        thresholded_crack_growth_rate(bad, CANONICAL_CYCLE, PANEL, LAW, T4)


def test_wrong_types_rejected():
    with pytest.raises(TypeError):
        thresholded_crack_growth_rate(A0, CANONICAL_CYCLE, PANEL, "law", T4)
    with pytest.raises(TypeError):
        thresholded_crack_growth_rate(A0, CANONICAL_CYCLE, PANEL, LAW, "threshold")
    with pytest.raises(TypeError):
        finite_width_threshold_crack_length(CANONICAL_CYCLE, UNIT, T4)
    with pytest.raises(TypeError):
        threshold_crack_length_constant_y(CANONICAL_CYCLE, PANEL, T4)
    with pytest.raises(TypeError):
        finite_width_threshold_crack_length("cycle", PANEL, T4)


# --------------------------------------------------------------------------
# Monotonicity of delta_K: once active, always active (item 16)
# --------------------------------------------------------------------------
def test_delta_k_increases_monotonically_with_crack_length():
    """For this geometry both sqrt(pi a) and Y(a) increase with a, so a crack
    that is above threshold stays above it as it grows."""
    lengths = [0.1e-3, 0.5e-3, 1.0e-3, 2.0e-3, 5.0e-3, 10.0e-3, 20.0e-3, 40.0e-3]
    values = [delta_k_for_geometry(a, CANONICAL_CYCLE, PANEL) for a in lengths]
    assert all(later > earlier for earlier, later in zip(values, values[1:]))


def test_growth_stays_active_above_the_threshold_crack_size():
    a_th = finite_width_threshold_crack_length(
        CANONICAL_CYCLE, PANEL, T4
    ).threshold_crack_length
    for a in (a_th * 1.01, a_th * 2, a_th * 10, a_th * 30):
        assert thresholded_crack_growth_rate(
            a, CANONICAL_CYCLE, PANEL, LAW, T4
        ).active_growth


def test_compression_can_still_inflate_the_algebraic_range():
    """No closure model: a compressive sigma_min widens delta_K and can lift a
    crack above threshold that a closure-aware model would not."""
    tension_only = StressCycle(sigma_max=60.0e6, sigma_min=20.0e6)
    with_compression = StressCycle(sigma_max=60.0e6, sigma_min=-60.0e6)
    a = 1.0e-3
    assert not thresholded_crack_growth_rate(
        a, tension_only, PANEL, LAW, T4
    ).active_growth
    assert thresholded_crack_growth_rate(
        a, with_compression, PANEL, LAW, T4
    ).active_growth
