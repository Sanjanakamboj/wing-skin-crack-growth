"""A-N, AR, AS. Finite-width geometry, conventions, Y(a), K and ΔK."""

from __future__ import annotations

import dataclasses
import math

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_PARIS_LAW,
    CANONICAL_PLATE_WIDTH,
    CrackGeometry,
    FiniteWidthCenterCrack,
    StressCycle,
    ThroughCrackGeometry,
    crack_growth_rate,
    crack_growth_rate_for_geometry,
    delta_k_for_geometry,
    delta_stress_intensity,
    geometry_factor_at,
    ligament_fraction,
    remaining_ligament,
    stress_intensity,
    stress_intensity_for_geometry,
    total_crack_length,
)

W = 0.1
PANEL = FiniteWidthCenterCrack(plate_width=W)
UNIT = ThroughCrackGeometry(geometry_factor=1.0)


# --------------------------------------------------------------------------
# A. geometry validation
# --------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [0.0, -0.1])
def test_non_positive_width_rejected(bad):
    with pytest.raises(ValueError, match="strictly positive"):
        FiniteWidthCenterCrack(plate_width=bad)


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_non_finite_width_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        FiniteWidthCenterCrack(plate_width=bad)


def test_geometry_is_immutable():
    with pytest.raises(dataclasses.FrozenInstanceError):
        PANEL.plate_width = 0.2


def test_canonical_panel_width_is_one_hundred_millimetres():
    assert CANONICAL_PLATE_WIDTH == 100.0e-3
    assert CANONICAL_FINITE_WIDTH_GEOMETRY.plate_width == 100.0e-3


def test_both_geometries_satisfy_the_protocol():
    assert isinstance(PANEL, CrackGeometry)
    assert isinstance(UNIT, CrackGeometry)


def test_constant_geometry_factor_at_is_crack_length_independent():
    """The Milestone 1 geometry gained geometry_factor_at without changing."""
    for a in (1.0e-3, 5.0e-3, 50.0e-3):
        assert UNIT.geometry_factor_at(a) == 1.0
    assert UNIT.geometry_factor == 1.0


# --------------------------------------------------------------------------
# B/C/D/E. half-crack convention and ligament diagnostics
# --------------------------------------------------------------------------
def test_a_is_the_half_crack_length_and_total_is_twice_it():
    a = 20.0e-3
    assert PANEL.total_crack_length(a) == pytest.approx(40.0e-3)
    assert total_crack_length(a) == pytest.approx(40.0e-3)


def test_ligament_is_width_minus_total_crack_length():
    a = 20.0e-3
    # W - 2a = 0.1 - 0.04 = 0.06 m
    assert PANEL.remaining_ligament(a) == pytest.approx(0.06)
    assert remaining_ligament(a, W) == pytest.approx(0.06)


def test_ligament_fraction_definition():
    a = 20.0e-3
    # 1 - 2a/W = 1 - 0.4 = 0.6
    assert PANEL.ligament_fraction(a) == pytest.approx(0.6)
    assert ligament_fraction(a, W) == pytest.approx(0.6)


@pytest.mark.parametrize("a", [1.0e-3, 10.0e-3, 24.9e-3, 49.0e-3])
def test_ligament_fraction_is_strictly_between_zero_and_one(a):
    fraction = PANEL.ligament_fraction(a)
    assert 0.0 < fraction < 1.0


def test_half_crack_convention_is_not_confused_with_total_length():
    """A 20 mm HALF crack fills 40 % of a 100 mm panel, not 20 %."""
    a = 20.0e-3
    assert PANEL.total_crack_length(a) / W == pytest.approx(0.4)
    assert a / W == pytest.approx(0.2)
    # Y depends on a/W, not on 2a/W: mistaking one for the other would give
    # Y(0.4) = 1.7989 instead of the correct Y(0.2) = 1.1118.
    assert PANEL.geometry_factor_at(a) == pytest.approx(1.1117859405028423, rel=1e-12)
    assert PANEL.geometry_factor_at(a) != pytest.approx(1.798907, rel=1e-4)


# --------------------------------------------------------------------------
# AS. admissibility bound
# --------------------------------------------------------------------------
def test_max_admissible_crack_length_is_half_the_width():
    assert PANEL.max_admissible_crack_length == pytest.approx(0.05)


@pytest.mark.parametrize("bad", [0.05, 0.06, 0.1])
def test_crack_at_or_beyond_half_width_rejected(bad):
    with pytest.raises(ValueError, match="not admissible"):
        PANEL.geometry_factor_at(bad)
    with pytest.raises(ValueError, match="not admissible"):
        PANEL.ligament_fraction(bad)
    with pytest.raises(ValueError, match="not admissible"):
        PANEL.total_crack_length(bad)


def test_crack_just_below_half_width_is_accepted():
    a = 0.05 * (1.0 - 1.0e-9)
    assert math.isfinite(PANEL.geometry_factor_at(a))
    assert PANEL.geometry_factor_at(a) > 100.0


# --------------------------------------------------------------------------
# AR. non-finite / non-positive rejection
# --------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [0.0, -1.0e-3])
def test_non_positive_crack_length_rejected(bad):
    with pytest.raises(ValueError, match="strictly positive"):
        PANEL.geometry_factor_at(bad)


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_non_finite_crack_length_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        PANEL.geometry_factor_at(bad)
    with pytest.raises(ValueError, match="finite"):
        stress_intensity_for_geometry(100.0e6, bad, PANEL)


# --------------------------------------------------------------------------
# F. Y hand calculation
# --------------------------------------------------------------------------
@pytest.mark.parametrize(
    "ratio, expected",
    [
        (0.01, 1.0002468111606977),
        (0.05, 1.0062132605904193),
        (0.10, 1.0254083207377767),
        (0.20, 1.1117859405028423),
        (0.30, 1.3043395327536766),
        (0.40, 1.7989074399478673),
    ],
)
def test_geometry_factor_hand_calculation(ratio, expected):
    # Y = 1 / sqrt(cos(pi * a / W)), with a = ratio * W
    assert PANEL.geometry_factor_at(ratio * W) == pytest.approx(expected, rel=1e-12)
    assert PANEL.geometry_factor_at(ratio * W) == pytest.approx(
        1.0 / math.sqrt(math.cos(math.pi * ratio)), rel=1e-14
    )


def test_geometry_factor_equals_sqrt_of_secant():
    """Y = sqrt(sec(pi a / W)) is the same expression as 1/sqrt(cos(...))."""
    for a in (1.0e-3, 10.0e-3, 20.0e-3, 40.0e-3):
        assert PANEL.geometry_factor_at(a) == pytest.approx(
            math.sqrt(1.0 / math.cos(math.pi * a / W)), rel=1e-14
        )


# --------------------------------------------------------------------------
# G/H/I/J. Y behaviour
# --------------------------------------------------------------------------
def test_geometry_factor_tends_to_one_for_small_crack_ratio():
    assert PANEL.geometry_factor_at(1.0e-9) == pytest.approx(1.0, abs=1e-12)
    assert PANEL.geometry_factor_at(1.0e-6) == pytest.approx(1.0, abs=1e-9)


@pytest.mark.parametrize("a", [0.5e-3, 1.0e-3, 5.0e-3, 20.0e-3, 40.0e-3])
def test_geometry_factor_is_greater_than_one(a):
    assert PANEL.geometry_factor_at(a) > 1.0


def test_geometry_factor_is_monotonic_in_crack_ratio():
    factors = [
        PANEL.geometry_factor_at(r * W)
        for r in (0.01, 0.05, 0.10, 0.20, 0.30, 0.40, 0.45)
    ]
    assert all(later > earlier for earlier, later in zip(factors, factors[1:]))


def test_geometry_factor_grows_strongly_near_half_width():
    assert PANEL.geometry_factor_at(0.45 * W) == pytest.approx(2.528330, rel=1e-5)
    assert PANEL.geometry_factor_at(0.49 * W) == pytest.approx(5.642360, rel=1e-5)
    assert PANEL.geometry_factor_at(0.499 * W) > 17.0
    assert PANEL.geometry_factor_at(0.4999 * W) > 56.0


def test_wider_panel_gives_lower_geometry_factor_at_fixed_crack_length():
    a = 10.0e-3
    factors = [
        FiniteWidthCenterCrack(plate_width=w).geometry_factor_at(a)
        for w in (0.05, 0.075, 0.1, 0.2, 0.5, 5.0)
    ]
    assert all(later < earlier for earlier, later in zip(factors, factors[1:]))
    assert factors[-1] == pytest.approx(1.0, abs=1e-5)


# --------------------------------------------------------------------------
# K/L. finite-width K
# --------------------------------------------------------------------------
def test_finite_width_stress_intensity_hand_calculation():
    # K = Y(0.2) * 120e6 * sqrt(pi * 0.02)
    #   = 1.1117859405028423 * 1.2e8 * 0.25066282746310002 = 3.3442008885619722e7
    k = stress_intensity_for_geometry(120.0e6, 20.0e-3, PANEL)
    assert k == pytest.approx(3.3442008885619722e7, rel=1e-12)


def test_finite_width_stress_intensity_exceeds_infinite_plate():
    a = 20.0e-3
    assert stress_intensity_for_geometry(120.0e6, a, PANEL) > stress_intensity(
        120.0e6, a, UNIT
    )
    assert stress_intensity_for_geometry(120.0e6, a, PANEL) == pytest.approx(
        stress_intensity(120.0e6, a, UNIT) * PANEL.geometry_factor_at(a), rel=1e-12
    )


def test_finite_width_stress_intensity_preserves_stress_sign():
    a = 20.0e-3
    tensile = stress_intensity_for_geometry(120.0e6, a, PANEL)
    compressive = stress_intensity_for_geometry(-120.0e6, a, PANEL)
    assert tensile > 0.0 > compressive
    assert compressive == pytest.approx(-tensile)
    assert stress_intensity_for_geometry(0.0, a, PANEL) == 0.0


def test_generic_helper_reproduces_milestone1_for_constant_geometry():
    """The geometry-aware helpers must reduce EXACTLY to the M1 functions."""
    for a in (1.0e-3, 5.0e-3, 10.0e-3):
        assert stress_intensity_for_geometry(120.0e6, a, UNIT) == pytest.approx(
            stress_intensity(120.0e6, a, UNIT), rel=1e-15
        )
        assert delta_k_for_geometry(a, CANONICAL_CYCLE, UNIT) == pytest.approx(
            delta_stress_intensity(a, CANONICAL_CYCLE, UNIT), rel=1e-15
        )
        assert crack_growth_rate_for_geometry(
            a, CANONICAL_CYCLE, UNIT, CANONICAL_PARIS_LAW
        ) == pytest.approx(
            crack_growth_rate(a, CANONICAL_CYCLE, UNIT, CANONICAL_PARIS_LAW),
            rel=1e-15,
        )


# --------------------------------------------------------------------------
# M/N. finite-width ΔK
# --------------------------------------------------------------------------
def test_finite_width_delta_k_hand_calculation():
    # dK = Y(0.2) * 100e6 * sqrt(pi * 0.02) = 2.7868340738016438e7
    dk = delta_k_for_geometry(20.0e-3, CANONICAL_CYCLE, PANEL)
    assert dk == pytest.approx(2.7868340738016438e7, rel=1e-12)


@pytest.mark.parametrize("a", [0.5e-3, 1.0e-3, 5.0e-3, 20.0e-3, 40.0e-3])
def test_finite_width_delta_k_is_never_below_infinite_plate(a):
    assert delta_k_for_geometry(a, CANONICAL_CYCLE, PANEL) >= delta_stress_intensity(
        a, CANONICAL_CYCLE, UNIT
    )


def test_finite_width_delta_k_approaches_infinite_plate_as_width_grows():
    a = 10.0e-3
    infinite = delta_stress_intensity(a, CANONICAL_CYCLE, UNIT)
    for width, tolerance in ((1.0, 1e-3), (10.0, 1e-5), (1000.0, 1e-9)):
        panel = FiniteWidthCenterCrack(plate_width=width)
        assert delta_k_for_geometry(a, CANONICAL_CYCLE, panel) == pytest.approx(
            infinite, rel=tolerance
        )


def test_finite_width_delta_k_still_uses_the_algebraic_range():
    """No crack closure: a compressive sigma_min still widens delta_K."""
    tension_only = StressCycle(sigma_max=100.0e6, sigma_min=0.0)
    with_compression = StressCycle(sigma_max=100.0e6, sigma_min=-100.0e6)
    a = 10.0e-3
    assert delta_k_for_geometry(a, with_compression, PANEL) == pytest.approx(
        2.0 * delta_k_for_geometry(a, tension_only, PANEL)
    )


def test_geometry_factor_at_dispatch_rejects_unsupported_objects():
    with pytest.raises(TypeError, match="geometry_factor_at"):
        geometry_factor_at(object(), 1.0e-3)


def test_growth_rate_for_geometry_rejects_bad_types():
    with pytest.raises(TypeError):
        crack_growth_rate_for_geometry(1.0e-3, CANONICAL_CYCLE, PANEL, "law")
    with pytest.raises(TypeError):
        delta_k_for_geometry(1.0e-3, "cycle", PANEL)
