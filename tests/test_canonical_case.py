"""The canonical wing-skin sanity case, verified end to end."""

from __future__ import annotations

import math

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_GEOMETRY,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CANONICAL_TARGET_CRACK_LENGTH,
    ILLUSTRATIVE_DISCLAIMER,
    analytical_cycles_to_crack_length,
    cycles_to_crack_length,
)


def test_canonical_inputs_are_as_documented():
    assert CANONICAL_GEOMETRY.geometry_factor == 1.0
    assert CANONICAL_CYCLE.sigma_max == 120.0e6
    assert CANONICAL_CYCLE.sigma_min == 20.0e6
    assert CANONICAL_CYCLE.stress_range == pytest.approx(100.0e6)
    assert CANONICAL_CYCLE.stress_ratio == pytest.approx(1.0 / 6.0)
    assert CANONICAL_INITIAL_CRACK_LENGTH == 1.0e-3
    assert CANONICAL_TARGET_CRACK_LENGTH == 10.0e-3
    assert CANONICAL_PARIS_LAW.m == 3.0
    assert CANONICAL_PARIS_LAW.C == pytest.approx(1.0e-29, rel=1e-12)
    assert CANONICAL_PARIS_LAW.c_on_mpa_basis == pytest.approx(1.0e-11, rel=1e-12)
    assert ILLUSTRATIVE_DISCLAIMER in CANONICAL_PARIS_LAW.source_note


def test_canonical_case_life_and_endpoints():
    result = cycles_to_crack_length(
        CANONICAL_INITIAL_CRACK_LENGTH,
        CANONICAL_TARGET_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    )
    analytical = analytical_cycles_to_crack_length(
        CANONICAL_INITIAL_CRACK_LENGTH,
        CANONICAL_TARGET_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    )

    assert analytical == pytest.approx(776634.4444503565, rel=1e-12)
    assert result.predicted_cycles == pytest.approx(analytical, rel=1e-8)

    # Endpoint delta_K in MPa*sqrt(m): a plausible mid-Paris band.
    assert result.initial_delta_k / 1.0e6 == pytest.approx(5.604991, rel=1e-6)
    assert result.final_delta_k / 1.0e6 == pytest.approx(17.724539, rel=1e-6)

    # Endpoint growth rates in m/cycle.
    assert result.initial_growth_rate == pytest.approx(1.76086e-9, rel=1e-5)
    assert result.final_growth_rate == pytest.approx(5.56833e-8, rel=1e-5)


def test_canonical_delta_k_stays_in_a_plausible_paris_band():
    """Sanity on the illustrative curve: well above a typical threshold and
    well below a typical aluminium toughness, so the Paris form is defensible."""
    result = cycles_to_crack_length(
        CANONICAL_INITIAL_CRACK_LENGTH,
        CANONICAL_TARGET_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    )
    assert 3.0e6 < result.initial_delta_k < 30.0e6
    assert 3.0e6 < result.final_delta_k < 30.0e6


def test_canonical_life_is_meaningful_but_inspectable():
    n = analytical_cycles_to_crack_length(
        CANONICAL_INITIAL_CRACK_LENGTH,
        CANONICAL_TARGET_CRACK_LENGTH,
        CANONICAL_CYCLE,
        CANONICAL_GEOMETRY,
        CANONICAL_PARIS_LAW,
    )
    assert 1.0e5 < n < 1.0e7
    assert math.isfinite(n)
