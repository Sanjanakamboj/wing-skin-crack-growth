"""Shared fixtures. Expected values in the tests are computed independently."""

from __future__ import annotations

import pytest

from crackgrowth import (
    INFINITE_PLATE_THROUGH_CRACK,
    ParisLaw,
    StressCycle,
    ThroughCrackGeometry,
)


@pytest.fixture
def unit_geometry() -> ThroughCrackGeometry:
    return INFINITE_PLATE_THROUGH_CRACK


@pytest.fixture
def cycle_100mpa_range() -> StressCycle:
    """sigma_max = 120 MPa, sigma_min = 20 MPa -> delta_sigma = 100 MPa."""
    return StressCycle(sigma_max=120.0e6, sigma_min=20.0e6)


@pytest.fixture
def simple_paris() -> ParisLaw:
    """C = 1e-29 (SI basis), m = 3.0 -- the canonical illustrative curve."""
    return ParisLaw(
        C=1.0e-29,
        m=3.0,
        name="test curve",
        source_note="ILLUSTRATIVE PARIS-LAW INPUT - NOT DESIGN ALLOWABLE (test fixture)",
    )
