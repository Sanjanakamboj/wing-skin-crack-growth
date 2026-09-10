"""A. Geometry validation."""

from __future__ import annotations

import dataclasses
import math

import pytest

from crackgrowth import INFINITE_PLATE_THROUGH_CRACK, ThroughCrackGeometry


def test_geometry_factor_is_stored():
    assert ThroughCrackGeometry(geometry_factor=1.12).geometry_factor == 1.12


def test_canonical_geometry_factor_is_unity():
    assert INFINITE_PLATE_THROUGH_CRACK.geometry_factor == 1.0


@pytest.mark.parametrize("bad", [0.0, -1.0, -0.5])
def test_non_positive_geometry_factor_rejected(bad):
    with pytest.raises(ValueError, match="strictly positive"):
        ThroughCrackGeometry(geometry_factor=bad)


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_non_finite_geometry_factor_rejected(bad):
    with pytest.raises(ValueError, match="finite"):
        ThroughCrackGeometry(geometry_factor=bad)


def test_geometry_is_immutable():
    geometry = ThroughCrackGeometry(geometry_factor=1.0)
    with pytest.raises(dataclasses.FrozenInstanceError):
        geometry.geometry_factor = 2.0


def test_integer_geometry_factor_is_coerced_to_float():
    geometry = ThroughCrackGeometry(geometry_factor=2)
    assert isinstance(geometry.geometry_factor, float)
