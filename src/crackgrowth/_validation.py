"""Shared input-validation primitives.

These helpers exist only to reject non-physical or non-finite input early and
loudly. They deliberately do NOT repair, clamp, or reinterpret values: every
engineering assumption in this package is made explicit at the call site rather
than hidden inside a convenience helper.
"""

from __future__ import annotations

import math

__all__ = ["require_finite", "require_positive", "require_non_negative"]


def require_finite(value: float, name: str) -> float:
    """Return ``value`` as a float, rejecting NaN and +/- infinity."""
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:  # pragma: no cover - defensive
        raise TypeError(f"{name} must be a real number, got {value!r}") from exc
    if not math.isfinite(numeric):
        raise ValueError(f"{name} must be finite, got {numeric!r}")
    return numeric


def require_positive(value: float, name: str) -> float:
    """Return ``value`` as a finite float, rejecting values <= 0."""
    numeric = require_finite(value, name)
    if numeric <= 0.0:
        raise ValueError(f"{name} must be strictly positive, got {numeric!r}")
    return numeric


def require_non_negative(value: float, name: str) -> float:
    """Return ``value`` as a finite float, rejecting values < 0."""
    numeric = require_finite(value, name)
    if numeric < 0.0:
        raise ValueError(f"{name} must be non-negative, got {numeric!r}")
    return numeric
