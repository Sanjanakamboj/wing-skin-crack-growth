"""Finite-width centre-cracked panel geometry.

Crack-length convention
-----------------------
Throughout this package ``a`` is the **HALF crack length**. A centre crack in a
panel of width ``W`` therefore has:

* total crack length ``2a``
* remaining total ligament ``W - 2a``
* ligament fraction ``1 - 2a/W``

and is physically admissible only while

    0 < 2a < W        equivalently    0 < a < W/2

This convention is stated everywhere it matters and is covered by dedicated
tests, because silently switching between ``a`` and ``2a`` is one of the easiest
and most damaging errors to make in this subject. The Milestone 1/2 constant-``Y``
model uses the same convention -- ``Y = 1`` is the *infinite-plate* centre crack
of total length ``2a`` -- so nothing about the meaning of ``a`` changes here.

Geometry factor
---------------
The standard secant finite-width correction for a centre-cracked panel:

    Y(a, W) = sqrt( sec( pi * a / W ) ) = 1 / sqrt( cos( pi * a / W ) )

Behaviour (all verified by test):

* ``Y -> 1`` as ``a/W -> 0`` (the infinite-plate limit);
* ``Y > 1`` for any finite non-zero ``a/W``;
* ``Y`` increases monotonically with ``a/W``;
* ``Y -> infinity`` as ``a -> W/2``, where the ligament vanishes.

Representative values: ``Y = 1.000247`` at ``a/W = 0.01``, ``1.0254`` at 0.10,
``1.1118`` at 0.20, ``1.3043`` at 0.30, ``1.7989`` at 0.40, ``5.6424`` at 0.49.

No other empirical correction factor is applied. The ligament quantities are
**diagnostic only**: this module contains no net-section-collapse model, and a
small ligament fraction is not treated as a failure criterion.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from ._validation import require_positive
from .loading import StressCycle
from .paris import ParisLaw

__all__ = [
    "CrackGeometry",
    "FiniteWidthCenterCrack",
    "total_crack_length",
    "remaining_ligament",
    "ligament_fraction",
    "geometry_factor_at",
    "stress_intensity_for_geometry",
    "delta_k_for_geometry",
    "crack_growth_rate_for_geometry",
]


@runtime_checkable
class CrackGeometry(Protocol):
    """Anything that can report a geometry factor at a given crack length.

    Both :class:`~crackgrowth.geometry.ThroughCrackGeometry` (constant ``Y``)
    and :class:`FiniteWidthCenterCrack` (``Y(a)``) satisfy this protocol, so the
    Milestone 3 geometry-aware routines accept either without special-casing.
    """

    def geometry_factor_at(self, crack_length: float) -> float:  # pragma: no cover
        ...


@dataclass(frozen=True)
class FiniteWidthCenterCrack:
    """A centre through crack in a finite-width panel.

    Parameters
    ----------
    plate_width:
        Panel width ``W`` [m]. Finite and strictly positive.
    description:
        Free-text note. Purely documentary; never affects a computed result.

    Notes
    -----
    ``a`` is the HALF crack length. Admissible crack lengths satisfy
    ``0 < a < W/2``; the exclusive upper bound is
    :attr:`max_admissible_crack_length`, at which ``Y`` diverges.
    """

    plate_width: float
    description: str = (
        "Centre through crack of total length 2a in a panel of width W; "
        "secant finite-width correction Y = 1/sqrt(cos(pi a / W))"
    )

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "plate_width", require_positive(self.plate_width, "plate_width")
        )

    @property
    def max_admissible_crack_length(self) -> float:
        """``W/2`` [m]: the EXCLUSIVE upper bound on the half crack length.

        At this value the ligament is zero and ``Y`` diverges; it is never a
        valid evaluation point.
        """
        return 0.5 * self.plate_width

    def _require_admissible(self, crack_length: float) -> float:
        a = require_positive(crack_length, "crack_length")
        if 2.0 * a >= self.plate_width:
            raise ValueError(
                "crack is not admissible in a finite-width panel: require "
                f"0 < 2a < W, but 2a={2.0 * a!r} m and W={self.plate_width!r} m "
                "(a is the HALF crack length)"
            )
        return a

    def geometry_factor_at(self, crack_length: float) -> float:
        """``Y(a) = 1 / sqrt(cos(pi * a / W))`` [-].

        Raises
        ------
        ValueError
            If ``2a >= W``. The bound is enforced rather than clamped: beyond it
            there is no ligament left and the correction is meaningless.
        """
        a = self._require_admissible(crack_length)
        return 1.0 / math.sqrt(math.cos(math.pi * a / self.plate_width))

    def total_crack_length(self, crack_length: float) -> float:
        """``2a`` [m] -- the full tip-to-tip crack length."""
        return 2.0 * self._require_admissible(crack_length)

    def remaining_ligament(self, crack_length: float) -> float:
        """``W - 2a`` [m]. Diagnostic only; not a collapse criterion."""
        return self.plate_width - 2.0 * self._require_admissible(crack_length)

    def ligament_fraction(self, crack_length: float) -> float:
        """``1 - 2a/W`` [-], strictly between 0 and 1. Diagnostic only."""
        a = self._require_admissible(crack_length)
        return 1.0 - 2.0 * a / self.plate_width


def total_crack_length(crack_length: float) -> float:
    """``2a`` [m] for the half-crack-length convention."""
    return 2.0 * require_positive(crack_length, "crack_length")


def remaining_ligament(crack_length: float, plate_width: float) -> float:
    """``W - 2a`` [m]. Diagnostic only."""
    a = require_positive(crack_length, "crack_length")
    w = require_positive(plate_width, "plate_width")
    return w - 2.0 * a


def ligament_fraction(crack_length: float, plate_width: float) -> float:
    """``1 - 2a/W`` [-]. Diagnostic only; not a net-section-collapse criterion."""
    a = require_positive(crack_length, "crack_length")
    w = require_positive(plate_width, "plate_width")
    return 1.0 - 2.0 * a / w


def geometry_factor_at(geometry: CrackGeometry, crack_length: float) -> float:
    """Geometry factor at ``crack_length`` for any supported geometry."""
    factor_at = getattr(geometry, "geometry_factor_at", None)
    if factor_at is None:
        raise TypeError(
            f"geometry must provide geometry_factor_at(crack_length), got {geometry!r}"
        )
    return factor_at(crack_length)


def stress_intensity_for_geometry(
    stress: float,
    crack_length: float,
    geometry: CrackGeometry,
) -> float:
    """Signed mode-I stress intensity ``K = Y(a) * sigma * sqrt(pi * a)``.

    This is the SAME equation as the Milestone 1 helper; only the source of
    ``Y`` differs. ``K`` remains signed: negative for compressive ``stress``, and
    no absolute value is taken.
    """
    from ._validation import require_finite

    sigma = require_finite(stress, "stress")
    a = require_positive(crack_length, "crack_length")
    return geometry_factor_at(geometry, a) * sigma * math.sqrt(math.pi * a)


def delta_k_for_geometry(
    crack_length: float,
    cycle: StressCycle,
    geometry: CrackGeometry,
) -> float:
    """``delta_K = Y(a) * delta_sigma * sqrt(pi * a)`` [Pa*sqrt(m)].

    The algebraic stress range is used, exactly as in Milestone 1: crack closure
    is still not modelled, so a compressive ``sigma_min`` contributes in full.
    Because ``Y(a) >= 1`` for a finite-width panel, this is never smaller than
    the infinite-plate value at the same crack length.
    """
    if not isinstance(cycle, StressCycle):
        raise TypeError(f"cycle must be a StressCycle, got {cycle!r}")
    a = require_positive(crack_length, "crack_length")
    return geometry_factor_at(geometry, a) * cycle.stress_range * math.sqrt(
        math.pi * a
    )


def crack_growth_rate_for_geometry(
    crack_length: float,
    cycle: StressCycle,
    geometry: CrackGeometry,
    paris_law: ParisLaw,
) -> float:
    """``da/dN = C * delta_K(a)**m`` [m/cycle] with a geometry-dependent ``Y``.

    The Paris law itself is the unchanged Milestone 1 model; only ``delta_K``
    now carries the crack-size dependence of ``Y``.
    """
    if not isinstance(paris_law, ParisLaw):
        raise TypeError(f"paris_law must be a ParisLaw, got {paris_law!r}")
    return paris_law.growth_rate(delta_k_for_geometry(crack_length, cycle, geometry))
