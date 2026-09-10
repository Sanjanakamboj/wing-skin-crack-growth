"""Crack geometry representation.

Milestone 1 scope
-----------------
One deliberately simple geometry is supported: an idealized through-thickness
crack characterised by a single constant dimensionless geometry factor ``Y``.

The mode-I stress-intensity factor is then

    K = Y * sigma * sqrt(pi * a)

with ``a`` the crack length in metres (SI throughout).

The canonical Milestone 1 value is ``Y = 1.0``, which corresponds to the
classical central through crack of length ``2a`` in an *infinite* plate. No
finite-width correction, no surface/corner-crack shape factor, and no
crack-length-dependent ``Y(a)`` function is modelled here; ``Y`` is a constant
supplied by the analyst and is an explicit engineering assumption.
"""

from __future__ import annotations

from dataclasses import dataclass

from ._validation import require_positive

__all__ = ["ThroughCrackGeometry", "INFINITE_PLATE_THROUGH_CRACK"]


@dataclass(frozen=True)
class ThroughCrackGeometry:
    """An idealized through crack with a constant geometry factor.

    Parameters
    ----------
    geometry_factor:
        The dimensionless mode-I geometry (shape) factor ``Y``. Must be finite
        and strictly positive. ``Y = 1.0`` is the infinite-plate through crack.
    description:
        Free-text note recording what the factor is meant to represent. Purely
        documentary; it never affects a computed result.
    """

    geometry_factor: float
    description: str = "Idealized through crack, constant geometry factor"

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "geometry_factor",
            require_positive(self.geometry_factor, "geometry_factor"),
        )


INFINITE_PLATE_THROUGH_CRACK = ThroughCrackGeometry(
    geometry_factor=1.0,
    description=(
        "Central through crack of length 2a in an infinite plate (Y = 1.0). "
        "No finite-width correction is applied in Milestone 1."
    ),
)
