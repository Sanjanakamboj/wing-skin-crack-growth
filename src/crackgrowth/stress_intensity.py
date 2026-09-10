"""Mode-I stress-intensity factor and stress-intensity range.

Governing relation (Milestone 1, constant geometry factor)
----------------------------------------------------------
    K = Y * sigma * sqrt(pi * a)

Units are SI throughout: ``sigma`` in Pa, ``a`` in m, ``K`` in Pa*sqrt(m).

Signs
-----
``stress_intensity`` is SIGNED and returns the algebraic product above. For
tensile ``sigma`` it is positive; for compressive ``sigma`` it is negative. No
absolute value is taken. A negative ``K`` is not physically an opening-mode
intensity -- a real crack closes and transmits compression -- so the negative
value should be read as a bookkeeping result, not as a driving force.

Stress-intensity range
----------------------
    delta_K = Y * delta_sigma * sqrt(pi * a),   delta_sigma = sigma_max - sigma_min

This is the intentionally simple algebraic range used by the canonical Paris-law
screen. Because Milestone 1 does NOT model crack closure, a compressive
``sigma_min`` still contributes in full to ``delta_sigma`` and therefore inflates
``delta_K``. Real aluminium behaviour would partially or wholly remove the
compressive part of the cycle. No effective-``delta_K`` model, no threshold
``delta_K_th``, and no fracture-toughness cutoff is applied here.
"""

from __future__ import annotations

import math

from ._validation import require_finite, require_positive
from .geometry import ThroughCrackGeometry
from .loading import StressCycle

__all__ = ["stress_intensity", "delta_stress_intensity"]


def stress_intensity(
    stress: float,
    crack_length: float,
    geometry: ThroughCrackGeometry,
) -> float:
    """Signed mode-I stress-intensity factor ``K = Y * sigma * sqrt(pi * a)``.

    Parameters
    ----------
    stress:
        Remote stress [Pa]. Tension positive, compression negative.
    crack_length:
        Crack length ``a`` [m]. Must be finite and strictly positive.
    geometry:
        Crack geometry supplying the constant factor ``Y``.

    Returns
    -------
    float
        Stress intensity [Pa*sqrt(m)], carrying the sign of ``stress``.
    """
    if not isinstance(geometry, ThroughCrackGeometry):
        raise TypeError(f"geometry must be a ThroughCrackGeometry, got {geometry!r}")
    sigma = require_finite(stress, "stress")
    a = require_positive(crack_length, "crack_length")
    return geometry.geometry_factor * sigma * math.sqrt(math.pi * a)


def delta_stress_intensity(
    crack_length: float,
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
) -> float:
    """Stress-intensity range ``delta_K = Y * delta_sigma * sqrt(pi * a)``.

    ``delta_sigma`` is the algebraic range ``sigma_max - sigma_min``; no closure
    correction is applied (see the module docstring). The result is always
    non-negative because :class:`~crackgrowth.loading.StressCycle` enforces
    ``sigma_min <= sigma_max``.

    Parameters
    ----------
    crack_length:
        Crack length ``a`` [m]. Must be finite and strictly positive.
    cycle:
        The constant-amplitude stress cycle.
    geometry:
        Crack geometry supplying the constant factor ``Y``.

    Returns
    -------
    float
        Stress-intensity range [Pa*sqrt(m)]. Zero when ``delta_sigma`` is zero.
    """
    if not isinstance(cycle, StressCycle):
        raise TypeError(f"cycle must be a StressCycle, got {cycle!r}")
    if not isinstance(geometry, ThroughCrackGeometry):
        raise TypeError(f"geometry must be a ThroughCrackGeometry, got {geometry!r}")
    a = require_positive(crack_length, "crack_length")
    return geometry.geometry_factor * cycle.stress_range * math.sqrt(math.pi * a)
