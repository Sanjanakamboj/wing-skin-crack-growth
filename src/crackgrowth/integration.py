"""Cycles-to-target crack growth: analytical reference and numerical quadrature.

Both routines answer the same question for constant-amplitude loading and a
CONSTANT geometry factor ``Y``:

    N = integral from a0 to af of  da / (da/dN)

with ``da/dN = C * (Y * delta_sigma * sqrt(pi * a)) ** m``.

Analytical reference (constant Y, constant delta_sigma)
-------------------------------------------------------
Separating the constants, ``da/dN = C * (Y * delta_sigma * sqrt(pi)) ** m * a ** (m/2)``,
so

    dN/da = 1 / [ C * (Y * delta_sigma * sqrt(pi)) ** m * a ** (m/2) ]

Integrating with ``p = 1 - m/2``:

    m != 2:   N = (af**p - a0**p) / ( p * C * (Y * delta_sigma * sqrt(pi)) ** m )
    m == 2:   N = ln(af / a0) / ( C * (Y * delta_sigma * sqrt(pi)) ** 2 )

Sign note for ``m > 2``: then ``p < 0``, so ``af**p < a0**p`` and the numerator is
negative while the denominator (which carries the factor ``p``) is also negative.
The quotient is positive, as a life must be. This is verified by test.

Numerical quadrature
--------------------
A fixed-grid composite Simpson rule is used. It is fully deterministic: the
analyst supplies the interval count, there is no adaptive refinement and no
hidden tolerance. Simpson's rule pairs sub-intervals, so an EVEN interval count
is required and an odd count is rejected rather than silently adjusted. The grid
spans exactly ``[a0, af]`` and never evaluates outside it.

Zero-range policy
-----------------
If ``delta_sigma == 0`` then ``delta_K == 0`` and ``da/dN == 0``: the crack never
reaches the target. Both routines return ``math.inf`` (an explicit, documented
no-growth result) rather than dividing by zero.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ._validation import require_positive
from .geometry import ThroughCrackGeometry
from .loading import StressCycle
from .paris import ParisLaw, crack_growth_rate
from .stress_intensity import delta_stress_intensity

__all__ = [
    "DEFAULT_INTEGRATION_INTERVALS",
    "CrackGrowthResult",
    "analytical_cycles_to_crack_length",
    "cycles_to_crack_length",
]

#: Default even interval count for the composite Simpson rule.
DEFAULT_INTEGRATION_INTERVALS = 1000


def _validate_crack_range(a_initial: float, a_final: float) -> tuple[float, float]:
    a0 = require_positive(a_initial, "a_initial")
    af = require_positive(a_final, "a_final")
    if af <= a0:
        raise ValueError(
            "a_final must be strictly greater than a_initial "
            f"(got a_initial={a0!r}, a_final={af!r})"
        )
    return a0, af


def _growth_constant(
    cycle: StressCycle, geometry: ThroughCrackGeometry, paris_law: ParisLaw
) -> float:
    """``C * (Y * delta_sigma * sqrt(pi)) ** m``, the a-independent rate factor."""
    return (
        paris_law.C
        * (geometry.geometry_factor * cycle.stress_range * math.sqrt(math.pi))
        ** paris_law.m
    )


def analytical_cycles_to_crack_length(
    a_initial: float,
    a_final: float,
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    paris_law: ParisLaw,
) -> float:
    """Closed-form cycles from ``a_initial`` to ``a_final`` [cycles].

    Valid only for a constant geometry factor and a constant stress range -- the
    Milestone 1 problem. This function is an INDEPENDENT reference: it evaluates
    the closed form directly and never calls the numerical integrator.

    Returns ``math.inf`` when the stress range is zero.
    """
    if not isinstance(cycle, StressCycle):
        raise TypeError(f"cycle must be a StressCycle, got {cycle!r}")
    if not isinstance(geometry, ThroughCrackGeometry):
        raise TypeError(f"geometry must be a ThroughCrackGeometry, got {geometry!r}")
    if not isinstance(paris_law, ParisLaw):
        raise TypeError(f"paris_law must be a ParisLaw, got {paris_law!r}")

    a0, af = _validate_crack_range(a_initial, a_final)
    if cycle.stress_range == 0.0:
        return math.inf

    constant = _growth_constant(cycle, geometry, paris_law)
    m = paris_law.m
    if m == 2.0:
        return math.log(af / a0) / constant
    p = 1.0 - 0.5 * m
    return (af**p - a0**p) / (p * constant)


@dataclass(frozen=True)
class CrackGrowthResult:
    """Deterministic cycles-to-target result with auditable endpoint values."""

    initial_crack_length: float
    final_crack_length: float
    predicted_cycles: float
    integration_intervals: int
    initial_delta_k: float
    final_delta_k: float
    initial_growth_rate: float
    final_growth_rate: float


def cycles_to_crack_length(
    a_initial: float,
    a_final: float,
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    paris_law: ParisLaw,
    intervals: int = DEFAULT_INTEGRATION_INTERVALS,
) -> CrackGrowthResult:
    """Numerically integrate crack growth from ``a_initial`` to ``a_final``.

    Uses a fixed-grid composite Simpson rule on ``1 / (da/dN)`` over
    ``[a_initial, a_final]``.

    Parameters
    ----------
    intervals:
        Number of Simpson sub-intervals. Must be an even integer >= 2. An odd
        count is rejected: Simpson's rule consumes sub-intervals in pairs, and
        silently incrementing the count would make the result depend on an
        undocumented adjustment.

    Returns
    -------
    CrackGrowthResult
        ``predicted_cycles`` is ``math.inf`` when the stress range is zero.
    """
    if not isinstance(cycle, StressCycle):
        raise TypeError(f"cycle must be a StressCycle, got {cycle!r}")
    if not isinstance(geometry, ThroughCrackGeometry):
        raise TypeError(f"geometry must be a ThroughCrackGeometry, got {geometry!r}")
    if not isinstance(paris_law, ParisLaw):
        raise TypeError(f"paris_law must be a ParisLaw, got {paris_law!r}")
    if isinstance(intervals, bool) or not isinstance(intervals, int):
        raise TypeError(f"intervals must be an int, got {intervals!r}")
    if intervals < 2:
        raise ValueError(f"intervals must be at least 2, got {intervals!r}")
    if intervals % 2 != 0:
        raise ValueError(
            f"intervals must be even for the composite Simpson rule, got {intervals!r}"
        )

    a0, af = _validate_crack_range(a_initial, a_final)

    initial_delta_k = delta_stress_intensity(a0, cycle, geometry)
    final_delta_k = delta_stress_intensity(af, cycle, geometry)
    initial_growth_rate = crack_growth_rate(a0, cycle, geometry, paris_law)
    final_growth_rate = crack_growth_rate(af, cycle, geometry, paris_law)

    if cycle.stress_range == 0.0:
        predicted_cycles = math.inf
    else:

        def integrand(a: float) -> float:
            # dN/da = 1 / (da/dN). The rate is obtained from the full modelling
            # chain (delta_K then Paris law) rather than from the closed-form
            # factored constant, so that this quadrature is genuinely
            # independent of analytical_cycles_to_crack_length and the two can
            # be cross-checked against each other.
            return 1.0 / crack_growth_rate(a, cycle, geometry, paris_law)

        h = (af - a0) / intervals
        total = integrand(a0) + integrand(af)
        for i in range(1, intervals):
            # Grid points are constructed from a0 + i*h and the endpoints are
            # used verbatim, so no abscissa ever falls outside [a0, af].
            a = a0 + i * h
            total += (4.0 if i % 2 == 1 else 2.0) * integrand(a)
        predicted_cycles = total * h / 3.0

    return CrackGrowthResult(
        initial_crack_length=a0,
        final_crack_length=af,
        predicted_cycles=predicted_cycles,
        integration_intervals=intervals,
        initial_delta_k=initial_delta_k,
        final_delta_k=final_delta_k,
        initial_growth_rate=initial_growth_rate,
        final_growth_rate=final_growth_rate,
    )
