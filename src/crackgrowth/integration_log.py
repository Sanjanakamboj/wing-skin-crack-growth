"""Geometry-aware crack-growth integration on a logarithmic crack grid.

Why a second integrator
-----------------------
The Milestone 1 integrator uses a uniform grid in ``a`` and is a verified
baseline; it is NOT modified, replaced, or re-defaulted here. Milestone 2 showed
that a uniform-``a`` grid loses efficiency when ``a_f/a0`` is large, because
``dN/da`` is steepest just above ``a0``. Milestone 3 adds a second, independent
integrator that is better suited to wide spans and to a crack-size-dependent
geometry factor.

Change of variable
------------------
With ``u = ln(a)``, ``a = exp(u)`` and ``da = a du``:

    N = integral[a0 -> af] da / (da/dN)
      = integral[ln a0 -> ln af] a / (da/dN) du

A uniform grid in ``u`` is geometrically spaced in ``a``, so it places points
where the integrand actually varies. The integrand ``a / (da/dN)`` is evaluated
through the full modelling chain at every point:

    a -> Y(a) -> delta_K(a) -> da/dN(a) -> a / (da/dN)

``Y`` is therefore re-evaluated at every abscissa and is never factored out of
the integral or frozen at ``a0``. For a constant-``Y`` geometry this reduces
exactly to the Milestone 1 problem and must agree with its closed form.

Determinism and bounds
----------------------
Fixed-grid composite Simpson in ``u``: the analyst supplies the interval count,
there is no adaptive refinement and no hidden tolerance. An even count is
required (Simpson consumes sub-intervals in pairs) and an odd count is rejected
rather than silently adjusted. The endpoints are used verbatim and interior
abscissae are ``exp(u0 + i*h)``, which lies strictly inside ``[a0, af]`` for
``0 < i < n``, so the grid never leaves the interval. No SciPy is required.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ._validation import require_positive
from .finite_width import (
    CrackGeometry,
    crack_growth_rate_for_geometry,
    delta_k_for_geometry,
    geometry_factor_at,
)
from .loading import StressCycle
from .paris import ParisLaw

__all__ = [
    "DEFAULT_LOG_INTERVALS",
    "LogGridGrowthResult",
    "integrate_crack_growth_log_grid",
]

#: Default even interval count for the log-grid composite Simpson rule.
DEFAULT_LOG_INTERVALS = 1000


@dataclass(frozen=True)
class LogGridGrowthResult:
    """Deterministic log-grid integration result with endpoint diagnostics."""

    initial_crack_length: float
    final_crack_length: float
    predicted_cycles: float
    integration_intervals: int
    integration_method: str
    initial_geometry_factor: float
    final_geometry_factor: float
    initial_delta_k: float
    final_delta_k: float
    initial_growth_rate: float
    final_growth_rate: float


def integrate_crack_growth_log_grid(
    a_initial: float,
    a_final: float,
    cycle: StressCycle,
    geometry: CrackGeometry,
    paris_law: ParisLaw,
    intervals: int = DEFAULT_LOG_INTERVALS,
) -> LogGridGrowthResult:
    """Integrate Paris growth from ``a_initial`` to ``a_final`` in ``ln(a)``.

    Accepts any geometry satisfying
    :class:`~crackgrowth.finite_width.CrackGeometry`, so it serves both the
    constant-``Y`` Milestone 1/2 model and the finite-width ``Y(a)`` model.

    Parameters
    ----------
    intervals:
        Number of Simpson sub-intervals in ``u = ln(a)``. Even integer >= 2.

    Returns
    -------
    LogGridGrowthResult
        ``predicted_cycles`` is ``math.inf`` when the stress range is zero
        (the unchanged Milestone 1 no-growth policy).
    """
    if not isinstance(cycle, StressCycle):
        raise TypeError(f"cycle must be a StressCycle, got {cycle!r}")
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

    a0 = require_positive(a_initial, "a_initial")
    af = require_positive(a_final, "a_final")
    if af <= a0:
        raise ValueError(
            "a_final must be strictly greater than a_initial "
            f"(got a_initial={a0!r}, a_final={af!r})"
        )

    initial_geometry_factor = geometry_factor_at(geometry, a0)
    final_geometry_factor = geometry_factor_at(geometry, af)
    initial_delta_k = delta_k_for_geometry(a0, cycle, geometry)
    final_delta_k = delta_k_for_geometry(af, cycle, geometry)
    initial_rate = crack_growth_rate_for_geometry(a0, cycle, geometry, paris_law)
    final_rate = crack_growth_rate_for_geometry(af, cycle, geometry, paris_law)

    if cycle.stress_range == 0.0:
        predicted_cycles = math.inf
    else:

        def integrand(a: float) -> float:
            # dN/du = a * dN/da = a / (da/dN). Y is re-evaluated here at every
            # abscissa, never hoisted out of the integral.
            return a / crack_growth_rate_for_geometry(a, cycle, geometry, paris_law)

        u0 = math.log(a0)
        h = (math.log(af) - u0) / intervals
        total = integrand(a0) + integrand(af)
        for i in range(1, intervals):
            total += (4.0 if i % 2 == 1 else 2.0) * integrand(math.exp(u0 + i * h))
        predicted_cycles = total * h / 3.0

    return LogGridGrowthResult(
        initial_crack_length=a0,
        final_crack_length=af,
        predicted_cycles=predicted_cycles,
        integration_intervals=intervals,
        integration_method="composite Simpson on a log(a) grid",
        initial_geometry_factor=initial_geometry_factor,
        final_geometry_factor=final_geometry_factor,
        initial_delta_k=initial_delta_k,
        final_delta_k=final_delta_k,
        initial_growth_rate=initial_rate,
        final_growth_rate=final_rate,
    )
