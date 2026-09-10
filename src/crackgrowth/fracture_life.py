"""Crack-growth life from an initial flaw to the toughness-derived fracture size.

This module composes verified pieces rather than restating any physics:

* the fracture boundary comes from :mod:`crackgrowth.fracture`;
* the growth integration and its closed-form reference are the UNCHANGED
  Milestone 1 routines :func:`~crackgrowth.integration.cycles_to_crack_length`
  and :func:`~crackgrowth.integration.analytical_cycles_to_crack_length`,
  evaluated with ``a_final = a_critical``.

No separate "fracture growth equation" exists or should exist: Milestone 2 only
moves the upper integration limit from an imposed value to a physics-derived one.

Initial-flaw admissibility
--------------------------
Before integrating, the initial flaw is compared with the critical size:

============================  ==================================================
``a0 < a_c``                  finite crack-growth life exists
``a0 == a_c``                 zero remaining cycles; already at the boundary
``a0 > a_c``                  flaw is already beyond the fracture boundary
no tensile boundary           ``sigma_max <= 0``; no mode-I fracture size exists
============================  ==================================================

The integrator is never asked to run backwards. For ``a0 >= a_c`` the reported
life is ``0.0`` and the status distinguishes the two cases: for ``a0 > a_c`` the
zero means "no crack-growth life remains -- the flaw is outside the LEFM
fracture boundary at the first application of ``sigma_max``", NOT "it grew to
critical in zero cycles".
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from ._validation import require_positive
from .fracture import (
    CriticalCrackResult,
    FractureAssessment,
    FractureBoundaryStatus,
    FractureToughness,
    assess_fracture,
    critical_crack_length,
)
from .geometry import ThroughCrackGeometry
from .integration import (
    DEFAULT_INTEGRATION_INTERVALS,
    CrackGrowthResult,
    analytical_cycles_to_crack_length,
    cycles_to_crack_length,
)
from .loading import StressCycle
from .paris import ParisLaw
from .stress_intensity import delta_stress_intensity

__all__ = [
    "FlawAdmissibility",
    "InitialFlawBeyondCriticalError",
    "LifeToFractureResult",
    "assess_initial_flaw",
    "cycles_to_critical_crack",
]


class FlawAdmissibility(Enum):
    """Where the initial flaw sits relative to the fracture boundary."""

    BELOW_CRITICAL = "initial flaw below critical size"
    AT_CRITICAL = "initial flaw exactly at critical size"
    ABOVE_CRITICAL = "initial flaw already beyond the fracture boundary"
    NO_TENSILE_BOUNDARY = "no tensile mode-I fracture boundary exists"


class InitialFlawBeyondCriticalError(ValueError):
    """Raised only in strict mode when ``a0 > a_critical``.

    The default behaviour of :func:`cycles_to_critical_crack` is to report this
    condition as a structured status instead of raising; pass ``strict=True`` to
    make it an error at the call site.
    """


def assess_initial_flaw(
    a_initial: float,
    critical: CriticalCrackResult,
) -> FlawAdmissibility:
    """Classify an initial flaw against a computed fracture boundary."""
    if not isinstance(critical, CriticalCrackResult):
        raise TypeError(f"critical must be a CriticalCrackResult, got {critical!r}")
    a0 = require_positive(a_initial, "a_initial")
    if critical.status is FractureBoundaryStatus.NO_TENSILE_BOUNDARY:
        return FlawAdmissibility.NO_TENSILE_BOUNDARY
    a_c = critical.critical_crack_length
    if a0 < a_c:
        return FlawAdmissibility.BELOW_CRITICAL
    if a0 == a_c:
        return FlawAdmissibility.AT_CRITICAL
    return FlawAdmissibility.ABOVE_CRITICAL


@dataclass(frozen=True)
class LifeToFractureResult:
    """Cycles from an initial flaw to the toughness-derived critical size."""

    initial_crack_length: float
    """Assumed initial flaw size ``a0`` [m]."""

    critical: CriticalCrackResult
    """The fracture boundary, including its status."""

    admissibility: FlawAdmissibility
    """Where ``a0`` sits relative to that boundary."""

    predicted_cycles: float
    """Numerically integrated life [cycles]. See the module docstring for the
    meaning of ``0.0`` under each admissibility status."""

    analytical_cycles: float
    """Independent closed-form life [cycles], from the Milestone 1 reference."""

    growth_result: CrackGrowthResult | None
    """Full Milestone 1 integration diagnostics, or ``None`` when no integration
    was performed (flaw at/beyond critical, or no tensile boundary)."""

    initial_assessment: FractureAssessment
    """Fracture screen at ``a0``."""

    critical_assessment: FractureAssessment | None
    """Fracture screen at ``a_critical``; ``None`` when no boundary exists."""

    delta_k_at_critical: float | None
    """``delta_K(a_critical)`` [Pa*sqrt(m)]; ``None`` when no boundary exists.

    Note this is generally NOT equal to ``K_IC``: for constant ``Y``,
    ``delta_K(a_c) = (delta_sigma / sigma_max) * K_IC``."""

    @property
    def critical_crack_length(self) -> float:
        """Convenience accessor for the critical crack length [m]."""
        return self.critical.critical_crack_length

    @property
    def has_finite_life(self) -> bool:
        return math.isfinite(self.predicted_cycles)


def cycles_to_critical_crack(
    a_initial: float,
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_INTEGRATION_INTERVALS,
    strict: bool = False,
) -> LifeToFractureResult:
    """Integrate Paris growth from ``a_initial`` to the critical crack size.

    Steps: compute ``a_c`` from ``K_IC`` and ``sigma_max``; classify the initial
    flaw; integrate the (unchanged) Milestone 1 growth model over
    ``[a0, a_c]``; report the fracture endpoint alongside the full integration
    diagnostics.

    Parameters
    ----------
    strict:
        When True, ``a_initial > a_critical`` raises
        :class:`InitialFlawBeyondCriticalError` instead of being reported as a
        structured status. Default False: prefer the structured result.

    Notes
    -----
    A zero stress range yields ``math.inf`` cycles (the Milestone 1 no-growth
    policy) even though a finite critical size exists: the crack simply never
    grows to it.
    """
    if not isinstance(paris_law, ParisLaw):
        raise TypeError(f"paris_law must be a ParisLaw, got {paris_law!r}")
    a0 = require_positive(a_initial, "a_initial")
    critical = critical_crack_length(cycle, geometry, toughness)
    admissibility = assess_initial_flaw(a0, critical)
    initial_assessment = assess_fracture(a0, cycle, geometry, toughness)

    if admissibility is FlawAdmissibility.NO_TENSILE_BOUNDARY:
        return LifeToFractureResult(
            initial_crack_length=a0,
            critical=critical,
            admissibility=admissibility,
            predicted_cycles=math.inf,
            analytical_cycles=math.inf,
            growth_result=None,
            initial_assessment=initial_assessment,
            critical_assessment=None,
            delta_k_at_critical=None,
        )

    a_c = critical.critical_crack_length
    critical_assessment = assess_fracture(a_c, cycle, geometry, toughness)
    delta_k_at_critical = delta_stress_intensity(a_c, cycle, geometry)

    if admissibility is FlawAdmissibility.ABOVE_CRITICAL and strict:
        raise InitialFlawBeyondCriticalError(
            f"initial flaw a0={a0!r} m already exceeds the critical crack length "
            f"a_c={a_c!r} m for sigma_max={cycle.sigma_max!r} Pa and "
            f"K_IC={toughness.k_ic!r} Pa*sqrt(m); no crack-growth life remains"
        )

    if admissibility in (
        FlawAdmissibility.AT_CRITICAL,
        FlawAdmissibility.ABOVE_CRITICAL,
    ):
        # No integration: the interval is empty or reversed. Reporting 0.0 here
        # keeps the Milestone 1 integrator's "a_final > a_initial" contract intact.
        return LifeToFractureResult(
            initial_crack_length=a0,
            critical=critical,
            admissibility=admissibility,
            predicted_cycles=0.0,
            analytical_cycles=0.0,
            growth_result=None,
            initial_assessment=initial_assessment,
            critical_assessment=critical_assessment,
            delta_k_at_critical=delta_k_at_critical,
        )

    growth_result = cycles_to_crack_length(
        a0, a_c, cycle, geometry, paris_law, intervals=intervals
    )
    analytical = analytical_cycles_to_crack_length(
        a0, a_c, cycle, geometry, paris_law
    )
    return LifeToFractureResult(
        initial_crack_length=a0,
        critical=critical,
        admissibility=admissibility,
        predicted_cycles=growth_result.predicted_cycles,
        analytical_cycles=analytical,
        growth_result=growth_result,
        initial_assessment=initial_assessment,
        critical_assessment=critical_assessment,
        delta_k_at_critical=delta_k_at_critical,
    )
