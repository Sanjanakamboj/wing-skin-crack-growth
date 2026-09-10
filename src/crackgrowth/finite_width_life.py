"""Crack-growth life to the finite-width fracture boundary.

Composes verified pieces:

* the fracture boundary from :mod:`crackgrowth.finite_width_fracture`
  (bisection on ``K_max(a) = K_IC``, since no closed form exists for ``Y(a)``);
* the growth integration from :mod:`crackgrowth.integration_log`
  (geometry-aware log-grid Simpson, with ``Y`` re-evaluated at every abscissa).

Initial-flaw admissibility follows the Milestone 2 policy exactly, reusing
:class:`~crackgrowth.fracture_life.FlawAdmissibility`: ``a0 < a_c`` gives a
finite life, ``a0 == a_c`` gives zero, ``a0 > a_c`` gives zero with an explicit
status (or raises under ``strict=True``), and a non-opening cycle gives
``math.inf``. The integrator is never asked to run backwards.

Frozen-Y diagnostic
-------------------
:func:`frozen_geometry_factor_life` evaluates the life with ``Y`` held at
``Y(a0)`` for the whole integration. It exists ONLY to quantify the error of
that approximation, and is labelled as such everywhere. It is never used to
produce a headline result, because freezing ``Y`` at the initial flaw ignores
exactly the amplification that matters as the crack grows.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ._validation import require_positive
from .finite_width import (
    FiniteWidthCenterCrack,
    delta_k_for_geometry,
    geometry_factor_at,
)
from .finite_width_fracture import (
    CriticalCrackStatus,
    FiniteWidthCriticalCrackResult,
    FiniteWidthFractureAssessment,
    assess_fracture_for_geometry,
    finite_width_critical_crack_length,
)
from .fracture import FractureToughness
from .fracture_life import FlawAdmissibility, InitialFlawBeyondCriticalError
from .geometry import ThroughCrackGeometry
from .integration_log import (
    DEFAULT_LOG_INTERVALS,
    LogGridGrowthResult,
    integrate_crack_growth_log_grid,
)
from .loading import StressCycle
from .paris import ParisLaw

__all__ = [
    "FiniteWidthLifeResult",
    "cycles_to_finite_width_fracture",
    "frozen_geometry_factor_life",
]


@dataclass(frozen=True)
class FiniteWidthLifeResult:
    """Life from an initial flaw to the finite-width fracture boundary."""

    initial_crack_length: float
    """Assumed initial HALF crack length ``a0`` [m]."""

    plate_width: float
    """Panel width ``W`` [m]."""

    critical: FiniteWidthCriticalCrackResult
    """The solved fracture boundary and its diagnostics."""

    admissibility: FlawAdmissibility
    """Where ``a0`` sits relative to that boundary."""

    predicted_cycles: float
    """Life to fracture [cycles]."""

    growth_result: LogGridGrowthResult | None
    """Full integration diagnostics, or ``None`` when no integration ran."""

    initial_geometry_factor: float
    """``Y(a0)`` [-]."""

    initial_assessment: FiniteWidthFractureAssessment
    """Fracture screen at ``a0``."""

    critical_assessment: FiniteWidthFractureAssessment | None
    """Fracture screen at ``a_c``; ``None`` when no boundary exists."""

    delta_k_at_critical: float | None
    """``delta_K(a_c)`` [Pa*sqrt(m)]. Generally NOT ``K_IC``: for any geometry,
    ``delta_K(a_c) = (delta_sigma/sigma_max) * K_IC``, because ``Y(a_c)``
    multiplies both ``delta_K`` and ``K_max`` and cancels in the ratio."""

    @property
    def critical_crack_length(self) -> float:
        return self.critical.critical_crack_length

    @property
    def ligament_fraction_at_critical(self) -> float:
        """Diagnostic only; not a collapse criterion."""
        return self.critical.ligament_fraction

    @property
    def integration_method(self) -> str:
        if self.growth_result is None:
            return "none (no integration performed)"
        return self.growth_result.integration_method


def cycles_to_finite_width_fracture(
    a_initial: float,
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_LOG_INTERVALS,
    strict: bool = False,
    **solver_options: float,
) -> FiniteWidthLifeResult:
    """Solve the finite-width boundary and integrate Paris growth up to it.

    Parameters
    ----------
    intervals:
        Even interval count for the log-grid Simpson integration.
    strict:
        When True, ``a0 > a_c`` raises
        :class:`~crackgrowth.fracture_life.InitialFlawBeyondCriticalError`.
    **solver_options:
        Passed through to
        :func:`~crackgrowth.finite_width_fracture.finite_width_critical_crack_length`
        (``lower_bound``, ``tolerance``, ``max_iterations``).
    """
    if not isinstance(geometry, FiniteWidthCenterCrack):
        raise TypeError(
            f"geometry must be a FiniteWidthCenterCrack, got {geometry!r}"
        )
    if not isinstance(paris_law, ParisLaw):
        raise TypeError(f"paris_law must be a ParisLaw, got {paris_law!r}")
    a0 = require_positive(a_initial, "a_initial")
    # Rejects a0 >= W/2 with the finite-width admissibility message.
    initial_geometry_factor = geometry.geometry_factor_at(a0)

    critical = finite_width_critical_crack_length(
        cycle, geometry, toughness, **solver_options
    )
    initial_assessment = assess_fracture_for_geometry(
        a0, cycle, geometry, toughness
    )

    if critical.status is not CriticalCrackStatus.ROOT_FOUND:
        admissibility = (
            FlawAdmissibility.NO_TENSILE_BOUNDARY
            if critical.status is CriticalCrackStatus.NO_TENSILE_BOUNDARY
            else FlawAdmissibility.ABOVE_CRITICAL
        )
        no_boundary = (
            critical.status is CriticalCrackStatus.NO_TENSILE_BOUNDARY
        )
        return FiniteWidthLifeResult(
            initial_crack_length=a0,
            plate_width=geometry.plate_width,
            critical=critical,
            admissibility=admissibility,
            predicted_cycles=math.inf if no_boundary else 0.0,
            growth_result=None,
            initial_geometry_factor=initial_geometry_factor,
            initial_assessment=initial_assessment,
            critical_assessment=None,
            delta_k_at_critical=None,
        )

    a_c = critical.critical_crack_length
    critical_assessment = assess_fracture_for_geometry(
        a_c, cycle, geometry, toughness
    )
    delta_k_at_critical = delta_k_for_geometry(a_c, cycle, geometry)

    if a0 < a_c:
        admissibility = FlawAdmissibility.BELOW_CRITICAL
    elif a0 == a_c:
        admissibility = FlawAdmissibility.AT_CRITICAL
    else:
        admissibility = FlawAdmissibility.ABOVE_CRITICAL

    if admissibility is FlawAdmissibility.ABOVE_CRITICAL and strict:
        raise InitialFlawBeyondCriticalError(
            f"initial flaw a0={a0!r} m already exceeds the finite-width critical "
            f"crack length a_c={a_c!r} m for W={geometry.plate_width!r} m, "
            f"sigma_max={cycle.sigma_max!r} Pa and K_IC={toughness.k_ic!r} "
            "Pa*sqrt(m); no crack-growth life remains"
        )

    if admissibility is not FlawAdmissibility.BELOW_CRITICAL:
        return FiniteWidthLifeResult(
            initial_crack_length=a0,
            plate_width=geometry.plate_width,
            critical=critical,
            admissibility=admissibility,
            predicted_cycles=0.0,
            growth_result=None,
            initial_geometry_factor=initial_geometry_factor,
            initial_assessment=initial_assessment,
            critical_assessment=critical_assessment,
            delta_k_at_critical=delta_k_at_critical,
        )

    growth_result = integrate_crack_growth_log_grid(
        a0, a_c, cycle, geometry, paris_law, intervals=intervals
    )
    return FiniteWidthLifeResult(
        initial_crack_length=a0,
        plate_width=geometry.plate_width,
        critical=critical,
        admissibility=admissibility,
        predicted_cycles=growth_result.predicted_cycles,
        growth_result=growth_result,
        initial_geometry_factor=initial_geometry_factor,
        initial_assessment=initial_assessment,
        critical_assessment=critical_assessment,
        delta_k_at_critical=delta_k_at_critical,
    )


def frozen_geometry_factor_life(
    a_initial: float,
    a_final: float,
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    intervals: int = DEFAULT_LOG_INTERVALS,
) -> float:
    """DIAGNOSTIC ONLY: life with ``Y`` frozen at ``Y(a0)`` for the whole growth.

    This deliberately WRONG approximation exists solely to quantify how much the
    crack-size dependence of ``Y`` matters. Because ``Y`` grows with ``a``,
    freezing it at the initial flaw understates the driving force over most of
    the growth and therefore OVERESTIMATES the life. Never use it for a headline
    result or for any selection decision.
    """
    a0 = require_positive(a_initial, "a_initial")
    frozen = ThroughCrackGeometry(
        geometry_factor=geometry_factor_at(geometry, a0),
        description=(
            "DIAGNOSTIC ONLY: finite-width geometry factor frozen at Y(a0); "
            "not a valid crack-growth model"
        ),
    )
    return integrate_crack_growth_log_grid(
        a0, a_final, cycle, frozen, paris_law, intervals=intervals
    ).predicted_cycles
