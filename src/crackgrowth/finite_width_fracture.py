"""Finite-width fracture boundary: critical crack size and residual strength.

With a crack-size-dependent geometry factor the Milestone 2 closed form

    a_c = (1/pi) * (K_IC / (Y * sigma_max))**2

is NO LONGER VALID and is never used here: ``Y`` itself depends on ``a``, so the
fracture condition

    K_max(a_c) = Y(a_c) * sigma_max * sqrt(pi * a_c) = K_IC

is implicit and must be solved numerically. The constant-``Y`` formula remains
available in :mod:`crackgrowth.fracture` as the infinite-plate reference.

Root solver
-----------
A plain bounded bisection on

    f(a) = K_max(a) - K_IC

over ``[a_lower, a_upper]`` with ``a_upper < W/2`` strictly. Bisection is chosen
for transparency: it is deterministic, needs no derivative, and cannot leave its
bracket. The bracket is never widened silently -- if the sign condition fails,
an explicit status is returned instead of a number.

``f`` is monotonically increasing on the admissible interval (both
``sqrt(pi a)`` and ``Y(a)`` increase with ``a``), so the root is unique.

Never evaluating at W/2
-----------------------
``Y`` diverges as ``a -> W/2``. The upper bracket is placed at
``(W/2) * (1 - UPPER_BOUND_MARGIN)`` and the solver only ever bisects strictly
inside its bracket, so ``a = W/2`` is never evaluated. Because ``Y -> infinity``
there, any tensile ``sigma_max > 0`` with finite ``K_IC`` normally yields a root
below that bound -- but the "no root in bracket" status is still implemented and
tested rather than assumed away.

Residual strength
-----------------
    sigma_residual(a) = K_IC / ( sqrt(pi * a) * Y(a) )

Because ``Y(a) >= 1``, this is never greater than the infinite-plate value at the
same crack length, and it falls toward zero as ``a -> W/2``.

The finite-width critical crack length remains an LEFM screening boundary, not a
certified residual-strength allowable. No net-section-collapse check, no
plastic-zone validity check, and no thickness correction is performed; the
ligament quantities reported here are diagnostic only.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from ._validation import require_finite, require_positive
from .finite_width import (
    CrackGeometry,
    FiniteWidthCenterCrack,
    delta_k_for_geometry,
    geometry_factor_at,
    stress_intensity_for_geometry,
)
from .fracture import FractureToughness
from .loading import StressCycle

__all__ = [
    "DEFAULT_ROOT_TOLERANCE",
    "UPPER_BOUND_MARGIN",
    "CriticalCrackStatus",
    "FiniteWidthCriticalCrackResult",
    "FiniteWidthFractureAssessment",
    "k_max_for_geometry",
    "residual_strength_for_geometry",
    "finite_width_critical_crack_length",
    "assess_fracture_for_geometry",
]

#: Default absolute tolerance on the critical crack length [m].
DEFAULT_ROOT_TOLERANCE = 1.0e-12

#: Fractional standoff of the upper bracket from ``W/2``, where ``Y`` diverges.
UPPER_BOUND_MARGIN = 1.0e-12


class CriticalCrackStatus(Enum):
    """Outcome of the finite-width critical-crack-size solve."""

    ROOT_FOUND = "finite critical crack size found"
    NO_TENSILE_BOUNDARY = "no tensile mode-I fracture boundary (sigma_max <= 0)"
    ALREADY_CRITICAL_AT_LOWER_BOUND = (
        "K_max already exceeds K_IC at the lower search bound"
    )
    NO_ROOT_IN_BRACKET = "no root within the admissible interval"


@dataclass(frozen=True)
class FiniteWidthCriticalCrackResult:
    """The finite-width fracture boundary and its ligament diagnostics."""

    critical_crack_length: float
    """Half crack length at fracture [m]; ``math.inf`` when no boundary exists."""

    plate_width: float
    """Panel width ``W`` [m]."""

    geometry_factor_at_critical: float
    """``Y(a_c)`` [-]."""

    k_max_at_critical: float
    """``K_max(a_c)`` [Pa*sqrt(m)]; equals ``K_IC`` within solver tolerance."""

    k_ic: float
    """Mode-I fracture toughness [Pa*sqrt(m)]."""

    sigma_max: float
    """Maximum cycle stress [Pa]."""

    total_crack_length_at_critical: float
    """``2 * a_c`` [m]."""

    residual_ligament: float
    """``W - 2*a_c`` [m]. Diagnostic only."""

    ligament_fraction: float
    """``1 - 2*a_c/W`` [-]. Diagnostic only, NOT a collapse criterion."""

    iterations: int
    """Bisection iterations performed."""

    tolerance: float
    """Absolute tolerance on ``a_c`` requested [m]."""

    status: CriticalCrackStatus
    """Solver outcome."""

    @property
    def has_finite_boundary(self) -> bool:
        return self.status is CriticalCrackStatus.ROOT_FOUND


@dataclass(frozen=True)
class FiniteWidthFractureAssessment:
    """A geometry-aware fracture screen at one crack length."""

    crack_length: float
    geometry_factor: float
    k_max: float
    k_ic: float
    residual_strength: float
    applied_sigma_max: float
    utilization: float
    margin: float
    ligament_fraction: float | None
    """Diagnostic only; ``None`` for a geometry with no defined width."""
    passes: bool


def k_max_for_geometry(
    crack_length: float,
    cycle: StressCycle,
    geometry: CrackGeometry,
) -> float:
    """``K_max(a) = Y(a) * sigma_max * sqrt(pi * a)`` [Pa*sqrt(m)], signed.

    The FRACTURE driving force. Growth still uses ``delta_K``; the two are never
    interchanged. Note that because ``Y(a)`` multiplies both, the identity
    ``delta_K(a)/K_max(a) = delta_sigma/sigma_max`` survives unchanged for any
    geometry -- including a finite-width one.
    """
    if not isinstance(cycle, StressCycle):
        raise TypeError(f"cycle must be a StressCycle, got {cycle!r}")
    return stress_intensity_for_geometry(cycle.sigma_max, crack_length, geometry)


def residual_strength_for_geometry(
    crack_length: float,
    geometry: CrackGeometry,
    toughness: FractureToughness,
) -> float:
    """``sigma_residual(a) = K_IC / ( sqrt(pi * a) * Y(a) )`` [Pa]."""
    if not isinstance(toughness, FractureToughness):
        raise TypeError(f"toughness must be a FractureToughness, got {toughness!r}")
    a = require_positive(crack_length, "crack_length")
    return toughness.k_ic / (math.sqrt(math.pi * a) * geometry_factor_at(geometry, a))


def _ligament_fraction_or_none(
    geometry: CrackGeometry, crack_length: float
) -> float | None:
    if isinstance(geometry, FiniteWidthCenterCrack):
        return geometry.ligament_fraction(crack_length)
    return None


def finite_width_critical_crack_length(
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
    toughness: FractureToughness,
    lower_bound: float = 1.0e-9,
    tolerance: float = DEFAULT_ROOT_TOLERANCE,
    max_iterations: int = 200,
) -> FiniteWidthCriticalCrackResult:
    """Solve ``K_max(a_c) = K_IC`` for the finite-width panel by bisection.

    Parameters
    ----------
    lower_bound:
        Lower end of the search bracket [m]. Must be admissible.
    tolerance:
        Absolute tolerance on ``a_c`` [m]. Exposed, not hidden.
    max_iterations:
        Hard iteration cap; bisection halves the bracket each step.

    Returns
    -------
    FiniteWidthCriticalCrackResult
        Carries an explicit :class:`CriticalCrackStatus`. The bracket is never
        widened: a failed sign condition is reported, not worked around.
    """
    if not isinstance(cycle, StressCycle):
        raise TypeError(f"cycle must be a StressCycle, got {cycle!r}")
    if not isinstance(geometry, FiniteWidthCenterCrack):
        raise TypeError(
            f"geometry must be a FiniteWidthCenterCrack, got {geometry!r}"
        )
    if not isinstance(toughness, FractureToughness):
        raise TypeError(f"toughness must be a FractureToughness, got {toughness!r}")
    tol = require_positive(tolerance, "tolerance")
    lo = require_positive(lower_bound, "lower_bound")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int):
        raise TypeError(f"max_iterations must be an int, got {max_iterations!r}")
    if max_iterations < 1:
        raise ValueError(f"max_iterations must be at least 1, got {max_iterations!r}")

    width = geometry.plate_width

    def _no_boundary(status: CriticalCrackStatus) -> FiniteWidthCriticalCrackResult:
        return FiniteWidthCriticalCrackResult(
            critical_crack_length=math.inf,
            plate_width=width,
            geometry_factor_at_critical=math.inf,
            k_max_at_critical=math.inf,
            k_ic=toughness.k_ic,
            sigma_max=cycle.sigma_max,
            total_crack_length_at_critical=math.inf,
            residual_ligament=-math.inf,
            ligament_fraction=-math.inf,
            iterations=0,
            tolerance=tol,
            status=status,
        )

    # Compression never opens the crack: abs(sigma_max) is not taken.
    if cycle.sigma_max <= 0.0:
        return _no_boundary(CriticalCrackStatus.NO_TENSILE_BOUNDARY)

    # Upper bracket sits strictly inside W/2; Y is never evaluated at W/2 itself.
    hi = 0.5 * width * (1.0 - UPPER_BOUND_MARGIN)
    if lo >= hi:
        raise ValueError(
            f"lower_bound={lo!r} m is not below the admissible upper bracket "
            f"{hi!r} m for a panel of width {width!r} m"
        )

    def f(a: float) -> float:
        return k_max_for_geometry(a, cycle, geometry) - toughness.k_ic

    f_lo, f_hi = f(lo), f(hi)
    if f_lo >= 0.0:
        return _no_boundary(CriticalCrackStatus.ALREADY_CRITICAL_AT_LOWER_BOUND)
    if f_hi <= 0.0:
        # Y -> infinity at W/2 makes this practically unreachable, but it is
        # reported honestly rather than assumed impossible.
        return _no_boundary(CriticalCrackStatus.NO_ROOT_IN_BRACKET)

    iterations = 0
    while iterations < max_iterations and (hi - lo) > tol:
        mid = 0.5 * (lo + hi)
        iterations += 1
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    a_c = 0.5 * (lo + hi)

    return FiniteWidthCriticalCrackResult(
        critical_crack_length=a_c,
        plate_width=width,
        geometry_factor_at_critical=geometry.geometry_factor_at(a_c),
        k_max_at_critical=k_max_for_geometry(a_c, cycle, geometry),
        k_ic=toughness.k_ic,
        sigma_max=cycle.sigma_max,
        total_crack_length_at_critical=geometry.total_crack_length(a_c),
        residual_ligament=geometry.remaining_ligament(a_c),
        ligament_fraction=geometry.ligament_fraction(a_c),
        iterations=iterations,
        tolerance=tol,
        status=CriticalCrackStatus.ROOT_FOUND,
    )


def assess_fracture_for_geometry(
    crack_length: float,
    cycle: StressCycle,
    geometry: CrackGeometry,
    toughness: FractureToughness,
) -> FiniteWidthFractureAssessment:
    """Screen one crack length against the fracture boundary, for any geometry.

    ``ligament_fraction`` is reported alongside the fracture margin but is kept
    strictly separate from it: it is a geometric diagnostic, never blended into
    the margin and never used as a pass/fail criterion. ``passes`` depends on
    ``K_max <= K_IC`` alone.
    """
    if not isinstance(toughness, FractureToughness):
        raise TypeError(f"toughness must be a FractureToughness, got {toughness!r}")
    a = require_positive(crack_length, "crack_length")
    k = k_max_for_geometry(a, cycle, geometry)
    utilization = 0.0 if k <= 0.0 else k / toughness.k_ic
    margin = math.inf if k <= 0.0 else toughness.k_ic / k - 1.0
    return FiniteWidthFractureAssessment(
        crack_length=a,
        geometry_factor=geometry_factor_at(geometry, a),
        k_max=k,
        k_ic=toughness.k_ic,
        residual_strength=residual_strength_for_geometry(a, geometry, toughness),
        applied_sigma_max=require_finite(cycle.sigma_max, "sigma_max"),
        utilization=utilization,
        margin=margin,
        ligament_fraction=_ligament_fraction_or_none(geometry, a),
        passes=k <= toughness.k_ic,
    )


# Kept importable for callers building their own diagnostics.
_ = delta_k_for_geometry
