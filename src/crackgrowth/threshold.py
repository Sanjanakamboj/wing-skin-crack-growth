"""Fatigue crack-growth threshold (delta_K_th) screening.

Policy: a deliberately simple HARD CUTOFF
------------------------------------------
    delta_K <= delta_K_th   ->   da/dN = 0
    delta_K >  delta_K_th   ->   da/dN = C * delta_K ** m      (unchanged Paris)

Two properties of this choice matter and are stated rather than buried:

* **Equality belongs to the no-growth side.** A crack exactly at threshold does
  not grow.
* **The Paris rate is NOT modified above threshold.** This model does not use
  ``C*(delta_K - delta_K_th)**m`` nor ``C*(delta_K**m - delta_K_th**m)``; those
  are different near-threshold laws and are out of scope. Above threshold the
  rate is bit-for-bit the Milestone 1 Paris rate.

The consequence is a **discontinuous** model: the predicted life is either
identical to the no-threshold prediction (initial crack above threshold) or
infinite (initial crack at or below threshold). There is no gradual life
extension. That discontinuity is a property of the chosen simplified model, not
of the material, and it is documented rather than smoothed away.

No crack closure
----------------
The stress-intensity range is unchanged from Milestone 1:

    delta_K = Y(a) * (sigma_max - sigma_min) * sqrt(pi * a)

It is NOT modified by ``K_min``, by the stress ratio ``R``, by compression, or
by any crack-opening level. This is a threshold model, not a closure model. A
compressive ``sigma_min`` therefore still inflates the algebraic ``delta_K`` and
can make a crack look active when a closure-aware model would not.

Threshold crack size
--------------------
For constant ``Y`` the crack size at which ``delta_K`` first reaches threshold
has a closed form:

    a_th = (1/pi) * ( delta_K_th / (Y * delta_sigma) ) ** 2

For a finite-width panel ``Y = Y(a)`` and the condition is implicit, so it is
solved by the same bounded bisection used for the Milestone 3 fracture boundary.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from ._validation import require_positive
from .finite_width import (
    CrackGeometry,
    FiniteWidthCenterCrack,
    delta_k_for_geometry,
    geometry_factor_at,
)
from .finite_width_fracture import (
    DEFAULT_ROOT_TOLERANCE,
    UPPER_BOUND_MARGIN,
)
from .fracture import mpa_sqrt_m_to_pa_sqrt_m
from .geometry import ThroughCrackGeometry
from .loading import StressCycle
from .paris import ParisLaw

__all__ = [
    "ILLUSTRATIVE_THRESHOLD_DISCLAIMER",
    "ILLUSTRATIVE_ALUMINIUM_LIKE_THRESHOLD",
    "CrackGrowthThreshold",
    "ThresholdedGrowthPoint",
    "ThresholdBoundaryStatus",
    "FiniteWidthThresholdResult",
    "thresholded_crack_growth_rate",
    "threshold_crack_length_constant_y",
    "finite_width_threshold_crack_length",
]

#: Mandatory provenance banner for any threshold that is not a verified dataset.
ILLUSTRATIVE_THRESHOLD_DISCLAIMER = (
    "ILLUSTRATIVE CRACK-GROWTH THRESHOLD INPUT - NOT DESIGN ALLOWABLE"
)


@dataclass(frozen=True)
class CrackGrowthThreshold:
    """Immutable constant-amplitude crack-growth threshold, on the SI basis.

    Parameters
    ----------
    name:
        Short identifier for the record.
    delta_k_threshold:
        ``delta_K_th`` in Pa*sqrt(m). Finite, > 0. Use
        :func:`~crackgrowth.fracture.mpa_sqrt_m_to_pa_sqrt_m` to convert a
        handbook value quoted in MPa*sqrt(m) -- the same conversion helper the
        fracture toughness uses, so there is one unit-conversion path, not two.
    source_note:
        Provenance statement. For any value not traceable to a verified dataset
        this must carry :data:`ILLUSTRATIVE_THRESHOLD_DISCLAIMER`.
    notes:
        Optional additional commentary.

    Notes
    -----
    The threshold is treated as a single CONSTANT. It carries no dependence on
    stress ratio ``R``, environment, temperature, load history, or crack size.
    """

    name: str
    delta_k_threshold: float
    source_note: str
    notes: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("name must be a non-empty string")
        object.__setattr__(
            self,
            "delta_k_threshold",
            require_positive(self.delta_k_threshold, "delta_k_threshold"),
        )
        if not isinstance(self.source_note, str) or not self.source_note.strip():
            raise ValueError("source_note must be a non-empty provenance statement")
        if not isinstance(self.notes, str):
            raise TypeError("notes must be a string")

    @property
    def delta_k_threshold_in_mpa_sqrt_m(self) -> float:
        """The threshold re-expressed in MPa*sqrt(m), for reporting."""
        return self.delta_k_threshold / 1.0e6

    def is_active(self, delta_k: float) -> bool:
        """True when ``delta_k`` strictly exceeds the threshold.

        Equality returns False: a crack exactly at threshold does not grow.
        """
        return delta_k > self.delta_k_threshold

    def threshold_ratio(self, delta_k: float) -> float:
        """``delta_K / delta_K_th`` [-]. Growth is active only when this is > 1."""
        return delta_k / self.delta_k_threshold

    def threshold_excess_ratio(self, delta_k: float) -> float:
        """``delta_K / delta_K_th - 1`` [-].

        Deliberately NOT called a margin of safety: its sign convention is the
        opposite of a strength margin. POSITIVE means the crack is being driven
        ABOVE threshold and therefore grows -- which is the adverse case, not
        the safe one. Negative means arrested under this simplified model.
        """
        return self.threshold_ratio(delta_k) - 1.0


@dataclass(frozen=True)
class ThresholdedGrowthPoint:
    """Auditable threshold-aware growth-rate evaluation at one crack length."""

    crack_length: float
    geometry_factor: float
    delta_k: float
    delta_k_threshold: float
    threshold_ratio: float
    threshold_excess_ratio: float
    active_growth: bool
    growth_rate: float
    """``da/dN`` [m/cycle]: exactly zero when not active, otherwise the
    unmodified Paris rate."""


class ThresholdBoundaryStatus(Enum):
    """Outcome of the threshold-crack-size solve."""

    FINITE_THRESHOLD_SIZE = "finite threshold crack size found"
    ZERO_STRESS_RANGE = "zero stress range: threshold is never attainable"
    ALREADY_ABOVE_AT_LOWER_BOUND = (
        "delta_K already exceeds the threshold at the lower search bound"
    )
    NO_ROOT_IN_BRACKET = (
        "threshold not reached anywhere in the admissible interval"
    )


@dataclass(frozen=True)
class FiniteWidthThresholdResult:
    """The finite-width crack size at which growth first becomes active."""

    threshold_crack_length: float
    """Half crack length at which ``delta_K`` reaches ``delta_K_th`` [m];
    ``math.inf`` when no such crack size exists in the admissible interval."""

    plate_width: float
    geometry_factor_at_threshold: float
    delta_k_at_threshold: float
    delta_k_threshold: float
    stress_range: float
    total_crack_length_at_threshold: float
    residual_ligament: float
    ligament_fraction: float
    """Diagnostic only; never a criterion."""
    iterations: int
    tolerance: float
    status: ThresholdBoundaryStatus

    @property
    def has_finite_boundary(self) -> bool:
        return self.status is ThresholdBoundaryStatus.FINITE_THRESHOLD_SIZE


def thresholded_crack_growth_rate(
    crack_length: float,
    cycle: StressCycle,
    geometry: CrackGeometry,
    paris_law: ParisLaw,
    threshold: CrackGrowthThreshold,
) -> ThresholdedGrowthPoint:
    """Threshold-aware growth rate at one crack length.

    Returns a structured result so the classification is auditable rather than
    implicit in a bare number. Above threshold the rate is the UNMODIFIED Paris
    rate; at or below it, exactly ``0.0``.

    This is a NEW path: the original
    :func:`~crackgrowth.paris.crack_growth_rate` and
    :func:`~crackgrowth.finite_width.crack_growth_rate_for_geometry` are
    untouched and still apply no threshold at all.
    """
    if not isinstance(paris_law, ParisLaw):
        raise TypeError(f"paris_law must be a ParisLaw, got {paris_law!r}")
    if not isinstance(threshold, CrackGrowthThreshold):
        raise TypeError(
            f"threshold must be a CrackGrowthThreshold, got {threshold!r}"
        )
    a = require_positive(crack_length, "crack_length")
    delta_k = delta_k_for_geometry(a, cycle, geometry)
    active = threshold.is_active(delta_k)
    return ThresholdedGrowthPoint(
        crack_length=a,
        geometry_factor=geometry_factor_at(geometry, a),
        delta_k=delta_k,
        delta_k_threshold=threshold.delta_k_threshold,
        threshold_ratio=threshold.threshold_ratio(delta_k),
        threshold_excess_ratio=threshold.threshold_excess_ratio(delta_k),
        active_growth=active,
        growth_rate=paris_law.growth_rate(delta_k) if active else 0.0,
    )


def threshold_crack_length_constant_y(
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    threshold: CrackGrowthThreshold,
) -> float:
    """Exact threshold crack size for a CONSTANT geometry factor [m].

        a_th = (1/pi) * ( delta_K_th / (Y * delta_sigma) ) ** 2

    Returns ``math.inf`` when the stress range is zero, since ``delta_K`` is
    then zero at every crack size and the threshold is never attained.

    Valid only for constant ``Y``; a finite-width panel must use
    :func:`finite_width_threshold_crack_length`.
    """
    if not isinstance(cycle, StressCycle):
        raise TypeError(f"cycle must be a StressCycle, got {cycle!r}")
    if not isinstance(geometry, ThroughCrackGeometry):
        raise TypeError(
            f"geometry must be a ThroughCrackGeometry, got {geometry!r}"
        )
    if not isinstance(threshold, CrackGrowthThreshold):
        raise TypeError(
            f"threshold must be a CrackGrowthThreshold, got {threshold!r}"
        )
    if cycle.stress_range == 0.0:
        return math.inf
    return (
        threshold.delta_k_threshold
        / (geometry.geometry_factor * cycle.stress_range)
    ) ** 2 / math.pi


def finite_width_threshold_crack_length(
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
    threshold: CrackGrowthThreshold,
    lower_bound: float = 1.0e-12,
    tolerance: float = DEFAULT_ROOT_TOLERANCE,
    max_iterations: int = 200,
) -> FiniteWidthThresholdResult:
    """Solve ``Y(a) * delta_sigma * sqrt(pi * a) = delta_K_th`` by bisection.

    Uses the same bounded, deterministic, bracket-safe philosophy as the
    Milestone 3 fracture solver: the bracket is never widened silently, the
    tolerance is exposed, and the upper bound sits strictly inside ``W/2`` so
    the divergence of ``Y`` there is never evaluated.

    ``delta_K(a)`` is monotonically increasing in ``a`` for this geometry (both
    ``sqrt(pi a)`` and ``Y(a)`` increase), so the root is unique -- and once a
    crack is above threshold it stays above it as it grows.
    """
    if not isinstance(cycle, StressCycle):
        raise TypeError(f"cycle must be a StressCycle, got {cycle!r}")
    if not isinstance(geometry, FiniteWidthCenterCrack):
        raise TypeError(
            f"geometry must be a FiniteWidthCenterCrack, got {geometry!r}"
        )
    if not isinstance(threshold, CrackGrowthThreshold):
        raise TypeError(
            f"threshold must be a CrackGrowthThreshold, got {threshold!r}"
        )
    tol = require_positive(tolerance, "tolerance")
    lo = require_positive(lower_bound, "lower_bound")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int):
        raise TypeError(f"max_iterations must be an int, got {max_iterations!r}")
    if max_iterations < 1:
        raise ValueError(f"max_iterations must be at least 1, got {max_iterations!r}")

    width = geometry.plate_width

    def _no_boundary(status: ThresholdBoundaryStatus) -> FiniteWidthThresholdResult:
        return FiniteWidthThresholdResult(
            threshold_crack_length=math.inf,
            plate_width=width,
            geometry_factor_at_threshold=math.inf,
            delta_k_at_threshold=math.inf,
            delta_k_threshold=threshold.delta_k_threshold,
            stress_range=cycle.stress_range,
            total_crack_length_at_threshold=math.inf,
            residual_ligament=-math.inf,
            ligament_fraction=-math.inf,
            iterations=0,
            tolerance=tol,
            status=status,
        )

    if cycle.stress_range == 0.0:
        return _no_boundary(ThresholdBoundaryStatus.ZERO_STRESS_RANGE)

    hi = 0.5 * width * (1.0 - UPPER_BOUND_MARGIN)
    if lo >= hi:
        raise ValueError(
            f"lower_bound={lo!r} m is not below the admissible upper bracket "
            f"{hi!r} m for a panel of width {width!r} m"
        )

    def f(a: float) -> float:
        return delta_k_for_geometry(a, cycle, geometry) - threshold.delta_k_threshold

    if f(lo) >= 0.0:
        return _no_boundary(ThresholdBoundaryStatus.ALREADY_ABOVE_AT_LOWER_BOUND)
    if f(hi) <= 0.0:
        return _no_boundary(ThresholdBoundaryStatus.NO_ROOT_IN_BRACKET)

    iterations = 0
    while iterations < max_iterations and (hi - lo) > tol:
        mid = 0.5 * (lo + hi)
        iterations += 1
        if f(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    a_th = 0.5 * (lo + hi)

    return FiniteWidthThresholdResult(
        threshold_crack_length=a_th,
        plate_width=width,
        geometry_factor_at_threshold=geometry.geometry_factor_at(a_th),
        delta_k_at_threshold=delta_k_for_geometry(a_th, cycle, geometry),
        delta_k_threshold=threshold.delta_k_threshold,
        stress_range=cycle.stress_range,
        total_crack_length_at_threshold=geometry.total_crack_length(a_th),
        residual_ligament=geometry.remaining_ligament(a_th),
        ligament_fraction=geometry.ligament_fraction(a_th),
        iterations=iterations,
        tolerance=tol,
        status=ThresholdBoundaryStatus.FINITE_THRESHOLD_SIZE,
    )


#: Canonical Milestone 4 threshold.
#:
#: ILLUSTRATIVE. Selected after auditing delta_K_th = 2, 3, 4, 5, 6, 8 and
#: 10 MPa*sqrt(m) against the canonical finite-width case, where
#: delta_K(a0) = 5.6064 MPa*sqrt(m). Values up to 5 leave the initial crack
#: active; 6 and above arrest it. 4 MPa*sqrt(m) is a round value in the range
#: usually quoted for aluminium alloys at moderate stress ratio, gives a
#: threshold ratio of 1.40 at the initial flaw, and places the threshold crack
#: size at 0.51 mm -- comfortably below both a0 and the fracture boundary, so a
#: genuine active-growth interval exists. It was NOT chosen to guarantee active
#: growth: 5 MPa*sqrt(m) would also have been active.
ILLUSTRATIVE_ALUMINIUM_LIKE_THRESHOLD = CrackGrowthThreshold(
    name="Illustrative aluminium-like crack-growth threshold",
    delta_k_threshold=mpa_sqrt_m_to_pa_sqrt_m(4.0),
    source_note=(
        ILLUSTRATIVE_THRESHOLD_DISCLAIMER
        + ": delta_K_th = 4 MPa*sqrt(m), a round value chosen after auditing the "
        "resulting classification and threshold crack size for the canonical "
        "finite-width case. Not measured, not traceable to a qualification "
        "dataset, and no specific alloy is claimed."
    ),
    notes=(
        "Treated as a single constant: no stress-ratio, environmental, "
        "temperature or load-history dependence, and no near-threshold growth "
        "law. Hard cutoff only. Audited at W = 100 mm, delta_sigma = 100 MPa: "
        "a_th = 0.5092 mm, delta_K(a0)/delta_K_th = 1.4016."
    ),
)
