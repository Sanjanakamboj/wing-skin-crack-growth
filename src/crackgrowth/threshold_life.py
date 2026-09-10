"""Threshold-aware life to the finite-width fracture boundary.

This module only classifies and dispatches; it introduces no new mechanics.
The active-growth branch calls the UNCHANGED Milestone 3 helper
:func:`~crackgrowth.finite_width_life.cycles_to_finite_width_fracture`, which in
turn uses the unchanged log-grid integrator and the unchanged Paris rate.

The central correctness point
-----------------------------
Under constant-amplitude loading a crack whose ``delta_K`` is at or below
``delta_K_th`` has ``da/dN = 0``. Because it does not grow at all, it can NEVER
reach the threshold crack size ``a_th`` on its own. The remaining crack-growth
life is therefore ``+infinity``.

It is emphatically NOT "the life from ``a_th`` to ``a_c``". The crack is never
advanced to ``a_th``, and no such life is ever reported. Reporting one would
silently teleport the crack past the arrest that the model just predicted.

"Infinite life" here means ONLY that this model predicts zero propagation under
this exact constant-amplitude cycle. It is not a statement about total
structural durability, and a threshold-arrested result is not a safe-life
certification result.

The discontinuity
-----------------
Because the hard cutoff does not modify the Paris rate above threshold, the
predicted life is either

* exactly the no-threshold Milestone 3 life (initial crack above threshold), or
* infinite (initial crack at or below threshold).

There is no intermediate regime and no gradual life extension. See
:mod:`crackgrowth.threshold`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from ._validation import require_positive
from .finite_width import FiniteWidthCenterCrack, delta_k_for_geometry
from .finite_width_fracture import (
    CriticalCrackStatus,
    FiniteWidthCriticalCrackResult,
    finite_width_critical_crack_length,
)
from .finite_width_life import cycles_to_finite_width_fracture
from .fracture import FractureToughness
from .integration_log import DEFAULT_LOG_INTERVALS, LogGridGrowthResult
from .loading import StressCycle
from .paris import ParisLaw
from .threshold import (
    CrackGrowthThreshold,
    FiniteWidthThresholdResult,
    ThresholdBoundaryStatus,
    finite_width_threshold_crack_length,
)

__all__ = [
    "GrowthState",
    "BoundaryOrdering",
    "ThresholdedLifeResult",
    "classify_growth_state",
    "compare_threshold_and_fracture_boundaries",
    "cycles_to_fracture_with_threshold",
]


class GrowthState(Enum):
    """What the model says the initial crack is actually doing."""

    ACTIVE_GROWTH = "growing: delta_K exceeds the threshold"
    AT_THRESHOLD = "arrested: delta_K is exactly at the threshold"
    ARRESTED_BELOW_THRESHOLD = "arrested: delta_K is below the threshold"
    AT_FRACTURE_BOUNDARY = "at the fracture boundary"
    ABOVE_FRACTURE_BOUNDARY = "already beyond the fracture boundary"
    NO_TENSILE_FRACTURE_BOUNDARY = "no tensile mode-I fracture boundary exists"

    @property
    def is_growing(self) -> bool:
        return self is GrowthState.ACTIVE_GROWTH

    @property
    def is_arrested(self) -> bool:
        return self in (
            GrowthState.AT_THRESHOLD,
            GrowthState.ARRESTED_BELOW_THRESHOLD,
        )


class BoundaryOrdering(Enum):
    """How the threshold crack size sits relative to the fracture boundary."""

    ACTIVE_INTERVAL_EXISTS = "a_th < a_c: a stable growth interval exists"
    COINCIDENT = "a_th == a_c: growth begins exactly at fracture"
    NO_ACTIVE_INTERVAL = "a_th > a_c: no growth interval before fracture"
    UNDEFINED = "one or both boundaries do not exist"


def compare_threshold_and_fracture_boundaries(
    threshold_result: FiniteWidthThresholdResult,
    critical_result: FiniteWidthCriticalCrackResult,
) -> BoundaryOrdering:
    """Diagnostic comparison of ``a_th`` and ``a_c``.

    When ``a_th > a_c`` there is no crack size at which this cycle both drives
    growth and remains below the fracture boundary. Under the hard-cutoff model
    a crack below ``a_c`` then stays arrested indefinitely: it cannot grow into
    the fracture boundary under the same constant amplitude. The model does not
    invent growth to bridge that gap.
    """
    if not threshold_result.has_finite_boundary or not critical_result.has_finite_boundary:
        return BoundaryOrdering.UNDEFINED
    a_th = threshold_result.threshold_crack_length
    a_c = critical_result.critical_crack_length
    if a_th < a_c:
        return BoundaryOrdering.ACTIVE_INTERVAL_EXISTS
    if a_th == a_c:
        return BoundaryOrdering.COINCIDENT
    return BoundaryOrdering.NO_ACTIVE_INTERVAL


def classify_growth_state(
    crack_length: float,
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
    threshold: CrackGrowthThreshold | None,
    critical_result: FiniteWidthCriticalCrackResult,
) -> GrowthState:
    """Classify one crack against the threshold and the fracture boundary.

    The fracture boundary is checked FIRST: a crack already at or beyond it is
    reported as such regardless of its threshold state, because no
    crack-growth life remains either way.
    """
    a = require_positive(crack_length, "crack_length")
    if critical_result.status is CriticalCrackStatus.NO_TENSILE_BOUNDARY:
        return GrowthState.NO_TENSILE_FRACTURE_BOUNDARY
    if not critical_result.has_finite_boundary:
        return GrowthState.ABOVE_FRACTURE_BOUNDARY
    a_c = critical_result.critical_crack_length
    if a == a_c:
        return GrowthState.AT_FRACTURE_BOUNDARY
    if a > a_c:
        return GrowthState.ABOVE_FRACTURE_BOUNDARY
    if threshold is None:
        return GrowthState.ACTIVE_GROWTH
    delta_k = delta_k_for_geometry(a, cycle, geometry)
    if delta_k == threshold.delta_k_threshold:
        return GrowthState.AT_THRESHOLD
    if delta_k < threshold.delta_k_threshold:
        return GrowthState.ARRESTED_BELOW_THRESHOLD
    return GrowthState.ACTIVE_GROWTH


@dataclass(frozen=True)
class ThresholdedLifeResult:
    """Threshold-aware life result with an unambiguous growth state."""

    initial_crack_length: float
    plate_width: float
    state: GrowthState
    """The unambiguous verdict for the initial crack."""

    predicted_cycles: float
    """``math.inf`` for an arrested crack or a cycle with no fracture boundary;
    ``0.0`` at or beyond the fracture boundary; otherwise the Milestone 3 life."""

    initial_delta_k: float
    delta_k_threshold: float | None
    """``None`` when the threshold screen is disabled."""

    threshold_ratio: float | None
    """``delta_K(a0)/delta_K_th``; ``None`` when disabled."""

    threshold_excess_ratio: float | None
    """``delta_K(a0)/delta_K_th - 1``; POSITIVE means growing, not safe."""

    critical: FiniteWidthCriticalCrackResult
    threshold_boundary: FiniteWidthThresholdResult | None
    boundary_ordering: BoundaryOrdering
    growth_result: LogGridGrowthResult | None
    """Milestone 3 integration diagnostics, or ``None`` when no integration ran."""

    delta_k_at_critical: float | None

    @property
    def critical_crack_length(self) -> float:
        return self.critical.critical_crack_length

    @property
    def threshold_crack_length(self) -> float | None:
        if self.threshold_boundary is None:
            return None
        return self.threshold_boundary.threshold_crack_length

    @property
    def threshold_to_toughness_ratio(self) -> float | None:
        """``delta_K_th / K_IC`` [-]. A DIAGNOSTIC only, never a criterion."""
        if self.delta_k_threshold is None:
            return None
        return self.delta_k_threshold / self.critical.k_ic


def cycles_to_fracture_with_threshold(
    a_initial: float,
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int = DEFAULT_LOG_INTERVALS,
    **solver_options: float,
) -> ThresholdedLifeResult:
    """Classify the initial crack, then integrate only if it is actually growing.

    Parameters
    ----------
    threshold:
        Pass ``None`` to DISABLE the threshold screen entirely, which reproduces
        the Milestone 3 result exactly. This is the honest way to express
        "no threshold", rather than constructing a zero-valued threshold object
        that the validation would reject.

    Notes
    -----
    An arrested crack returns ``math.inf`` and performs no integration. The
    crack is never advanced to ``a_th`` first -- see the module docstring.
    """
    if not isinstance(geometry, FiniteWidthCenterCrack):
        raise TypeError(
            f"geometry must be a FiniteWidthCenterCrack, got {geometry!r}"
        )
    if not isinstance(paris_law, ParisLaw):
        raise TypeError(f"paris_law must be a ParisLaw, got {paris_law!r}")
    if threshold is not None and not isinstance(threshold, CrackGrowthThreshold):
        raise TypeError(
            "threshold must be a CrackGrowthThreshold or None, got "
            f"{threshold!r}"
        )
    a0 = require_positive(a_initial, "a_initial")
    geometry.geometry_factor_at(a0)  # rejects a0 >= W/2

    critical = finite_width_critical_crack_length(
        cycle, geometry, toughness, **solver_options
    )
    initial_delta_k = delta_k_for_geometry(a0, cycle, geometry)
    state = classify_growth_state(a0, cycle, geometry, threshold, critical)

    threshold_boundary = (
        None
        if threshold is None
        else finite_width_threshold_crack_length(
            cycle, geometry, threshold, **solver_options
        )
    )
    ordering = (
        BoundaryOrdering.UNDEFINED
        if threshold_boundary is None
        else compare_threshold_and_fracture_boundaries(threshold_boundary, critical)
    )
    delta_k_at_critical = (
        delta_k_for_geometry(critical.critical_crack_length, cycle, geometry)
        if critical.has_finite_boundary
        else None
    )

    def _result(
        predicted_cycles: float, growth_result: LogGridGrowthResult | None
    ) -> ThresholdedLifeResult:
        return ThresholdedLifeResult(
            initial_crack_length=a0,
            plate_width=geometry.plate_width,
            state=state,
            predicted_cycles=predicted_cycles,
            initial_delta_k=initial_delta_k,
            delta_k_threshold=(
                None if threshold is None else threshold.delta_k_threshold
            ),
            threshold_ratio=(
                None
                if threshold is None
                else threshold.threshold_ratio(initial_delta_k)
            ),
            threshold_excess_ratio=(
                None
                if threshold is None
                else threshold.threshold_excess_ratio(initial_delta_k)
            ),
            critical=critical,
            threshold_boundary=threshold_boundary,
            boundary_ordering=ordering,
            growth_result=growth_result,
            delta_k_at_critical=delta_k_at_critical,
        )

    if state is GrowthState.NO_TENSILE_FRACTURE_BOUNDARY:
        # delta_K may still be algebraically above threshold here, but with no
        # tensile fracture boundary there is no target to grow to, so no finite
        # cycles-to-fracture is invented.
        return _result(math.inf, None)

    if state in (
        GrowthState.AT_FRACTURE_BOUNDARY,
        GrowthState.ABOVE_FRACTURE_BOUNDARY,
    ):
        return _result(0.0, None)

    if state.is_arrested:
        # No integration: da/dN is exactly zero, so the crack never advances.
        return _result(math.inf, None)

    milestone3 = cycles_to_finite_width_fracture(
        a0, cycle, geometry, paris_law, toughness, intervals=intervals, **solver_options
    )
    return _result(milestone3.predicted_cycles, milestone3.growth_result)
