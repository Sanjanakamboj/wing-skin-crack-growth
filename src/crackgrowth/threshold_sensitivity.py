"""Threshold-aware sensitivity sweeps.

Each sweep varies one quantity and states what is held fixed. No new physics:
every row is a complete, independently inspectable threshold-aware analysis.

A ``None`` entry in a threshold sweep means the screen is DISABLED for that row
(reproducing the Milestone 3 result), which is the honest representation of
"zero threshold" -- a zero-valued threshold object would be rejected by
validation, and rightly so.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ._validation import require_positive
from .finite_width import FiniteWidthCenterCrack
from .fracture import FractureToughness
from .integration_log import DEFAULT_LOG_INTERVALS
from .loading import StressCycle
from .paris import ParisLaw
from .threshold import CrackGrowthThreshold
from .threshold_life import (
    BoundaryOrdering,
    GrowthState,
    cycles_to_fracture_with_threshold,
)

__all__ = [
    "ThresholdPoint",
    "threshold_sensitivity",
    "initial_crack_sensitivity_with_threshold",
    "stress_range_sensitivity_with_threshold",
    "max_stress_sensitivity_with_threshold",
    "width_sensitivity_with_threshold",
    "toughness_sensitivity_with_threshold",
]


@dataclass(frozen=True)
class ThresholdPoint:
    """One swept threshold-aware case."""

    parameter_name: str
    parameter_value: float
    plate_width: float
    initial_crack_length: float
    initial_geometry_factor: float
    initial_delta_k: float
    delta_k_threshold: float | None
    threshold_ratio: float | None
    state: GrowthState
    threshold_crack_length: float | None
    critical_crack_length: float
    boundary_ordering: BoundaryOrdering
    sigma_max: float
    stress_range: float
    k_ic: float
    predicted_cycles: float

    @property
    def active_growth(self) -> bool:
        return self.state.is_growing

    @property
    def arrested(self) -> bool:
        return self.state.is_arrested


def _point(
    parameter_name: str,
    parameter_value: float,
    a_initial: float,
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int,
) -> ThresholdPoint:
    result = cycles_to_fracture_with_threshold(
        a_initial,
        cycle,
        geometry,
        paris_law,
        toughness,
        threshold,
        intervals=intervals,
    )
    return ThresholdPoint(
        parameter_name=parameter_name,
        parameter_value=float(parameter_value),
        plate_width=geometry.plate_width,
        initial_crack_length=a_initial,
        initial_geometry_factor=geometry.geometry_factor_at(a_initial),
        initial_delta_k=result.initial_delta_k,
        delta_k_threshold=result.delta_k_threshold,
        threshold_ratio=result.threshold_ratio,
        state=result.state,
        threshold_crack_length=result.threshold_crack_length,
        critical_crack_length=result.critical_crack_length,
        boundary_ordering=result.boundary_ordering,
        sigma_max=cycle.sigma_max,
        stress_range=cycle.stress_range,
        k_ic=toughness.k_ic,
        predicted_cycles=result.predicted_cycles,
    )


def threshold_sensitivity(
    thresholds: Sequence[CrackGrowthThreshold | None],
    a_initial: float,
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_LOG_INTERVALS,
) -> tuple[ThresholdPoint, ...]:
    """Vary ``delta_K_th``; everything else fixed.

    A ``None`` entry disables the screen for that row. Because the hard cutoff
    does not modify the Paris rate, every active row returns the SAME life --
    the identical Milestone 3 value -- and every arrested row returns infinity.
    The sweep therefore shows a step, not a trend.
    """
    return tuple(
        _point(
            "delta_k_threshold",
            0.0 if threshold is None else threshold.delta_k_threshold,
            a_initial,
            cycle,
            geometry,
            paris_law,
            toughness,
            threshold,
            intervals,
        )
        for threshold in thresholds
    )


def initial_crack_sensitivity_with_threshold(
    initial_crack_lengths: Sequence[float],
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int = DEFAULT_LOG_INTERVALS,
) -> tuple[ThresholdPoint, ...]:
    """Vary the assumed initial flaw against a fixed panel, cycle and threshold.

    Small flaws may fall below threshold and arrest; larger ones grow. Whether a
    transition appears depends on the inputs and is reported, not engineered.
    """
    return tuple(
        _point(
            "initial_crack_length",
            a0,
            a0,
            cycle,
            geometry,
            paris_law,
            toughness,
            threshold,
            intervals,
        )
        for a0 in initial_crack_lengths
    )


def stress_range_sensitivity_with_threshold(
    sigma_max_values: Sequence[float],
    sigma_min: float,
    a_initial: float,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int = DEFAULT_LOG_INTERVALS,
) -> tuple[ThresholdPoint, ...]:
    """Vary ``sigma_max`` at FIXED ``sigma_min``, so ``delta_sigma`` moves too.

    Both the growth driving force and the fracture boundary change, and low
    ranges may fall below threshold entirely.
    """
    return tuple(
        _point(
            "sigma_max",
            sigma_max,
            a_initial,
            StressCycle(sigma_max=sigma_max, sigma_min=sigma_min),
            geometry,
            paris_law,
            toughness,
            threshold,
            intervals,
        )
        for sigma_max in sigma_max_values
    )


def max_stress_sensitivity_with_threshold(
    sigma_max_values: Sequence[float],
    stress_range: float,
    a_initial: float,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int = DEFAULT_LOG_INTERVALS,
) -> tuple[ThresholdPoint, ...]:
    """Vary ``sigma_max`` at FIXED ``delta_sigma`` (``sigma_min`` follows).

    ``delta_K`` is unchanged across the sweep, so the threshold classification
    is unchanged too -- the threshold sees only ``delta_sigma``. Only the
    fracture endpoint moves, so for active rows the life changes purely through
    ``a_c``. This separates the threshold and fracture criteria cleanly.
    """
    delta_sigma = require_positive(stress_range, "stress_range")
    return tuple(
        _point(
            "sigma_max",
            sigma_max,
            a_initial,
            StressCycle(sigma_max=sigma_max, sigma_min=sigma_max - delta_sigma),
            geometry,
            paris_law,
            toughness,
            threshold,
            intervals,
        )
        for sigma_max in sigma_max_values
    )


def width_sensitivity_with_threshold(
    plate_widths: Sequence[float],
    a_initial: float,
    cycle: StressCycle,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int = DEFAULT_LOG_INTERVALS,
) -> tuple[ThresholdPoint, ...]:
    """Vary the panel width ``W``.

    A narrower panel raises ``Y(a0)`` and therefore ``delta_K(a0)`` slightly, so
    width can in principle move a marginal crack across the threshold. The
    effect is usually small at ``a0 << W`` and is measured, not assumed.
    """
    points = []
    for width in plate_widths:
        w = require_positive(width, "plate_width")
        if 2.0 * a_initial >= w:
            raise ValueError(
                f"initial crack is not admissible at W={w!r} m: require "
                f"2*a0 < W but 2*a0={2.0 * a_initial!r} m (a0 is the HALF "
                "crack length)"
            )
        points.append(
            _point(
                "plate_width",
                w,
                a_initial,
                cycle,
                FiniteWidthCenterCrack(plate_width=w),
                paris_law,
                toughness,
                threshold,
                intervals,
            )
        )
    return tuple(points)


def toughness_sensitivity_with_threshold(
    k_ic_values: Sequence[float],
    a_initial: float,
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int = DEFAULT_LOG_INTERVALS,
) -> tuple[ThresholdPoint, ...]:
    """Vary ``K_IC``; everything else fixed.

    ``K_IC`` does not enter ``delta_K``, so the threshold classification at the
    initial crack is INDEPENDENT of toughness. Only the fracture endpoint -- and
    hence the active-growth life -- changes. This separation is verified by test.
    """
    points = []
    for k_ic in k_ic_values:
        swept = FractureToughness(
            name=f"{toughness.name} (swept K_IC)",
            k_ic=k_ic,
            source_note=toughness.source_note,
            notes=toughness.notes,
        )
        points.append(
            _point(
                "k_ic",
                k_ic,
                a_initial,
                cycle,
                geometry,
                paris_law,
                swept,
                threshold,
                intervals,
            )
        )
    return tuple(points)
