"""Sensitivity sweeps around the toughness-derived fracture boundary.

Each sweep varies exactly ONE quantity and states what is held constant. Two of
them are deliberately distinct and must not be conflated:

* :func:`max_stress_sensitivity_at_fixed_range` varies ``sigma_max`` while
  holding ``delta_sigma`` FIXED (so ``sigma_min = sigma_max - delta_sigma`` and
  ``R`` changes). Paris growth at a given crack length is therefore UNCHANGED,
  and only the fracture endpoint moves. This isolates ``K_max``.
* :func:`stress_range_sensitivity_at_fixed_min` varies ``delta_sigma`` while
  holding ``sigma_min`` FIXED (so ``sigma_max`` and hence the fracture boundary
  move too). This changes both driving forces at once.

As in Milestone 1 these helpers add no physics: every row is a complete,
independently inspectable analysis.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ._validation import require_positive
from .fracture import (
    FractureToughness,
    assess_fracture,
    critical_crack_length,
    residual_strength,
)
from .fracture_life import cycles_to_critical_crack
from .geometry import ThroughCrackGeometry
from .integration import DEFAULT_INTEGRATION_INTERVALS
from .loading import StressCycle
from .paris import ParisLaw, crack_growth_rate
from .stress_intensity import delta_stress_intensity

__all__ = [
    "FracturePoint",
    "ResidualStrengthRow",
    "toughness_sensitivity",
    "max_stress_sensitivity_at_fixed_range",
    "stress_range_sensitivity_at_fixed_min",
    "geometry_factor_sensitivity",
    "initial_crack_sensitivity_to_fracture",
    "residual_strength_table",
]


@dataclass(frozen=True)
class FracturePoint:
    """One swept case: the varied quantity, the boundary, and the life."""

    parameter_name: str
    parameter_value: float
    initial_crack_length: float
    critical_crack_length: float
    sigma_max: float
    stress_range: float
    geometry_factor: float
    k_ic: float
    initial_delta_k: float
    initial_growth_rate: float
    initial_utilization: float
    initial_margin: float
    numerical_cycles: float
    analytical_cycles: float

    @property
    def initial_crack_fraction(self) -> float:
        """``a0 / a_critical`` [-]."""
        return self.initial_crack_length / self.critical_crack_length


def _point(
    parameter_name: str,
    parameter_value: float,
    a_initial: float,
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int,
) -> FracturePoint:
    result = cycles_to_critical_crack(
        a_initial, cycle, geometry, paris_law, toughness, intervals=intervals
    )
    return FracturePoint(
        parameter_name=parameter_name,
        parameter_value=float(parameter_value),
        initial_crack_length=a_initial,
        critical_crack_length=result.critical_crack_length,
        sigma_max=cycle.sigma_max,
        stress_range=cycle.stress_range,
        geometry_factor=geometry.geometry_factor,
        k_ic=toughness.k_ic,
        initial_delta_k=delta_stress_intensity(a_initial, cycle, geometry),
        initial_growth_rate=crack_growth_rate(
            a_initial, cycle, geometry, paris_law
        ),
        initial_utilization=result.initial_assessment.utilization,
        initial_margin=result.initial_assessment.margin,
        numerical_cycles=result.predicted_cycles,
        analytical_cycles=result.analytical_cycles,
    )


def toughness_sensitivity(
    k_ic_values: Sequence[float],
    a_initial: float,
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_INTEGRATION_INTERVALS,
) -> tuple[FracturePoint, ...]:
    """Vary ``K_IC`` [Pa*sqrt(m)]; everything else fixed.

    The provenance note of ``toughness`` is carried onto each swept record, so a
    sweep of illustrative values stays flagged as illustrative.
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
                "k_ic", k_ic, a_initial, cycle, geometry, paris_law, swept, intervals
            )
        )
    return tuple(points)


def max_stress_sensitivity_at_fixed_range(
    sigma_max_values: Sequence[float],
    stress_range: float,
    a_initial: float,
    geometry: ThroughCrackGeometry,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_INTEGRATION_INTERVALS,
) -> tuple[FracturePoint, ...]:
    """Vary ``sigma_max`` at FIXED ``delta_sigma`` (``sigma_min`` follows, ``R`` changes).

    Because ``delta_sigma`` is held constant, ``delta_K(a)`` and therefore the
    Paris growth rate at any given crack length are unchanged across the sweep.
    Only the fracture boundary ``a_c`` moves -- inward as ``sigma_max`` rises --
    so the life falls purely because the endpoint moves, not because the crack
    grows faster. This isolates ``K_max`` from ``delta_K``.
    """
    delta_sigma = require_positive(stress_range, "stress_range")
    points = []
    for sigma_max in sigma_max_values:
        cycle = StressCycle(
            sigma_max=sigma_max, sigma_min=sigma_max - delta_sigma
        )
        points.append(
            _point(
                "sigma_max",
                sigma_max,
                a_initial,
                cycle,
                geometry,
                paris_law,
                toughness,
                intervals,
            )
        )
    return tuple(points)


def stress_range_sensitivity_at_fixed_min(
    stress_ranges: Sequence[float],
    sigma_min: float,
    a_initial: float,
    geometry: ThroughCrackGeometry,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_INTEGRATION_INTERVALS,
) -> tuple[FracturePoint, ...]:
    """Vary ``delta_sigma`` at FIXED ``sigma_min`` (so ``sigma_max`` moves too).

    Distinct from :func:`max_stress_sensitivity_at_fixed_range`: here BOTH the
    growth driving force ``delta_K`` and the fracture boundary ``a_c`` change,
    so the life trend is the combination of the two effects rather than an
    isolated endpoint effect.
    """
    points = []
    for delta_sigma in stress_ranges:
        cycle = StressCycle(
            sigma_max=sigma_min + delta_sigma, sigma_min=sigma_min
        )
        points.append(
            _point(
                "stress_range",
                delta_sigma,
                a_initial,
                cycle,
                geometry,
                paris_law,
                toughness,
                intervals,
            )
        )
    return tuple(points)


def geometry_factor_sensitivity(
    geometry_factors: Sequence[float],
    a_initial: float,
    cycle: StressCycle,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_INTEGRATION_INTERVALS,
) -> tuple[FracturePoint, ...]:
    """Vary ``Y``; the cycle, flaw size, Paris curve and toughness are fixed.

    ``Y`` raises the growth rate (as ``Y**m``) AND pulls the fracture boundary in
    (as ``Y**-2``), so the life falls on both counts.
    """
    points = []
    for y in geometry_factors:
        points.append(
            _point(
                "geometry_factor",
                y,
                a_initial,
                cycle,
                ThroughCrackGeometry(geometry_factor=y),
                paris_law,
                toughness,
                intervals,
            )
        )
    return tuple(points)


def initial_crack_sensitivity_to_fracture(
    initial_crack_lengths: Sequence[float],
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_INTEGRATION_INTERVALS,
) -> tuple[FracturePoint, ...]:
    """Vary the assumed initial flaw size against a fixed fracture boundary."""
    return tuple(
        _point(
            "initial_crack_length",
            a0,
            a0,
            cycle,
            geometry,
            paris_law,
            toughness,
            intervals,
        )
        for a0 in initial_crack_lengths
    )


@dataclass(frozen=True)
class ResidualStrengthRow:
    """One crack length screened against the fracture boundary."""

    crack_length: float
    k_max: float
    k_ic: float
    residual_strength: float
    applied_sigma_max: float
    utilization: float
    margin: float
    passes: bool


def residual_strength_table(
    crack_lengths: Sequence[float],
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    toughness: FractureToughness,
) -> tuple[ResidualStrengthRow, ...]:
    """Screen a list of crack lengths; no new mechanics, only reporting."""
    rows = []
    for a in crack_lengths:
        assessment = assess_fracture(a, cycle, geometry, toughness)
        rows.append(
            ResidualStrengthRow(
                crack_length=assessment.crack_length,
                k_max=assessment.k_max,
                k_ic=assessment.k_ic,
                residual_strength=residual_strength(a, geometry, toughness),
                applied_sigma_max=assessment.applied_sigma_max,
                utilization=assessment.utilization,
                margin=assessment.margin,
                passes=assessment.passes,
            )
        )
    return tuple(rows)


# Re-exported for callers that want the boundary without the sweep machinery.
_ = critical_crack_length
