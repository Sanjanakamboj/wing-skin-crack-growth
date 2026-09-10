"""Sensitivity sweeps for the finite-width centre-cracked panel.

Each sweep varies exactly one quantity and states what is held fixed. As in
earlier milestones these helpers add no physics: every row is a complete,
independently inspectable analysis.

Two stress sweeps are deliberately distinct and must not be conflated:

* :func:`max_stress_sensitivity_at_fixed_range` holds ``delta_sigma`` fixed, so
  ``delta_K(a)`` and the growth rate at any crack length are UNCHANGED and only
  the fracture boundary moves;
* :func:`stress_range_sensitivity_at_fixed_min` holds ``sigma_min`` fixed, so
  both the growth driving force and the boundary move together.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ._validation import require_positive
from .finite_width import (
    FiniteWidthCenterCrack,
    delta_k_for_geometry,
    ligament_fraction,
)
from .finite_width_fracture import k_max_for_geometry
from .finite_width_life import cycles_to_finite_width_fracture
from .fracture import FractureToughness
from .integration_log import DEFAULT_LOG_INTERVALS
from .loading import StressCycle
from .paris import ParisLaw, crack_growth_rate

__all__ = [
    "FiniteWidthPoint",
    "GeometryAmplificationRow",
    "width_sensitivity",
    "initial_crack_sensitivity_finite_width",
    "max_stress_sensitivity_at_fixed_range",
    "stress_range_sensitivity_at_fixed_min",
    "toughness_sensitivity_finite_width",
    "geometry_amplification_table",
]


@dataclass(frozen=True)
class FiniteWidthPoint:
    """One swept finite-width case."""

    parameter_name: str
    parameter_value: float
    plate_width: float
    initial_crack_length: float
    initial_geometry_factor: float
    critical_crack_length: float
    geometry_factor_at_critical: float
    ligament_fraction_at_critical: float
    sigma_max: float
    stress_range: float
    k_ic: float
    initial_delta_k: float
    initial_growth_rate: float
    initial_utilization: float
    predicted_cycles: float

    @property
    def initial_crack_width_ratio(self) -> float:
        """``a0 / W`` [-]."""
        return self.initial_crack_length / self.plate_width


def _point(
    parameter_name: str,
    parameter_value: float,
    a_initial: float,
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int,
) -> FiniteWidthPoint:
    result = cycles_to_finite_width_fracture(
        a_initial, cycle, geometry, paris_law, toughness, intervals=intervals
    )
    growth = result.growth_result
    return FiniteWidthPoint(
        parameter_name=parameter_name,
        parameter_value=float(parameter_value),
        plate_width=geometry.plate_width,
        initial_crack_length=a_initial,
        initial_geometry_factor=result.initial_geometry_factor,
        critical_crack_length=result.critical_crack_length,
        geometry_factor_at_critical=result.critical.geometry_factor_at_critical,
        ligament_fraction_at_critical=result.critical.ligament_fraction,
        sigma_max=cycle.sigma_max,
        stress_range=cycle.stress_range,
        k_ic=toughness.k_ic,
        initial_delta_k=delta_k_for_geometry(a_initial, cycle, geometry),
        initial_growth_rate=(
            growth.initial_growth_rate if growth is not None else 0.0
        ),
        initial_utilization=result.initial_assessment.utilization,
        predicted_cycles=result.predicted_cycles,
    )


def width_sensitivity(
    plate_widths: Sequence[float],
    a_initial: float,
    cycle: StressCycle,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_LOG_INTERVALS,
) -> tuple[FiniteWidthPoint, ...]:
    """Vary the panel width ``W``; everything else fixed.

    Requires ``2*a_initial < W`` for every width in the sweep. Wider panels give
    a lower ``Y`` at any fixed crack length, a larger critical crack, and a
    longer life, converging on the infinite-plate result.
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
                intervals,
            )
        )
    return tuple(points)


def initial_crack_sensitivity_finite_width(
    initial_crack_lengths: Sequence[float],
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_LOG_INTERVALS,
) -> tuple[FiniteWidthPoint, ...]:
    """Vary the assumed initial flaw against a fixed panel and boundary."""
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


def max_stress_sensitivity_at_fixed_range(
    sigma_max_values: Sequence[float],
    stress_range: float,
    a_initial: float,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_LOG_INTERVALS,
) -> tuple[FiniteWidthPoint, ...]:
    """Vary ``sigma_max`` at FIXED ``delta_sigma`` (``sigma_min`` follows).

    ``delta_K(a)`` is identical across the sweep, so the growth rate at any
    crack length is unchanged and the life falls purely because the fracture
    boundary moves inward -- the same ``K_max`` / ``delta_K`` separation
    demonstrated in Milestone 2, now with a crack-size-dependent ``Y``.
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
            intervals,
        )
        for sigma_max in sigma_max_values
    )


def stress_range_sensitivity_at_fixed_min(
    sigma_max_values: Sequence[float],
    sigma_min: float,
    a_initial: float,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_LOG_INTERVALS,
) -> tuple[FiniteWidthPoint, ...]:
    """Vary ``sigma_max`` at FIXED ``sigma_min``, so ``delta_sigma`` moves too.

    Distinct from :func:`max_stress_sensitivity_at_fixed_range`: here BOTH the
    growth driving force and the fracture boundary change.
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
            intervals,
        )
        for sigma_max in sigma_max_values
    )


def toughness_sensitivity_finite_width(
    k_ic_values: Sequence[float],
    a_initial: float,
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_LOG_INTERVALS,
) -> tuple[FiniteWidthPoint, ...]:
    """Vary ``K_IC``; everything else fixed.

    Note that the Milestone 2 exact scaling ``a_c ~ K_IC**2`` does NOT survive
    finite width: ``Y(a_c)`` grows with ``a_c``, so a larger toughness buys
    proportionally less critical crack length. This is verified by test.
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
                intervals,
            )
        )
    return tuple(points)


@dataclass(frozen=True)
class GeometryAmplificationRow:
    """How strongly finite width amplifies K at one crack length."""

    crack_length: float
    plate_width: float
    crack_width_ratio: float
    """``a / W`` [-]."""
    total_crack_length: float
    ligament_fraction: float
    geometry_factor: float
    k_max: float
    delta_k: float
    k_max_infinite_plate: float
    """``K_max`` for the same crack with ``Y = 1``, for comparison."""

    @property
    def amplification(self) -> float:
        """``K_max / K_max(Y=1)`` [-], which equals ``Y(a)``."""
        return self.k_max / self.k_max_infinite_plate


def geometry_amplification_table(
    crack_lengths: Sequence[float],
    cycle: StressCycle,
    geometry: FiniteWidthCenterCrack,
) -> tuple[GeometryAmplificationRow, ...]:
    """Tabulate ``Y(a)``, ``K_max`` and ``delta_K`` over a range of crack sizes.

    Crack lengths at or beyond ``W/2`` are rejected rather than skipped, so an
    inadmissible request is never silently dropped from the table.
    """
    from .geometry import ThroughCrackGeometry
    from .stress_intensity import stress_intensity

    infinite = ThroughCrackGeometry(geometry_factor=1.0)
    rows = []
    for a in crack_lengths:
        rows.append(
            GeometryAmplificationRow(
                crack_length=a,
                plate_width=geometry.plate_width,
                crack_width_ratio=a / geometry.plate_width,
                total_crack_length=geometry.total_crack_length(a),
                ligament_fraction=ligament_fraction(a, geometry.plate_width),
                geometry_factor=geometry.geometry_factor_at(a),
                k_max=k_max_for_geometry(a, cycle, geometry),
                delta_k=delta_k_for_geometry(a, cycle, geometry),
                k_max_infinite_plate=stress_intensity(
                    cycle.sigma_max, a, infinite
                ),
            )
        )
    return tuple(rows)


# Kept importable for callers comparing against the constant-Y growth rate.
_ = crack_growth_rate
