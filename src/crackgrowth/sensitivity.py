"""Deterministic single-parameter sensitivity sweeps.

These helpers only re-run the Milestone 1 life calculation over a list of
explicitly supplied parameter values. They introduce no new physics, apply no
correction, and hide no assumption: each swept case is a complete, independently
inspectable analysis in its own right.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from ._validation import require_positive
from .geometry import ThroughCrackGeometry
from .integration import (
    DEFAULT_INTEGRATION_INTERVALS,
    analytical_cycles_to_crack_length,
    cycles_to_crack_length,
)
from .loading import StressCycle
from .paris import ParisLaw

__all__ = [
    "SensitivityPoint",
    "stress_range_sensitivity",
    "initial_crack_length_sensitivity",
]


@dataclass(frozen=True)
class SensitivityPoint:
    """One swept case: the varied quantity and the resulting life."""

    parameter_name: str
    parameter_value: float
    initial_crack_length: float
    final_crack_length: float
    stress_range: float
    numerical_cycles: float
    analytical_cycles: float


def stress_range_sensitivity(
    scales: Sequence[float],
    a_initial: float,
    a_final: float,
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    paris_law: ParisLaw,
    intervals: int = DEFAULT_INTEGRATION_INTERVALS,
) -> tuple[SensitivityPoint, ...]:
    """Scale the stress cycle by each factor in ``scales`` and report the life.

    Both ``sigma_max`` and ``sigma_min`` are multiplied by the scale factor, so
    the stress ratio ``R`` is preserved and only ``delta_sigma`` changes. For the
    constant-Y Paris model this must give ``N`` proportional to
    ``delta_sigma ** -m``.
    """
    points = []
    for scale in scales:
        s = require_positive(scale, "scale")
        scaled = StressCycle(
            sigma_max=cycle.sigma_max * s, sigma_min=cycle.sigma_min * s
        )
        result = cycles_to_crack_length(
            a_initial, a_final, scaled, geometry, paris_law, intervals=intervals
        )
        points.append(
            SensitivityPoint(
                parameter_name="stress_range_scale",
                parameter_value=s,
                initial_crack_length=a_initial,
                final_crack_length=a_final,
                stress_range=scaled.stress_range,
                numerical_cycles=result.predicted_cycles,
                analytical_cycles=analytical_cycles_to_crack_length(
                    a_initial, a_final, scaled, geometry, paris_law
                ),
            )
        )
    return tuple(points)


def initial_crack_length_sensitivity(
    initial_crack_lengths: Sequence[float],
    a_final: float,
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    paris_law: ParisLaw,
    intervals: int = DEFAULT_INTEGRATION_INTERVALS,
) -> tuple[SensitivityPoint, ...]:
    """Vary the initial crack length at a fixed target and report the life."""
    points = []
    for a0 in initial_crack_lengths:
        result = cycles_to_crack_length(
            a0, a_final, cycle, geometry, paris_law, intervals=intervals
        )
        points.append(
            SensitivityPoint(
                parameter_name="initial_crack_length",
                parameter_value=float(a0),
                initial_crack_length=float(a0),
                final_crack_length=a_final,
                stress_range=cycle.stress_range,
                numerical_cycles=result.predicted_cycles,
                analytical_cycles=analytical_cycles_to_crack_length(
                    a0, a_final, cycle, geometry, paris_law
                ),
            )
        )
    return tuple(points)
