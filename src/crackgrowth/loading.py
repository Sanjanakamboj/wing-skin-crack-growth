"""Constant-amplitude cyclic stress representation.

Sign convention
---------------
* Tensile stress is POSITIVE.
* Compressive stress is NEGATIVE.
* Mode-I crack opening is driven by tensile (positive) stress.

Milestone 1 models a single constant-amplitude cycle defined by its maximum and
minimum remote stress. A compressive ``sigma_min`` is accepted and is *not*
clamped to zero at the cycle-object level: whether compression contributes to
the effective driving force is a crack-closure modelling decision, and closure
is explicitly out of Milestone 1 scope. See :mod:`crackgrowth.stress_intensity`.
"""

from __future__ import annotations

from dataclasses import dataclass

from ._validation import require_finite

__all__ = ["StressCycle"]


@dataclass(frozen=True)
class StressCycle:
    """A constant-amplitude remote stress cycle, in pascals.

    Parameters
    ----------
    sigma_max:
        Maximum remote stress [Pa]. Tension positive.
    sigma_min:
        Minimum remote stress [Pa]. Tension positive, compression negative.
        Must satisfy ``sigma_min <= sigma_max``.
    """

    sigma_max: float
    sigma_min: float

    def __post_init__(self) -> None:
        sigma_max = require_finite(self.sigma_max, "sigma_max")
        sigma_min = require_finite(self.sigma_min, "sigma_min")
        if sigma_min > sigma_max:
            raise ValueError(
                "sigma_min must not exceed sigma_max "
                f"(got sigma_min={sigma_min!r}, sigma_max={sigma_max!r})"
            )
        object.__setattr__(self, "sigma_max", sigma_max)
        object.__setattr__(self, "sigma_min", sigma_min)

    @property
    def stress_range(self) -> float:
        """Algebraic stress range ``delta_sigma = sigma_max - sigma_min`` [Pa].

        This is the raw algebraic range. A compressive ``sigma_min`` therefore
        increases it. Milestone 1 does not model crack closure, so no effective
        (closure-corrected) range is offered.
        """
        return self.sigma_max - self.sigma_min

    @property
    def mean_stress(self) -> float:
        """Mean stress ``(sigma_max + sigma_min) / 2`` [Pa]."""
        return 0.5 * (self.sigma_max + self.sigma_min)

    @property
    def alternating_stress(self) -> float:
        """Alternating stress ``delta_sigma / 2`` [Pa]."""
        return 0.5 * self.stress_range

    @property
    def stress_ratio(self) -> float:
        """Stress ratio ``R = sigma_min / sigma_max`` [-].

        Raises
        ------
        ZeroDivisionError
            If ``sigma_max == 0``, for which ``R`` is undefined. Milestone 1
            deliberately raises rather than substituting a placeholder value.
        """
        if self.sigma_max == 0.0:
            raise ZeroDivisionError(
                "stress ratio R is undefined when sigma_max is zero"
            )
        return self.sigma_min / self.sigma_max
