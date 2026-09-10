"""Paris-law material representation and crack-growth rate.

Governing relation
------------------
    da/dN = C * (delta_K) ** m

Units
-----
This package stores and evaluates the Paris coefficient ``C`` on an SI basis,
i.e. with ``delta_K`` expressed in Pa*sqrt(m) and ``da/dN`` in m/cycle. The
exponent ``m`` is dimensionless, so ``C`` carries units of
``m/cycle / (Pa*sqrt(m))**m`` and its numerical value therefore depends on ``m``.

Handbook Paris coefficients are almost always quoted on an MPa*sqrt(m) basis:

    da/dN [m/cycle] = C_MPa * (delta_K [MPa*sqrt(m)]) ** m

Substituting ``delta_K[MPa*sqrt(m)] = delta_K[Pa*sqrt(m)] / 1e6`` gives

    da/dN = C_MPa * (delta_K[Pa*sqrt(m)] / 1e6) ** m
          = (C_MPa / (1e6) ** m) * (delta_K[Pa*sqrt(m)]) ** m

hence the exact conversion

    C_SI = C_MPa / (1e6) ** m

which is implemented by :func:`paris_c_from_mpa_basis`. Mixing the two bases
silently is a classic and expensive error, so this package never guesses: a
:class:`ParisLaw` value is always interpreted on the SI basis, and conversion
from the MPa basis must be requested explicitly.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ._validation import require_finite, require_positive
from .geometry import ThroughCrackGeometry
from .loading import StressCycle
from .stress_intensity import delta_stress_intensity

__all__ = [
    "ParisLaw",
    "CrackGrowthRatePoint",
    "MPA_ROOT_M_IN_PA_ROOT_M",
    "paris_c_from_mpa_basis",
    "paris_c_to_mpa_basis",
    "paris_law_from_mpa_basis",
    "crack_growth_rate",
    "crack_growth_rate_point",
    "ILLUSTRATIVE_DISCLAIMER",
    "ILLUSTRATIVE_ALUMINIUM_LIKE_PARIS",
]

#: 1 MPa*sqrt(m) expressed in Pa*sqrt(m).
MPA_ROOT_M_IN_PA_ROOT_M = 1.0e6

#: Mandatory provenance banner for any curve that is not a verified dataset.
ILLUSTRATIVE_DISCLAIMER = "ILLUSTRATIVE PARIS-LAW INPUT - NOT DESIGN ALLOWABLE"


@dataclass(frozen=True)
class ParisLaw:
    """Immutable Paris-law material record on the SI (Pa*sqrt(m)) basis.

    Parameters
    ----------
    C:
        Paris coefficient in ``m/cycle / (Pa*sqrt(m))**m``. Finite, > 0.
    m:
        Paris exponent [-]. Finite, > 0.
    name:
        Short identifier for the curve.
    source_note:
        Provenance statement. For any curve that is not traceable to a verified
        dataset this must carry :data:`ILLUSTRATIVE_DISCLAIMER`.
    notes:
        Optional additional commentary (validity range, caveats).
    """

    C: float
    m: float
    name: str
    source_note: str
    notes: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "C", require_positive(self.C, "C"))
        object.__setattr__(self, "m", require_positive(self.m, "m"))
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("name must be a non-empty string")
        if not isinstance(self.source_note, str) or not self.source_note.strip():
            raise ValueError("source_note must be a non-empty provenance statement")
        if not isinstance(self.notes, str):
            raise TypeError("notes must be a string")

    @property
    def c_on_mpa_basis(self) -> float:
        """This curve's coefficient re-expressed on the MPa*sqrt(m) basis."""
        return paris_c_to_mpa_basis(self.C, self.m)

    def growth_rate(self, delta_k: float) -> float:
        """Evaluate ``da/dN = C * delta_K**m`` for a given range [Pa*sqrt(m)].

        ``delta_k`` must be finite and non-negative. A zero range returns
        exactly zero growth (see :func:`crack_growth_rate`).
        """
        dk = require_finite(delta_k, "delta_k")
        if dk < 0.0:
            raise ValueError(f"delta_k must be non-negative, got {dk!r}")
        if dk == 0.0:
            return 0.0
        return self.C * dk**self.m


def paris_c_from_mpa_basis(c_mpa: float, m: float) -> float:
    """Convert a Paris coefficient from the MPa*sqrt(m) basis to the SI basis.

    ``C_SI = C_MPa / (1e6) ** m`` -- see the module docstring for the derivation.
    """
    c = require_positive(c_mpa, "c_mpa")
    exponent = require_positive(m, "m")
    return c / MPA_ROOT_M_IN_PA_ROOT_M**exponent


def paris_c_to_mpa_basis(c_si: float, m: float) -> float:
    """Inverse of :func:`paris_c_from_mpa_basis`: ``C_MPa = C_SI * (1e6) ** m``."""
    c = require_positive(c_si, "c_si")
    exponent = require_positive(m, "m")
    return c * MPA_ROOT_M_IN_PA_ROOT_M**exponent


def paris_law_from_mpa_basis(
    c_mpa: float,
    m: float,
    name: str,
    source_note: str,
    notes: str = "",
) -> ParisLaw:
    """Build a :class:`ParisLaw` from a handbook-style MPa*sqrt(m) coefficient.

    The conversion is explicit and named; nothing about the units is inferred.
    """
    return ParisLaw(
        C=paris_c_from_mpa_basis(c_mpa, m),
        m=m,
        name=name,
        source_note=source_note,
        notes=notes,
    )


@dataclass(frozen=True)
class CrackGrowthRatePoint:
    """Auditable single-point crack-growth-rate evaluation."""

    crack_length: float
    delta_stress: float
    geometry_factor: float
    delta_k: float
    growth_rate: float


def crack_growth_rate(
    crack_length: float,
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    paris_law: ParisLaw,
) -> float:
    """Paris-law crack-growth rate ``da/dN`` [m/cycle] at crack length ``a``.

    Evaluates ``delta_K = Y * delta_sigma * sqrt(pi * a)`` and then
    ``da/dN = C * delta_K ** m``.

    Zero-range policy: if ``delta_sigma == 0`` then ``delta_K == 0`` and the
    growth rate is exactly ``0.0``. No threshold ``delta_K_th`` and no
    fracture-toughness cutoff is applied, so every non-zero range produces some
    growth, however small -- that is a Milestone 1 modelling limitation, not a
    physical claim.
    """
    if not isinstance(paris_law, ParisLaw):
        raise TypeError(f"paris_law must be a ParisLaw, got {paris_law!r}")
    delta_k = delta_stress_intensity(crack_length, cycle, geometry)
    return paris_law.growth_rate(delta_k)


def crack_growth_rate_point(
    crack_length: float,
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    paris_law: ParisLaw,
) -> CrackGrowthRatePoint:
    """As :func:`crack_growth_rate`, returning the intermediate values too."""
    if not isinstance(paris_law, ParisLaw):
        raise TypeError(f"paris_law must be a ParisLaw, got {paris_law!r}")
    delta_k = delta_stress_intensity(crack_length, cycle, geometry)
    return CrackGrowthRatePoint(
        crack_length=float(crack_length),
        delta_stress=cycle.stress_range,
        geometry_factor=geometry.geometry_factor,
        delta_k=delta_k,
        growth_rate=paris_law.growth_rate(delta_k),
    )


#: Canonical Milestone 1 curve.
#:
#: These coefficients are ILLUSTRATIVE. They were chosen after auditing the
#: resulting delta_K range (5.6 - 17.7 MPa*sqrt(m)), growth rates
#: (1.8e-9 - 5.6e-8 m/cycle) and life (~7.8e5 cycles) for the canonical wing-skin
#: case, so that the model sits in a numerically useful and physically plausible
#: mid-Paris regime. They are NOT traceable to a measured dataset and no specific
#: aerospace alloy is claimed.
ILLUSTRATIVE_ALUMINIUM_LIKE_PARIS = paris_law_from_mpa_basis(
    c_mpa=1.0e-11,
    m=3.0,
    name="Illustrative aluminium-like Paris curve",
    source_note=(
        ILLUSTRATIVE_DISCLAIMER
        + ": C = 1.0e-11 m/cycle on the MPa*sqrt(m) basis with m = 3.0, chosen to "
        "place the canonical wing-skin case in a plausible mid-Paris regime. Not "
        "measured, not traceable to a qualification dataset, and no specific alloy "
        "is claimed."
    ),
    notes=(
        "Round illustrative values. Audited over a = 1-10 mm at delta_sigma = 100 "
        "MPa: delta_K = 5.6-17.7 MPa*sqrt(m), da/dN = 1.8e-9 to 5.6e-8 m/cycle. "
        "Coefficients were not tuned to produce any particular target life."
    ),
)
