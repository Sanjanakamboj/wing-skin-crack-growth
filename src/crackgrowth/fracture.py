"""Fracture-toughness screening: critical crack size and residual strength.

Milestone 2 replaces Milestone 1's imposed final crack length with a fracture
boundary derived from the material's mode-I plane-strain fracture toughness.

Two distinct driving forces
---------------------------
This module rests on a distinction that must not be blurred:

* **Fatigue crack GROWTH** is driven by the stress-intensity RANGE ``delta_K``
  (Milestone 1). It depends on ``delta_sigma = sigma_max - sigma_min``.
* **FRACTURE** is driven by the MAXIMUM stress intensity in the cycle,
  ``K_max``. It depends on ``sigma_max`` alone.

    K_max(a) = Y * sigma_max * sqrt(pi * a)

    fracture when  K_max = K_IC

Using ``delta_K`` as a fracture criterion is a real and expensive error, and is
guarded against by explicit tests. In general ``delta_K != K_IC`` at fracture;
for constant ``Y`` the two are related by

    delta_K(a) / K_max(a) = delta_sigma / sigma_max

so at the critical crack size ``delta_K(a_c) = (delta_sigma / sigma_max) * K_IC``.

Critical crack size
-------------------
Solving ``K_IC = Y * sigma_max * sqrt(pi * a_c)`` for ``a_c``:

    a_c = (1 / pi) * ( K_IC / (Y * sigma_max) ) ** 2      (sigma_max > 0)

Residual strength
-----------------
Inverting the same equation for the allowable maximum tensile stress at a given
crack size:

    sigma_residual(a) = K_IC / ( Y * sqrt(pi * a) )

Maximum-stress policy
---------------------
Mode-I fracture requires a tensile opening stress. If ``sigma_max <= 0`` the
cycle never opens the crack, no tensile mode-I fracture boundary exists within
this simple model, and the critical crack size is reported as ``math.inf`` with
an explicit :class:`FractureBoundaryStatus`. ``abs(sigma_max)`` is NEVER taken:
compression must not manufacture a fictitious opening-fracture boundary.

Scope
-----
Pure LEFM. No plastic-zone correction is applied -- see
:mod:`crackgrowth` limitations. If small-scale yielding is violated (a large
plastic zone relative to crack length, ligament, or thickness), a K_IC-based
critical crack size becomes questionable and this screen is not valid. No
thickness or state-of-stress validity check, no plane-stress/plane-strain
transition, no net-section collapse check, and no finite-width or ligament
correction is performed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from ._validation import require_finite, require_positive
from .geometry import ThroughCrackGeometry
from .loading import StressCycle
from .paris import MPA_ROOT_M_IN_PA_ROOT_M
from .stress_intensity import stress_intensity

__all__ = [
    "FractureToughness",
    "FractureBoundaryStatus",
    "CriticalCrackResult",
    "FractureAssessment",
    "ILLUSTRATIVE_TOUGHNESS_DISCLAIMER",
    "ILLUSTRATIVE_ALUMINIUM_LIKE_TOUGHNESS",
    "mpa_sqrt_m_to_pa_sqrt_m",
    "pa_sqrt_m_to_mpa_sqrt_m",
    "k_max",
    "critical_crack_length",
    "residual_strength",
    "fracture_utilization",
    "fracture_margin",
    "assess_fracture",
]

#: Mandatory provenance banner for any toughness that is not a verified dataset.
ILLUSTRATIVE_TOUGHNESS_DISCLAIMER = (
    "ILLUSTRATIVE FRACTURE-TOUGHNESS INPUT - NOT DESIGN ALLOWABLE"
)


def mpa_sqrt_m_to_pa_sqrt_m(value_mpa_sqrt_m: float) -> float:
    """Convert a stress intensity from MPa*sqrt(m) to Pa*sqrt(m).

    ``1 MPa*sqrt(m) = 1e6 Pa*sqrt(m)``. Handbook toughness values are quoted in
    MPa*sqrt(m); this package stores them in Pa*sqrt(m). The conversion is named
    and explicit so the basis is never inferred.
    """
    return require_finite(value_mpa_sqrt_m, "value_mpa_sqrt_m") * MPA_ROOT_M_IN_PA_ROOT_M


def pa_sqrt_m_to_mpa_sqrt_m(value_pa_sqrt_m: float) -> float:
    """Inverse of :func:`mpa_sqrt_m_to_pa_sqrt_m`."""
    return require_finite(value_pa_sqrt_m, "value_pa_sqrt_m") / MPA_ROOT_M_IN_PA_ROOT_M


@dataclass(frozen=True)
class FractureToughness:
    """Immutable mode-I fracture-toughness record, stored on the SI basis.

    Parameters
    ----------
    name:
        Short identifier for the material record.
    k_ic:
        Mode-I fracture toughness ``K_IC`` in Pa*sqrt(m). Finite, > 0. Use
        :func:`mpa_sqrt_m_to_pa_sqrt_m` to convert a handbook value.
    source_note:
        Provenance statement. For any value that is not traceable to a verified
        dataset this must carry :data:`ILLUSTRATIVE_TOUGHNESS_DISCLAIMER`.
    notes:
        Optional additional commentary.
    """

    name: str
    k_ic: float
    source_note: str
    notes: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("name must be a non-empty string")
        object.__setattr__(self, "k_ic", require_positive(self.k_ic, "k_ic"))
        if not isinstance(self.source_note, str) or not self.source_note.strip():
            raise ValueError("source_note must be a non-empty provenance statement")
        if not isinstance(self.notes, str):
            raise TypeError("notes must be a string")

    @property
    def k_ic_in_mpa_sqrt_m(self) -> float:
        """This toughness re-expressed in MPa*sqrt(m), for reporting."""
        return pa_sqrt_m_to_mpa_sqrt_m(self.k_ic)


class FractureBoundaryStatus(Enum):
    """Whether a tensile mode-I fracture boundary exists for a given cycle."""

    #: ``sigma_max > 0``: a finite critical crack size exists.
    FINITE_CRITICAL_SIZE = "finite critical crack size"
    #: ``sigma_max <= 0``: the cycle never opens the crack in mode I.
    NO_TENSILE_BOUNDARY = "no tensile mode-I fracture boundary"


@dataclass(frozen=True)
class CriticalCrackResult:
    """The toughness-derived fracture boundary for one cycle and geometry."""

    critical_crack_length: float
    """Critical crack length [m]; ``math.inf`` when no tensile boundary exists."""

    sigma_max: float
    """Maximum cycle stress used for the fracture criterion [Pa]."""

    geometry_factor: float
    """Constant geometry factor ``Y`` [-]."""

    k_ic: float
    """Mode-I fracture toughness [Pa*sqrt(m)]."""

    status: FractureBoundaryStatus
    """Whether the critical size is finite or no tensile boundary exists."""

    @property
    def has_finite_boundary(self) -> bool:
        return self.status is FractureBoundaryStatus.FINITE_CRITICAL_SIZE


@dataclass(frozen=True)
class FractureAssessment:
    """A residual-strength / fracture screen at one crack length."""

    crack_length: float
    """Crack length ``a`` [m]."""

    applied_sigma_max: float
    """Maximum applied tensile stress [Pa]."""

    k_max: float
    """Signed maximum stress intensity ``Y*sigma_max*sqrt(pi*a)`` [Pa*sqrt(m)]."""

    k_ic: float
    """Mode-I fracture toughness [Pa*sqrt(m)]."""

    residual_strength: float
    """Allowable maximum tensile stress at this crack length [Pa]."""

    utilization: float
    """``K_max / K_IC`` [-]. Zero for a non-opening (compressive) cycle."""

    margin: float
    """``K_IC / K_max - 1`` [-]. ``math.inf`` for a non-opening cycle."""

    passes: bool
    """True when ``K_max <= K_IC``. The boundary itself passes."""

    status: FractureBoundaryStatus
    """Whether the cycle opens the crack in mode I at all."""


def k_max(
    crack_length: float,
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
) -> float:
    """Maximum mode-I stress intensity in the cycle [Pa*sqrt(m)].

    ``K_max = Y * sigma_max * sqrt(pi * a)``, evaluated through the verified
    Milestone 1 :func:`~crackgrowth.stress_intensity.stress_intensity` helper and
    therefore SIGNED: negative when ``sigma_max`` is compressive. This is the
    FRACTURE driving force, and is not ``delta_K``.
    """
    if not isinstance(cycle, StressCycle):
        raise TypeError(f"cycle must be a StressCycle, got {cycle!r}")
    return stress_intensity(cycle.sigma_max, crack_length, geometry)


def critical_crack_length(
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    toughness: FractureToughness,
) -> CriticalCrackResult:
    """Crack length at which ``K_max`` reaches ``K_IC``.

        a_c = (1 / pi) * ( K_IC / (Y * sigma_max) ) ** 2

    Returns a structured result. When ``sigma_max <= 0`` there is no tensile
    mode-I fracture boundary: the critical length is ``math.inf`` and the status
    says so. ``abs(sigma_max)`` is never taken.
    """
    _require_types(cycle, geometry, toughness)
    if cycle.sigma_max <= 0.0:
        return CriticalCrackResult(
            critical_crack_length=math.inf,
            sigma_max=cycle.sigma_max,
            geometry_factor=geometry.geometry_factor,
            k_ic=toughness.k_ic,
            status=FractureBoundaryStatus.NO_TENSILE_BOUNDARY,
        )
    a_c = (
        toughness.k_ic / (geometry.geometry_factor * cycle.sigma_max)
    ) ** 2 / math.pi
    return CriticalCrackResult(
        critical_crack_length=a_c,
        sigma_max=cycle.sigma_max,
        geometry_factor=geometry.geometry_factor,
        k_ic=toughness.k_ic,
        status=FractureBoundaryStatus.FINITE_CRITICAL_SIZE,
    )


def residual_strength(
    crack_length: float,
    geometry: ThroughCrackGeometry,
    toughness: FractureToughness,
) -> float:
    """Allowable maximum tensile stress at a given crack length [Pa].

        sigma_residual(a) = K_IC / ( Y * sqrt(pi * a) )

    This is the same fracture equation as :func:`critical_crack_length`, solved
    for stress instead of crack size, so evaluating it at ``a_c`` returns
    ``sigma_max`` exactly. It is an LEFM screening value, NOT a certified
    residual-strength allowable: no net-section collapse, ligament, or
    finite-width check is performed, and for large cracks or low toughness the
    net section may yield long before this stress is reached.
    """
    if not isinstance(geometry, ThroughCrackGeometry):
        raise TypeError(f"geometry must be a ThroughCrackGeometry, got {geometry!r}")
    if not isinstance(toughness, FractureToughness):
        raise TypeError(f"toughness must be a FractureToughness, got {toughness!r}")
    a = require_positive(crack_length, "crack_length")
    return toughness.k_ic / (geometry.geometry_factor * math.sqrt(math.pi * a))


def fracture_utilization(k_max_value: float, toughness: FractureToughness) -> float:
    """``K_max / K_IC`` [-].

    A non-opening cycle (``K_max <= 0``) has no mode-I fracture demand and
    returns ``0.0``: the criterion is non-governing, not merely satisfied.
    """
    if not isinstance(toughness, FractureToughness):
        raise TypeError(f"toughness must be a FractureToughness, got {toughness!r}")
    k = require_finite(k_max_value, "k_max_value")
    if k <= 0.0:
        return 0.0
    return k / toughness.k_ic


def fracture_margin(k_max_value: float, toughness: FractureToughness) -> float:
    """Preliminary fracture margin ``MS_K = K_IC / K_max - 1`` [-].

    Positive means ``K_max`` is below toughness; exactly zero at the boundary.
    A non-opening cycle (``K_max <= 0``) returns ``math.inf``, documenting that
    the mode-I fracture criterion is non-governing rather than infinitely safe.

    This is a screening indicator only. It is NOT a certification margin of
    safety: no load factor, no scatter allowance, no material variability, no
    thickness/state-of-stress validity check and no net-section check are
    included.
    """
    if not isinstance(toughness, FractureToughness):
        raise TypeError(f"toughness must be a FractureToughness, got {toughness!r}")
    k = require_finite(k_max_value, "k_max_value")
    if k <= 0.0:
        return math.inf
    return toughness.k_ic / k - 1.0


def assess_fracture(
    crack_length: float,
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    toughness: FractureToughness,
) -> FractureAssessment:
    """Screen one crack length against the mode-I fracture boundary.

    ``passes`` is True when ``K_max <= K_IC``; the boundary itself passes.
    """
    _require_types(cycle, geometry, toughness)
    a = require_positive(crack_length, "crack_length")
    k = k_max(a, cycle, geometry)
    opening = cycle.sigma_max > 0.0
    return FractureAssessment(
        crack_length=a,
        applied_sigma_max=cycle.sigma_max,
        k_max=k,
        k_ic=toughness.k_ic,
        residual_strength=residual_strength(a, geometry, toughness),
        utilization=fracture_utilization(k, toughness),
        margin=fracture_margin(k, toughness),
        passes=k <= toughness.k_ic,
        status=(
            FractureBoundaryStatus.FINITE_CRITICAL_SIZE
            if opening
            else FractureBoundaryStatus.NO_TENSILE_BOUNDARY
        ),
    )


def _require_types(
    cycle: StressCycle,
    geometry: ThroughCrackGeometry,
    toughness: FractureToughness,
) -> None:
    if not isinstance(cycle, StressCycle):
        raise TypeError(f"cycle must be a StressCycle, got {cycle!r}")
    if not isinstance(geometry, ThroughCrackGeometry):
        raise TypeError(f"geometry must be a ThroughCrackGeometry, got {geometry!r}")
    if not isinstance(toughness, FractureToughness):
        raise TypeError(f"toughness must be a FractureToughness, got {toughness!r}")


#: Canonical Milestone 2 toughness.
#:
#: ILLUSTRATIVE. Selected after auditing K_IC = 20, 25, 30, 35, 40, 50 and
#: 60 MPa*sqrt(m) against the canonical cycle: 30 MPa*sqrt(m) places the critical
#: crack size at 19.89 mm -- comfortably above the 1 mm initial flaw, within a
#: plausible inspectable crack-size range, and deliberately NOT equal to
#: Milestone 1's arbitrary 10 mm target. Not measured, not traceable to a
#: qualification dataset, and no specific alloy is claimed.
ILLUSTRATIVE_ALUMINIUM_LIKE_TOUGHNESS = FractureToughness(
    name="Illustrative aluminium-like plane-strain toughness",
    k_ic=mpa_sqrt_m_to_pa_sqrt_m(30.0),
    source_note=(
        ILLUSTRATIVE_TOUGHNESS_DISCLAIMER
        + ": K_IC = 30 MPa*sqrt(m), a round value chosen after auditing the "
        "resulting critical crack size and life for the canonical wing-skin case. "
        "Not measured, not traceable to a qualification dataset, and no specific "
        "alloy is claimed."
    ),
    notes=(
        "Audited at sigma_max = 120 MPa, Y = 1.0: a_c = 19.89 mm, "
        "delta_K(a_c) = 25.0 MPa*sqrt(m), life from 1 mm = 8.81e5 cycles. "
        "Treated as a single constant value: no thickness or state-of-stress "
        "validity check, and no plane-stress/plane-strain transition is modelled."
    ),
)
