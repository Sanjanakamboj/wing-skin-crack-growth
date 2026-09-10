"""Paris-law fatigue crack growth for an idealized aircraft wing-skin crack.

Milestone 1: linear-elastic fracture mechanics with a constant geometry factor,
constant-amplitude loading, and a Paris-law growth model, integrated from an
initial crack length to an IMPOSED target crack length.

Milestone 2: a mode-I fracture-toughness screen -- critical crack size from
``K_IC``, residual strength, fracture margin -- and crack-growth life integrated
to that physics-derived boundary instead of an imposed one. Growth is driven by
``delta_K``; fracture is driven by ``K_max``. The two are never interchanged.

All quantities are SI internally:

===============================  ==========================================
crack length ``a``               m
stress ``sigma``                 Pa (tension positive)
stress intensity ``K``           Pa*sqrt(m)
Paris coefficient ``C``          m/cycle / (Pa*sqrt(m))**m
crack-growth rate ``da/dN``      m/cycle
life ``N``                       cycles
===============================  ==========================================
"""

from __future__ import annotations

from .canonical import (
    CANONICAL_CYCLE,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GEOMETRY,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CANONICAL_TARGET_CRACK_LENGTH,
)
from .fracture import (
    ILLUSTRATIVE_ALUMINIUM_LIKE_TOUGHNESS,
    ILLUSTRATIVE_TOUGHNESS_DISCLAIMER,
    CriticalCrackResult,
    FractureAssessment,
    FractureBoundaryStatus,
    FractureToughness,
    assess_fracture,
    critical_crack_length,
    fracture_margin,
    fracture_utilization,
    k_max,
    mpa_sqrt_m_to_pa_sqrt_m,
    pa_sqrt_m_to_mpa_sqrt_m,
    residual_strength,
)
from .fracture_life import (
    FlawAdmissibility,
    InitialFlawBeyondCriticalError,
    LifeToFractureResult,
    assess_initial_flaw,
    cycles_to_critical_crack,
)
from .fracture_sensitivity import (
    FracturePoint,
    ResidualStrengthRow,
    geometry_factor_sensitivity,
    initial_crack_sensitivity_to_fracture,
    max_stress_sensitivity_at_fixed_range,
    residual_strength_table,
    stress_range_sensitivity_at_fixed_min,
    toughness_sensitivity,
)
from .geometry import INFINITE_PLATE_THROUGH_CRACK, ThroughCrackGeometry
from .integration import (
    DEFAULT_INTEGRATION_INTERVALS,
    CrackGrowthResult,
    analytical_cycles_to_crack_length,
    cycles_to_crack_length,
)
from .loading import StressCycle
from .paris import (
    ILLUSTRATIVE_ALUMINIUM_LIKE_PARIS,
    ILLUSTRATIVE_DISCLAIMER,
    MPA_ROOT_M_IN_PA_ROOT_M,
    CrackGrowthRatePoint,
    ParisLaw,
    crack_growth_rate,
    crack_growth_rate_point,
    paris_c_from_mpa_basis,
    paris_c_to_mpa_basis,
    paris_law_from_mpa_basis,
)
from .sensitivity import (
    SensitivityPoint,
    initial_crack_length_sensitivity,
    stress_range_sensitivity,
)
from .stress_intensity import delta_stress_intensity, stress_intensity

__version__ = "0.2.0"

__all__ = [
    "__version__",
    # geometry
    "ThroughCrackGeometry",
    "INFINITE_PLATE_THROUGH_CRACK",
    # loading
    "StressCycle",
    # stress intensity
    "stress_intensity",
    "delta_stress_intensity",
    # paris
    "ParisLaw",
    "CrackGrowthRatePoint",
    "crack_growth_rate",
    "crack_growth_rate_point",
    "paris_c_from_mpa_basis",
    "paris_c_to_mpa_basis",
    "paris_law_from_mpa_basis",
    "MPA_ROOT_M_IN_PA_ROOT_M",
    "ILLUSTRATIVE_DISCLAIMER",
    "ILLUSTRATIVE_ALUMINIUM_LIKE_PARIS",
    # integration
    "CrackGrowthResult",
    "cycles_to_crack_length",
    "analytical_cycles_to_crack_length",
    "DEFAULT_INTEGRATION_INTERVALS",
    # sensitivity
    "SensitivityPoint",
    "stress_range_sensitivity",
    "initial_crack_length_sensitivity",
    # canonical case
    "CANONICAL_GEOMETRY",
    "CANONICAL_CYCLE",
    "CANONICAL_PARIS_LAW",
    "CANONICAL_INITIAL_CRACK_LENGTH",
    "CANONICAL_TARGET_CRACK_LENGTH",
    "CANONICAL_FRACTURE_TOUGHNESS",
    # fracture toughness and residual strength (Milestone 2)
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
    # life to the fracture boundary (Milestone 2)
    "FlawAdmissibility",
    "InitialFlawBeyondCriticalError",
    "LifeToFractureResult",
    "assess_initial_flaw",
    "cycles_to_critical_crack",
    # fracture sensitivity (Milestone 2)
    "FracturePoint",
    "ResidualStrengthRow",
    "toughness_sensitivity",
    "max_stress_sensitivity_at_fixed_range",
    "stress_range_sensitivity_at_fixed_min",
    "geometry_factor_sensitivity",
    "initial_crack_sensitivity_to_fracture",
    "residual_strength_table",
]
