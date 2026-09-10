"""Paris-law fatigue crack growth for an idealized aircraft wing-skin crack.

Milestone 1: linear-elastic fracture mechanics with a constant geometry factor,
constant-amplitude loading, and a Paris-law growth model, integrated from an
initial crack length to an IMPOSED target crack length.

Milestone 2: a mode-I fracture-toughness screen -- critical crack size from
``K_IC``, residual strength, fracture margin -- and crack-growth life integrated
to that physics-derived boundary instead of an imposed one. Growth is driven by
``delta_K``; fracture is driven by ``K_max``. The two are never interchanged.

Milestone 3: a finite-width centre-cracked panel with a crack-size-dependent
geometry factor ``Y(a) = 1/sqrt(cos(pi a / W))``, a numerically solved fracture
boundary (no closed form survives ``Y(a)``), and a geometry-aware log-grid
integrator. Throughout, ``a`` is the HALF crack length and the total crack
length is ``2a``.

Milestone 4: a constant-amplitude crack-growth threshold ``delta_K_th`` applied
as a HARD CUTOFF, distinguishing a crack that is actively growing from one that
is arrested below threshold or already at the fracture boundary. The
no-threshold Paris path remains available and unchanged.

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
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GEOMETRY,
    CANONICAL_GROWTH_THRESHOLD,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CANONICAL_PLATE_WIDTH,
    CANONICAL_TARGET_CRACK_LENGTH,
)
from .threshold import (
    ILLUSTRATIVE_ALUMINIUM_LIKE_THRESHOLD,
    ILLUSTRATIVE_THRESHOLD_DISCLAIMER,
    CrackGrowthThreshold,
    FiniteWidthThresholdResult,
    ThresholdBoundaryStatus,
    ThresholdedGrowthPoint,
    finite_width_threshold_crack_length,
    threshold_crack_length_constant_y,
    thresholded_crack_growth_rate,
)
from .threshold_life import (
    BoundaryOrdering,
    GrowthState,
    ThresholdedLifeResult,
    classify_growth_state,
    compare_threshold_and_fracture_boundaries,
    cycles_to_fracture_with_threshold,
)
from .threshold_sensitivity import (
    ThresholdPoint,
    initial_crack_sensitivity_with_threshold,
    max_stress_sensitivity_with_threshold,
    stress_range_sensitivity_with_threshold,
    threshold_sensitivity,
    toughness_sensitivity_with_threshold,
    width_sensitivity_with_threshold,
)
from .finite_width import (
    CrackGeometry,
    FiniteWidthCenterCrack,
    crack_growth_rate_for_geometry,
    delta_k_for_geometry,
    geometry_factor_at,
    ligament_fraction,
    remaining_ligament,
    stress_intensity_for_geometry,
    total_crack_length,
)
from .finite_width_fracture import (
    DEFAULT_ROOT_TOLERANCE,
    UPPER_BOUND_MARGIN,
    CriticalCrackStatus,
    FiniteWidthCriticalCrackResult,
    FiniteWidthFractureAssessment,
    assess_fracture_for_geometry,
    finite_width_critical_crack_length,
    k_max_for_geometry,
    residual_strength_for_geometry,
)
from .finite_width_life import (
    FiniteWidthLifeResult,
    cycles_to_finite_width_fracture,
    frozen_geometry_factor_life,
)
from .finite_width_sensitivity import (
    FiniteWidthPoint,
    GeometryAmplificationRow,
    geometry_amplification_table,
    initial_crack_sensitivity_finite_width,
    toughness_sensitivity_finite_width,
    width_sensitivity,
)
from .finite_width_sensitivity import (
    max_stress_sensitivity_at_fixed_range as finite_width_max_stress_sensitivity,
)
from .finite_width_sensitivity import (
    stress_range_sensitivity_at_fixed_min as finite_width_stress_range_sensitivity,
)
from .integration_log import (
    DEFAULT_LOG_INTERVALS,
    LogGridGrowthResult,
    integrate_crack_growth_log_grid,
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

__version__ = "0.4.0"

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
    # finite-width geometry (Milestone 3)
    "CrackGeometry",
    "FiniteWidthCenterCrack",
    "CANONICAL_PLATE_WIDTH",
    "CANONICAL_FINITE_WIDTH_GEOMETRY",
    "total_crack_length",
    "remaining_ligament",
    "ligament_fraction",
    "geometry_factor_at",
    "stress_intensity_for_geometry",
    "delta_k_for_geometry",
    "crack_growth_rate_for_geometry",
    # finite-width fracture boundary (Milestone 3)
    "CriticalCrackStatus",
    "FiniteWidthCriticalCrackResult",
    "FiniteWidthFractureAssessment",
    "k_max_for_geometry",
    "residual_strength_for_geometry",
    "finite_width_critical_crack_length",
    "assess_fracture_for_geometry",
    "DEFAULT_ROOT_TOLERANCE",
    "UPPER_BOUND_MARGIN",
    # geometry-aware log-grid integration (Milestone 3)
    "LogGridGrowthResult",
    "integrate_crack_growth_log_grid",
    "DEFAULT_LOG_INTERVALS",
    # finite-width life (Milestone 3)
    "FiniteWidthLifeResult",
    "cycles_to_finite_width_fracture",
    "frozen_geometry_factor_life",
    # finite-width sensitivity (Milestone 3)
    "FiniteWidthPoint",
    "GeometryAmplificationRow",
    "width_sensitivity",
    "initial_crack_sensitivity_finite_width",
    "finite_width_max_stress_sensitivity",
    "finite_width_stress_range_sensitivity",
    "toughness_sensitivity_finite_width",
    "geometry_amplification_table",
    # crack-growth threshold (Milestone 4)
    "CrackGrowthThreshold",
    "CANONICAL_GROWTH_THRESHOLD",
    "ILLUSTRATIVE_THRESHOLD_DISCLAIMER",
    "ILLUSTRATIVE_ALUMINIUM_LIKE_THRESHOLD",
    "ThresholdedGrowthPoint",
    "ThresholdBoundaryStatus",
    "FiniteWidthThresholdResult",
    "thresholded_crack_growth_rate",
    "threshold_crack_length_constant_y",
    "finite_width_threshold_crack_length",
    # threshold-aware life (Milestone 4)
    "GrowthState",
    "BoundaryOrdering",
    "ThresholdedLifeResult",
    "classify_growth_state",
    "compare_threshold_and_fracture_boundaries",
    "cycles_to_fracture_with_threshold",
    # threshold sensitivity (Milestone 4)
    "ThresholdPoint",
    "threshold_sensitivity",
    "initial_crack_sensitivity_with_threshold",
    "stress_range_sensitivity_with_threshold",
    "max_stress_sensitivity_with_threshold",
    "width_sensitivity_with_threshold",
    "toughness_sensitivity_with_threshold",
]
