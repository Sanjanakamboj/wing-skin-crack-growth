"""Milestone 4: crack-growth threshold screening for a finite-width wing skin.

Run with::

    python examples/crack_growth_threshold.py

All internal computation is SI. Engineering units (mm, MPa, MPa*sqrt(m)) appear
only in the printed output, where the conversion is explicit.

``a`` is the HALF crack length throughout; the total crack length is ``2a``.
This example adds a threshold screen on top of the Milestone 3 finite-width
model, which is left unchanged and still runs on its own.
"""

from __future__ import annotations

import math

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GROWTH_THRESHOLD,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CANONICAL_PLATE_WIDTH,
    CrackGrowthThreshold,
    cycles_to_finite_width_fracture,
    cycles_to_fracture_with_threshold,
    delta_k_for_geometry,
    initial_crack_sensitivity_with_threshold,
    k_max_for_geometry,
    max_stress_sensitivity_with_threshold,
    mpa_sqrt_m_to_pa_sqrt_m,
    stress_range_sensitivity_with_threshold,
    threshold_sensitivity,
    toughness_sensitivity_with_threshold,
    width_sensitivity_with_threshold,
)

MM_PER_M = 1.0e3
MPA_PER_PA = 1.0e-6

THRESHOLD_SWEEP_MPA = (None, 2.0, 3.0, 4.0, 5.0, 5.5, 5.6, 5.7, 6.0, 8.0, 10.0)
INITIAL_CRACKS_M = (0.1e-3, 0.25e-3, 0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3, 8.0e-3)
SIGMA_MAX_FIXED_MIN_MPA = (40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0)
SIGMA_MAX_FIXED_RANGE_MPA = (100.0, 120.0, 140.0, 160.0, 180.0)
WIDTHS_M = (40.0e-3, 50.0e-3, 75.0e-3, 100.0e-3, 150.0e-3, 200.0e-3, 500.0e-3)
K_IC_SWEEP_MPA = (20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0)


def _rule(title: str) -> None:
    print()
    print(title)
    print("-" * len(title))


def _cycles(value: float) -> str:
    return "arrested (inf)" if math.isinf(value) else f"{value:,.1f}"


def _state(value) -> str:
    return value.name


def _threshold(mpa: float) -> CrackGrowthThreshold:
    return CrackGrowthThreshold(
        name=f"swept {mpa} MPa*sqrt(m)",
        delta_k_threshold=mpa_sqrt_m_to_pa_sqrt_m(mpa),
        source_note=CANONICAL_GROWTH_THRESHOLD.source_note,
    )


def main() -> None:
    panel = CANONICAL_FINITE_WIDTH_GEOMETRY
    cycle = CANONICAL_CYCLE
    law = CANONICAL_PARIS_LAW
    toughness = CANONICAL_FRACTURE_TOUGHNESS
    threshold = CANONICAL_GROWTH_THRESHOLD
    a0 = CANONICAL_INITIAL_CRACK_LENGTH
    width = CANONICAL_PLATE_WIDTH

    print("=" * 78)
    print("CRACK-GROWTH THRESHOLD SCREENING (MILESTONE 4)")
    print("=" * 78)

    # ------------------------------------------------------------------
    _rule("STUDY BASIS")
    print("  Crack convention      : a = HALF crack length; total crack = 2a")
    print("  Geometry              : centre through crack, finite-width panel")
    print("  Geometry factor       : Y(a) = 1 / sqrt(cos(pi * a / W))")
    print(f"  Plate width W         : {width * MM_PER_M:10.4f}  mm")
    print(f"  a0 (assumed flaw)     : {a0 * MM_PER_M:10.4f}  mm")
    print(f"  sigma_max             : {cycle.sigma_max * MPA_PER_PA:10.4f}  MPa")
    print(f"  sigma_min             : {cycle.sigma_min * MPA_PER_PA:10.4f}  MPa")
    print(f"  delta_sigma           : {cycle.stress_range * MPA_PER_PA:10.4f}  MPa")
    print(f"  stress ratio R        : {cycle.stress_ratio:10.4f}  [-]")
    print(f"  Paris C  (SI basis)   : {law.C:10.4e}  m/cycle/(Pa*sqrt(m))**m")
    print(f"  Paris m               : {law.m:10.4f}  [-]")
    print("  K_IC                  : "
          f"{toughness.k_ic_in_mpa_sqrt_m:10.4f}  MPa*sqrt(m)")
    print("  delta_K_th            : "
          f"{threshold.delta_k_threshold_in_mpa_sqrt_m:10.4f}  MPa*sqrt(m)")
    print("  Threshold policy      : HARD CUTOFF -- da/dN = 0 for")
    print("                          delta_K <= delta_K_th, and the UNMODIFIED")
    print("                          Paris rate above it. Equality arrests.")
    print()
    print("  PARIS PROVENANCE")
    print(f"    {law.source_note}")
    print()
    print("  TOUGHNESS PROVENANCE")
    print(f"    {toughness.source_note}")
    print()
    print("  THRESHOLD PROVENANCE")
    print(f"    {threshold.source_note}")

    # ------------------------------------------------------------------
    result = cycles_to_fracture_with_threshold(
        a0, cycle, panel, law, toughness, threshold
    )
    milestone3 = cycles_to_finite_width_fracture(
        a0, cycle, panel, law, toughness
    )

    _rule("INITIAL CRACK STATE")
    print(f"  Y(a0)                 : {panel.geometry_factor_at(a0):10.6f}  [-]")
    print("  K_max(a0)             : "
          f"{k_max_for_geometry(a0, cycle, panel) * MPA_PER_PA:10.4f}  MPa*sqrt(m)")
    print("  delta_K(a0)           : "
          f"{result.initial_delta_k * MPA_PER_PA:10.4f}  MPa*sqrt(m)")
    print(f"  delta_K(a0)/delta_K_th: {result.threshold_ratio:10.6f}  [-]")
    print("  Threshold excess      : "
          f"{result.threshold_excess_ratio:10.6f}  [-]  (POSITIVE = growing,")
    print("                                            NOT a safety margin)")
    print(f"  STATE                 : {result.state.name}")
    print(f"                          {result.state.value}")
    print("  Fracture utilization  : "
          f"{result.initial_delta_k / result.delta_k_at_critical:10.6f}  "
          "[-]  (delta_K basis)")

    # ------------------------------------------------------------------
    a_th = result.threshold_crack_length
    a_c = result.critical_crack_length
    boundary = result.threshold_boundary

    _rule("BOUNDARIES")
    print(f"  Threshold crack a_th  : {a_th * MM_PER_M:10.4f}  mm")
    print(f"  Fracture crack a_c    : {a_c * MM_PER_M:10.4f}  mm")
    print(f"  a_th / a_c            : {a_th / a_c:10.6f}  [-]")
    print(f"  Ordering              : {result.boundary_ordering.name}")
    print(f"                          {result.boundary_ordering.value}")
    print("  Ligament frac at a_th : "
          f"{boundary.ligament_fraction:10.6f}  [-]  (diagnostic only)")
    print("  Ligament frac at a_c  : "
          f"{result.critical.ligament_fraction:10.6f}  [-]  (diagnostic only)")
    print("  delta_K(a_th)         : "
          f"{boundary.delta_k_at_threshold * MPA_PER_PA:10.4f}  MPa*sqrt(m)"
          "  (= delta_K_th)")
    print("  delta_K(a_c)          : "
          f"{result.delta_k_at_critical * MPA_PER_PA:10.4f}  MPa*sqrt(m)")
    print("  delta_K(a_c)/K_IC     : "
          f"{result.delta_k_at_critical / toughness.k_ic:10.6f}  [-]  "
          f"(= delta_sigma/sigma_max = {cycle.stress_range / cycle.sigma_max:.6f})")
    print("  delta_K_th / K_IC     : "
          f"{result.threshold_to_toughness_ratio:10.6f}  [-]  (DIAGNOSTIC only,")
    print("                                            never a criterion)")
    print("  A growth interval exists only while delta_K_th < delta_K(a_c) for")
    print("  this monotonic geometry; here 4.00 < 25.00 MPa*sqrt(m).")

    # ------------------------------------------------------------------
    _rule("LIFE")
    print(f"  Threshold-aware life  : {_cycles(result.predicted_cycles):>18}  cycles")
    print(f"  M3 no-threshold life  : "
          f"{_cycles(milestone3.predicted_cycles):>18}  cycles")
    difference = result.predicted_cycles - milestone3.predicted_cycles
    print(f"  Difference            : {difference:18.6f}  cycles")
    print(f"  Bit-for-bit identical : "
          f"{str(result.predicted_cycles == milestone3.predicted_cycles):>18}")
    print(f"  Classification        : {result.state.name}")
    print()
    print("  KEY RESULT: because the hard cutoff does NOT modify the Paris rate")
    print("  above threshold, an already-growing crack has EXACTLY the M3 life.")
    print("  The threshold changes the verdict only when the initial crack is")
    print("  at or below it -- and then the life is infinite, not merely longer.")
    print("  There is no gradual life extension; the model is discontinuous.")

    # ------------------------------------------------------------------
    _rule("THRESHOLD SENSITIVITY (a0, W and cycle fixed)")
    print(f"  {'dK_th':>9}  {'ratio@a0':>9}  {'a_th [mm]':>10}  "
          f"{'state':>24}  {'life':>15}")
    for point in threshold_sensitivity(
        tuple(None if v is None else _threshold(v) for v in THRESHOLD_SWEEP_MPA),
        a0, cycle, panel, law, toughness,
    ):
        label = "disabled" if point.delta_k_threshold is None else (
            f"{point.delta_k_threshold * MPA_PER_PA:.1f}"
        )
        ratio = "-" if point.threshold_ratio is None else f"{point.threshold_ratio:.6f}"
        size = "-" if point.threshold_crack_length is None else (
            f"{point.threshold_crack_length * MM_PER_M:.4f}"
        )
        print(f"  {label:>9}  {ratio:>9}  {size:>10}  {_state(point.state):>24}  "
              f"{_cycles(point.predicted_cycles):>15}")
    print("  The life takes exactly TWO values across the whole sweep. The")
    print("  crossing sits between 5.6 and 5.7, i.e. at delta_K(a0) = 5.6064.")

    # ------------------------------------------------------------------
    _rule("INITIAL-CRACK SENSITIVITY (delta_K_th = 4 MPa*sqrt(m))")
    print(f"  {'a0 [mm]':>8}  {'dK(a0)':>9}  {'ratio':>8}  {'state':>24}  "
          f"{'life':>15}")
    for point in initial_crack_sensitivity_with_threshold(
        INITIAL_CRACKS_M, cycle, panel, law, toughness, threshold
    ):
        print(f"  {point.parameter_value * MM_PER_M:8.2f}  "
              f"{point.initial_delta_k * MPA_PER_PA:9.4f}  "
              f"{point.threshold_ratio:8.4f}  {_state(point.state):>24}  "
              f"{_cycles(point.predicted_cycles):>15}")
    print("  A real transition: flaws at or below 0.5 mm are arrested, 1 mm and")
    print("  above grow. This was reported, not engineered.")

    # ------------------------------------------------------------------
    _rule("STRESS-RANGE SENSITIVITY (sigma_min = 20 MPa fixed)")
    print(f"  {'s_max':>7}  {'d_sigma':>8}  {'dK(a0)':>9}  {'ratio':>8}  "
          f"{'a_c [mm]':>9}  {'state':>24}  {'life':>15}")
    for point in stress_range_sensitivity_with_threshold(
        tuple(v * 1.0e6 for v in SIGMA_MAX_FIXED_MIN_MPA),
        20.0e6, a0, panel, law, toughness, threshold,
    ):
        print(f"  {point.sigma_max * MPA_PER_PA:7.1f}  "
              f"{point.stress_range * MPA_PER_PA:8.1f}  "
              f"{point.initial_delta_k * MPA_PER_PA:9.4f}  "
              f"{point.threshold_ratio:8.4f}  "
              f"{point.critical_crack_length * MM_PER_M:9.4f}  "
              f"{_state(point.state):>24}  {_cycles(point.predicted_cycles):>15}")
    print("  A second real transition, between sigma_max = 80 and 100 MPa. Both")
    print("  the driving range and the fracture boundary move in this sweep.")

    # ------------------------------------------------------------------
    _rule("FIXED delta_sigma = 100 MPa (threshold state frozen)")
    print(f"  {'s_max':>7}  {'dK(a0)':>10}  {'ratio':>9}  {'a_th [mm]':>10}  "
          f"{'a_c [mm]':>9}  {'life':>13}")
    for point in max_stress_sensitivity_with_threshold(
        tuple(v * 1.0e6 for v in SIGMA_MAX_FIXED_RANGE_MPA),
        100.0e6, a0, panel, law, toughness, threshold,
    ):
        print(f"  {point.sigma_max * MPA_PER_PA:7.1f}  "
              f"{point.initial_delta_k * MPA_PER_PA:10.6f}  "
              f"{point.threshold_ratio:9.6f}  "
              f"{point.threshold_crack_length * MM_PER_M:10.4f}  "
              f"{point.critical_crack_length * MM_PER_M:9.4f}  "
              f"{point.predicted_cycles:13.1f}")
    print("  delta_K(a0), the ratio and a_th are IDENTICAL down the column: the")
    print("  threshold sees only delta_sigma. The life changes purely because")
    print("  K_max moves the fracture endpoint.")

    # ------------------------------------------------------------------
    _rule("WIDTH SENSITIVITY")
    print(f"  {'W [mm]':>8}  {'Y(a0)':>9}  {'dK(a0)':>10}  {'ratio':>9}  "
          f"{'a_th [mm]':>10}  {'a_c [mm]':>9}  {'life':>13}")
    for point in width_sensitivity_with_threshold(
        WIDTHS_M, a0, cycle, law, toughness, threshold
    ):
        print(f"  {point.plate_width * MM_PER_M:8.1f}  "
              f"{point.initial_geometry_factor:9.6f}  "
              f"{point.initial_delta_k * MPA_PER_PA:10.6f}  "
              f"{point.threshold_ratio:9.6f}  "
              f"{point.threshold_crack_length * MM_PER_M:10.4f}  "
              f"{point.critical_crack_length * MM_PER_M:9.4f}  "
              f"{point.predicted_cycles:13.1f}")
    print("  Width raises Y(a0) and so delta_K(a0), but at a0/W <= 0.025 the")
    print("  ratio moves only ~0.15 %, so no width here crosses the threshold.")
    print("  The life still moves strongly, through a_c.")

    # ------------------------------------------------------------------
    _rule("TOUGHNESS SENSITIVITY (threshold state must not move)")
    print(f"  {'K_IC':>7}  {'dK(a0)':>10}  {'ratio':>9}  {'a_th [mm]':>10}  "
          f"{'a_c [mm]':>9}  {'life':>13}")
    for point in toughness_sensitivity_with_threshold(
        tuple(mpa_sqrt_m_to_pa_sqrt_m(v) for v in K_IC_SWEEP_MPA),
        a0, cycle, panel, law, toughness, threshold,
    ):
        print(f"  {point.k_ic * MPA_PER_PA:7.1f}  "
              f"{point.initial_delta_k * MPA_PER_PA:10.6f}  "
              f"{point.threshold_ratio:9.6f}  "
              f"{point.threshold_crack_length * MM_PER_M:10.4f}  "
              f"{point.critical_crack_length * MM_PER_M:9.4f}  "
              f"{point.predicted_cycles:13.1f}")
    print("  K_IC does not enter delta_K, so the threshold verdict and a_th are")
    print("  unchanged. Only the fracture endpoint, and hence the life, moves.")

    # ------------------------------------------------------------------
    _rule("INTERPRETATION")
    print("  * A threshold can arrest a small crack ENTIRELY in this")
    print("    constant-amplitude model. An arrested crack has da/dN = 0, so it")
    print("    can never grow up to a_th on its own: the reported life is")
    print("    infinite, and the crack is never advanced to a_th first.")
    print("  * Once the initial crack is above threshold, this hard-cutoff")
    print("    model gives EXACTLY the Milestone 3 Paris life. The threshold")
    print("    buys no extra cycles -- it only changes the verdict.")
    print("  * The model is therefore discontinuous: infinite life or the full")
    print("    no-threshold life, with nothing in between. That is a property")
    print("    of the chosen simplified cutoff, not of the material.")
    print("  * No crack closure is modelled. delta_K is still the algebraic")
    print("    range, so a compressive sigma_min inflates it and can make a")
    print("    crack look active that a closure-aware model would arrest.")
    print("  * The threshold is a single constant: no R-ratio dependence, no")
    print("    environmental, temperature or load-history dependence, and no")
    print("    near-threshold growth law.")
    print("  * delta_K_th is ILLUSTRATIVE, not a design allowable, and no")
    print("    specific alloy is claimed.")
    print("  * 'Infinite life' means ONLY that this model predicts zero")
    print("    propagation under this exact constant-amplitude cycle. It is not")
    print("    total structural durability, it says nothing about initiation,")
    print("    corrosion or fretting, and a threshold-arrested result is NOT a")
    print("    safe-life certification result.")
    print()


if __name__ == "__main__":
    main()
