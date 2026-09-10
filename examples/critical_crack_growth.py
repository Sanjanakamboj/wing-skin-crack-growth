"""Milestone 2: crack growth to a fracture-toughness-derived critical size.

Run with::

    python examples/critical_crack_growth.py

All internal computation is SI. Engineering units (mm, MPa, MPa*sqrt(m)) appear
only in the printed output, where the conversion is explicit.

This example supersedes the imposed 10 mm endpoint of
``examples/wing_skin_crack_growth.py`` with a boundary derived from ``K_IC``.
That Milestone 1 example is left unchanged and still runs.
"""

from __future__ import annotations

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GEOMETRY,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CANONICAL_TARGET_CRACK_LENGTH,
    analytical_cycles_to_crack_length,
    critical_crack_length,
    cycles_to_critical_crack,
    delta_stress_intensity,
    geometry_factor_sensitivity,
    initial_crack_sensitivity_to_fracture,
    k_max,
    max_stress_sensitivity_at_fixed_range,
    mpa_sqrt_m_to_pa_sqrt_m,
    residual_strength,
    residual_strength_table,
    stress_range_sensitivity_at_fixed_min,
    toughness_sensitivity,
)

MM_PER_M = 1.0e3
MPA_PER_PA = 1.0e-6

K_IC_SWEEP_MPA = (20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0)
SIGMA_MAX_SWEEP_MPA = (100.0, 120.0, 140.0, 160.0, 180.0)
STRESS_RANGE_SWEEP_MPA = (50.0, 75.0, 100.0, 125.0, 150.0)
GEOMETRY_FACTORS = (0.8, 0.9, 1.0, 1.1, 1.2)
INITIAL_CRACK_LENGTHS_M = (0.25e-3, 0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3)
TABLE_CRACK_LENGTHS_M = (0.5e-3, 1.0e-3, 2.0e-3, 5.0e-3, 10.0e-3, 20.0e-3)

# Wide sweeps reach a_c/a0 ratios near 80, where the fixed uniform Simpson grid
# needs more intervals to stay at full precision. See the README span note.
WIDE_SPAN_INTERVALS = 8000


def _rule(title: str) -> None:
    print()
    print(title)
    print("-" * len(title))


def main() -> None:
    geometry = CANONICAL_GEOMETRY
    cycle = CANONICAL_CYCLE
    law = CANONICAL_PARIS_LAW
    toughness = CANONICAL_FRACTURE_TOUGHNESS
    a0 = CANONICAL_INITIAL_CRACK_LENGTH

    print("=" * 78)
    print("WING-SKIN CRACK GROWTH TO A CRITICAL CRACK SIZE (MILESTONE 2)")
    print("=" * 78)

    # ------------------------------------------------------------------
    _rule("STUDY BASIS")
    print(f"  Geometry convention   : {geometry.description}")
    print("  K formulation         : K = Y * sigma * sqrt(pi * a)")
    print(f"  Geometry factor Y     : {geometry.geometry_factor:10.4f}  [-]")
    print(f"  sigma_max             : {cycle.sigma_max * MPA_PER_PA:10.4f}  MPa")
    print(f"  sigma_min             : {cycle.sigma_min * MPA_PER_PA:10.4f}  MPa")
    print(f"  delta_sigma           : {cycle.stress_range * MPA_PER_PA:10.4f}  MPa")
    print(f"  stress ratio R        : {cycle.stress_ratio:10.4f}  [-]")
    print(f"  Paris C  (SI basis)   : {law.C:10.4e}  m/cycle/(Pa*sqrt(m))**m")
    print(f"  Paris m               : {law.m:10.4f}  [-]")
    print("  K_IC                  : "
          f"{toughness.k_ic_in_mpa_sqrt_m:10.4f}  MPa*sqrt(m)")
    print(f"  a0 (assumed flaw)     : {a0 * MM_PER_M:10.4f}  mm")
    print()
    print("  PARIS PROVENANCE")
    print(f"    {law.source_note}")
    print()
    print("  TOUGHNESS PROVENANCE")
    print(f"    {toughness.source_note}")

    # ------------------------------------------------------------------
    result = cycles_to_critical_crack(a0, cycle, geometry, law, toughness)
    a_c = result.critical_crack_length

    _rule("FRACTURE BOUNDARY")
    print("  Criterion             : K_max = K_IC   (fracture uses K_max,")
    print("                          NOT delta_K -- growth uses delta_K)")
    print(f"  Critical crack a_c    : {a_c * MM_PER_M:10.4f}  mm")
    print(f"  Admissibility         : {result.admissibility.value}")
    print("  K_max(a0)             : "
          f"{k_max(a0, cycle, geometry) * MPA_PER_PA:10.4f}  MPa*sqrt(m)")
    print("  K_max(a_c)            : "
          f"{k_max(a_c, cycle, geometry) * MPA_PER_PA:10.4f}  MPa*sqrt(m)")
    print("  delta_K(a_c)          : "
          f"{result.delta_k_at_critical * MPA_PER_PA:10.4f}  MPa*sqrt(m)")
    print("  delta_K(a_c)/K_IC     : "
          f"{result.delta_k_at_critical / toughness.k_ic:10.4f}  "
          f"[-]  (= delta_sigma/sigma_max = "
          f"{cycle.stress_range / cycle.sigma_max:.4f})")
    print("  Residual strength(a0) : "
          f"{residual_strength(a0, geometry, toughness) * MPA_PER_PA:10.4f}  MPa")
    print("  Residual strength(a_c): "
          f"{residual_strength(a_c, geometry, toughness) * MPA_PER_PA:10.4f}  MPa"
          "   (= sigma_max, by construction)")
    print("  Initial utilization   : "
          f"{result.initial_assessment.utilization:10.4f}  [-]")
    print("  Initial margin MS_K   : "
          f"{result.initial_assessment.margin:10.4f}  [-]  (screening only)")

    # ------------------------------------------------------------------
    _rule("LIFE TO FRACTURE")
    growth = result.growth_result
    difference = abs(result.predicted_cycles - result.analytical_cycles)
    print(f"  Integration method    : composite Simpson, "
          f"{growth.integration_intervals} intervals (fixed grid)")
    print(f"  Numerical cycles      : {result.predicted_cycles:16.4f}  cycles")
    print(f"  Analytical cycles     : {result.analytical_cycles:16.4f}  cycles")
    print(f"  Relative difference   : "
          f"{difference / result.analytical_cycles:16.4e}  [-]")
    print(f"  da/dN(a0)             : {growth.initial_growth_rate:16.4e}  m/cycle")
    print(f"  da/dN(a_c)            : {growth.final_growth_rate:16.4e}  m/cycle")
    print("  Rate ratio            : "
          f"{growth.final_growth_rate / growth.initial_growth_rate:16.2f}  [-]")

    # ------------------------------------------------------------------
    _rule("COMPARISON TO MILESTONE 1")
    m1_life = analytical_cycles_to_crack_length(
        a0, CANONICAL_TARGET_CRACK_LENGTH, cycle, geometry, law
    )
    print(f"  M1 imposed target     : "
          f"{CANONICAL_TARGET_CRACK_LENGTH * MM_PER_M:10.4f}  mm  (arbitrary)")
    print(f"  M1 life to target     : {m1_life:16.1f}  cycles")
    print(f"  M2 critical size a_c  : {a_c * MM_PER_M:10.4f}  mm  "
          "(toughness-derived)")
    print(f"  M2 life to a_c        : {result.analytical_cycles:16.1f}  cycles")
    print(f"  Additional life       : {result.analytical_cycles - m1_life:16.1f}"
          "  cycles")
    print("  The M1 target was conservative here only by accident: it is an")
    print("  imposed number, and nothing guaranteed it fell below a_c.")

    # ------------------------------------------------------------------
    _rule("SENSITIVITY -- FRACTURE TOUGHNESS (all else fixed)")
    print(f"  {'K_IC':>10}  {'a_c [mm]':>10}  {'cycles':>14}  {'MS_K(a0)':>10}")
    for point in toughness_sensitivity(
        tuple(mpa_sqrt_m_to_pa_sqrt_m(v) for v in K_IC_SWEEP_MPA),
        a0, cycle, geometry, law, toughness, intervals=WIDE_SPAN_INTERVALS,
    ):
        print(f"  {point.k_ic * MPA_PER_PA:8.1f}    "
              f"{point.critical_crack_length * MM_PER_M:10.3f}  "
              f"{point.numerical_cycles:14.1f}  {point.initial_margin:10.3f}")
    print("  Higher toughness -> larger a_c -> longer life, but with strongly")
    print("  diminishing returns: doubling K_IC from 25 to 50 MPa*sqrt(m) adds")
    print("  only ~18 % life, because most cycles are spent while a is small.")

    # ------------------------------------------------------------------
    _rule("SENSITIVITY -- sigma_max AT FIXED delta_sigma = 100 MPa")
    print(f"  {'s_max':>7}  {'s_min':>7}  {'R':>7}  {'da/dN(a0)':>11}  "
          f"{'a_c [mm]':>9}  {'cycles':>13}")
    for point in max_stress_sensitivity_at_fixed_range(
        tuple(v * 1.0e6 for v in SIGMA_MAX_SWEEP_MPA),
        100.0e6, a0, geometry, law, toughness,
    ):
        s_min = point.sigma_max - point.stress_range
        print(f"  {point.sigma_max * MPA_PER_PA:7.1f}  {s_min * MPA_PER_PA:7.1f}  "
              f"{s_min / point.sigma_max:7.4f}  "
              f"{point.initial_growth_rate:11.4e}  "
              f"{point.critical_crack_length * MM_PER_M:9.3f}  "
              f"{point.numerical_cycles:13.1f}")
    print("  da/dN(a0) is IDENTICAL down the column: delta_K depends only on")
    print("  delta_sigma. The life falls purely because K_max rises and pulls")
    print("  the fracture boundary inward. This is the K_max / delta_K split.")

    # ------------------------------------------------------------------
    _rule("SENSITIVITY -- delta_sigma AT FIXED sigma_min = 20 MPa")
    print(f"  {'d_sigma':>8}  {'s_max':>7}  {'dK(a0)':>9}  {'da/dN(a0)':>11}  "
          f"{'a_c [mm]':>9}  {'cycles':>13}")
    for point in stress_range_sensitivity_at_fixed_min(
        tuple(v * 1.0e6 for v in STRESS_RANGE_SWEEP_MPA),
        20.0e6, a0, geometry, law, toughness,
    ):
        print(f"  {point.stress_range * MPA_PER_PA:8.1f}  "
              f"{point.sigma_max * MPA_PER_PA:7.1f}  "
              f"{point.initial_delta_k * MPA_PER_PA:9.4f}  "
              f"{point.initial_growth_rate:11.4e}  "
              f"{point.critical_crack_length * MM_PER_M:9.3f}  "
              f"{point.numerical_cycles:13.1f}")
    print("  Distinct from the sweep above: here sigma_min is held instead, so")
    print("  BOTH the growth rate and the fracture boundary move together.")

    # ------------------------------------------------------------------
    _rule("SENSITIVITY -- GEOMETRY FACTOR Y")
    print(f"  {'Y':>6}  {'da/dN(a0)':>11}  {'a_c [mm]':>9}  {'cycles':>13}")
    for point in geometry_factor_sensitivity(
        GEOMETRY_FACTORS, a0, cycle, law, toughness
    ):
        print(f"  {point.parameter_value:6.2f}  "
              f"{point.initial_growth_rate:11.4e}  "
              f"{point.critical_crack_length * MM_PER_M:9.3f}  "
              f"{point.numerical_cycles:13.1f}")
    print("  Y hurts twice: growth rate rises as Y**m while a_c falls as Y**-2.")

    # ------------------------------------------------------------------
    _rule("SENSITIVITY -- INITIAL FLAW SIZE")
    print(f"  {'a0 [mm]':>8}  {'a0/a_c':>8}  {'utilization':>12}  {'cycles':>14}")
    for point in initial_crack_sensitivity_to_fracture(
        INITIAL_CRACK_LENGTHS_M, cycle, geometry, law, toughness,
        intervals=WIDE_SPAN_INTERVALS,
    ):
        print(f"  {point.parameter_value * MM_PER_M:8.2f}  "
              f"{point.initial_crack_fraction:8.4f}  "
              f"{point.initial_utilization:12.4f}  "
              f"{point.numerical_cycles:14.1f}")
    print("  The assumed initial flaw dominates the answer far more than the")
    print("  toughness does.")

    # ------------------------------------------------------------------
    _rule("RESIDUAL-STRENGTH SCREEN (sigma_max = 120 MPa)")
    print(f"  {'a [mm]':>8}  {'K_max':>12}  {'sigma_res':>11}  "
          f"{'util':>7}  {'MS_K':>9}  {'result':>7}")
    for row in residual_strength_table(
        TABLE_CRACK_LENGTHS_M, cycle, geometry, toughness
    ):
        print(f"  {row.crack_length * MM_PER_M:8.2f}  "
              f"{row.k_max * MPA_PER_PA:8.4f} MPa.rm  "
              f"{row.residual_strength * MPA_PER_PA:7.1f} MPa  "
              f"{row.utilization:7.4f}  {row.margin:9.4f}  "
              f"{'PASS' if row.passes else 'FAIL':>7}")
    print(f"  The screen crosses from PASS to FAIL at a_c = {a_c * MM_PER_M:.3f} mm.")

    # ------------------------------------------------------------------
    _rule("INTERPRETATION")
    print("  * Crack growth accelerates with a, because delta_K is proportional")
    print("    to sqrt(a) and da/dN to delta_K**m.")
    print("  * The fracture endpoint is now derived from physics within LEFM")
    print("    rather than imposed: a_c is where K_max reaches K_IC.")
    print("  * Growth and fracture are driven by DIFFERENT quantities. At a_c,")
    print("    K_max = K_IC exactly, but delta_K = (delta_sigma/sigma_max)*K_IC,")
    print("    here 5/6 of K_IC. Using delta_K as a fracture criterion would")
    print("    overestimate the critical crack size substantially.")
    print("  * The result remains ILLUSTRATIVE: both the Paris curve and K_IC")
    print("    are illustrative inputs, not design allowables, and no specific")
    print("    alloy is claimed.")
    print("  * The toughness-derived critical crack size is an LEFM screening")
    print("    boundary, NOT a certified residual-strength allowable.")
    print("  * Still absent: delta_K threshold, crack closure, plastic-zone")
    print("    correction, net-section collapse, finite-width correction,")
    print("    R-ratio-dependent Paris coefficients, spectrum loading, scatter.")
    print()


if __name__ == "__main__":
    main()
