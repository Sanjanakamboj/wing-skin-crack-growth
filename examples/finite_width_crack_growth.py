"""Milestone 3: crack growth in a finite-width centre-cracked wing skin.

Run with::

    python examples/finite_width_crack_growth.py

All internal computation is SI. Engineering units (mm, MPa, MPa*sqrt(m)) appear
only in the printed output, where the conversion is explicit.

Throughout, ``a`` is the HALF crack length; the total tip-to-tip crack length is
``2a``. This example supersedes the infinite-plate endpoint of
``examples/critical_crack_growth.py``, which is left unchanged and still runs.
"""

from __future__ import annotations

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GEOMETRY,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CANONICAL_PLATE_WIDTH,
    cycles_to_critical_crack,
    cycles_to_finite_width_fracture,
    finite_width_max_stress_sensitivity,
    finite_width_stress_range_sensitivity,
    frozen_geometry_factor_life,
    geometry_amplification_table,
    initial_crack_sensitivity_finite_width,
    mpa_sqrt_m_to_pa_sqrt_m,
    residual_strength_for_geometry,
    toughness_sensitivity_finite_width,
    width_sensitivity,
)

MM_PER_M = 1.0e3
MPA_PER_PA = 1.0e-6

WIDTHS_M = (40.0e-3, 50.0e-3, 75.0e-3, 100.0e-3, 150.0e-3, 200.0e-3, 500.0e-3)
INITIAL_CRACKS_M = (0.25e-3, 0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3, 8.0e-3)
SIGMA_MAX_FIXED_RANGE_MPA = (100.0, 120.0, 140.0, 160.0, 180.0)
SIGMA_MAX_FIXED_MIN_MPA = (80.0, 100.0, 120.0, 140.0, 160.0, 180.0)
K_IC_SWEEP_MPA = (20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0)
AMPLIFICATION_CRACKS_M = (
    0.5e-3, 1.0e-3, 2.0e-3, 5.0e-3, 10.0e-3, 15.0e-3, 20.0e-3, 30.0e-3, 40.0e-3,
)
CONVERGENCE_INTERVALS = (50, 100, 200, 500, 1000, 2000, 4000)


def _rule(title: str) -> None:
    print()
    print(title)
    print("-" * len(title))


def main() -> None:
    panel = CANONICAL_FINITE_WIDTH_GEOMETRY
    cycle = CANONICAL_CYCLE
    law = CANONICAL_PARIS_LAW
    toughness = CANONICAL_FRACTURE_TOUGHNESS
    a0 = CANONICAL_INITIAL_CRACK_LENGTH
    width = CANONICAL_PLATE_WIDTH

    print("=" * 78)
    print("FINITE-WIDTH WING-SKIN CRACK GROWTH (MILESTONE 3)")
    print("=" * 78)

    # ------------------------------------------------------------------
    _rule("STUDY BASIS")
    print("  Crack convention      : a = HALF crack length; total crack = 2a")
    print("  Geometry              : centre through crack in a finite-width panel")
    print("  Geometry factor       : Y(a) = 1 / sqrt(cos(pi * a / W))")
    print(f"  Plate width W         : {width * MM_PER_M:10.4f}  mm")
    print(f"  a0 (assumed flaw)     : {a0 * MM_PER_M:10.4f}  mm  "
          f"(total 2a0 = {2 * a0 * MM_PER_M:.4f} mm)")
    print(f"  sigma_max             : {cycle.sigma_max * MPA_PER_PA:10.4f}  MPa")
    print(f"  sigma_min             : {cycle.sigma_min * MPA_PER_PA:10.4f}  MPa")
    print(f"  delta_sigma           : {cycle.stress_range * MPA_PER_PA:10.4f}  MPa")
    print(f"  stress ratio R        : {cycle.stress_ratio:10.4f}  [-]")
    print(f"  Paris C  (SI basis)   : {law.C:10.4e}  m/cycle/(Pa*sqrt(m))**m")
    print(f"  Paris m               : {law.m:10.4f}  [-]")
    print("  K_IC                  : "
          f"{toughness.k_ic_in_mpa_sqrt_m:10.4f}  MPa*sqrt(m)")
    print("  Integration           : composite Simpson on a log(a) grid,")
    print("                          Y re-evaluated at every abscissa")
    print("  Critical-size solver  : bounded bisection on K_max(a) - K_IC")
    print()
    print("  PARIS PROVENANCE")
    print(f"    {law.source_note}")
    print()
    print("  TOUGHNESS PROVENANCE")
    print(f"    {toughness.source_note}")

    # ------------------------------------------------------------------
    result = cycles_to_finite_width_fracture(a0, cycle, panel, law, toughness)
    critical = result.critical
    a_c = result.critical_crack_length

    _rule("FINITE-WIDTH GEOMETRY AT THE INITIAL FLAW")
    print(f"  a0 / W                : {a0 / width:10.6f}  [-]")
    print(f"  Y(a0)                 : {result.initial_geometry_factor:10.6f}  [-]")
    print("  Remaining ligament    : "
          f"{panel.remaining_ligament(a0) * MM_PER_M:10.4f}  mm")
    print(f"  Ligament fraction     : {panel.ligament_fraction(a0):10.6f}  [-]"
          "   (diagnostic only)")

    # ------------------------------------------------------------------
    _rule("FRACTURE BOUNDARY")
    print(f"  Solver status         : {critical.status.value}")
    print(f"  Bisection iterations  : {critical.iterations:10d}")
    print(f"  Tolerance on a_c      : {critical.tolerance:10.2e}  m")
    print(f"  Critical half crack   : {a_c * MM_PER_M:10.4f}  mm")
    print("  Total crack 2a_c      : "
          f"{critical.total_crack_length_at_critical * MM_PER_M:10.4f}  mm")
    print(f"  Y(a_c)                : {critical.geometry_factor_at_critical:10.6f}  [-]")
    print("  K_max(a_c)            : "
          f"{critical.k_max_at_critical * MPA_PER_PA:10.4f}  MPa*sqrt(m)  (= K_IC)")
    print("  delta_K(a_c)          : "
          f"{result.delta_k_at_critical * MPA_PER_PA:10.4f}  MPa*sqrt(m)")
    print("  delta_K(a_c)/K_IC     : "
          f"{result.delta_k_at_critical / toughness.k_ic:10.6f}  [-]  "
          f"(= delta_sigma/sigma_max = {cycle.stress_range / cycle.sigma_max:.6f})")
    print("  Residual ligament     : "
          f"{critical.residual_ligament * MM_PER_M:10.4f}  mm")
    print(f"  Ligament fraction     : {critical.ligament_fraction:10.6f}  [-]"
          "   (diagnostic only)")
    print("  Residual strength(a_c): "
          f"{residual_strength_for_geometry(a_c, panel, toughness) * MPA_PER_PA:10.4f}"
          "  MPa   (= sigma_max)")

    # ------------------------------------------------------------------
    infinite = cycles_to_critical_crack(
        a0, cycle, CANONICAL_GEOMETRY, law, toughness
    )
    growth = result.growth_result

    _rule("LIFE TO FRACTURE")
    print(f"  Integration method    : {result.integration_method}")
    print(f"  Intervals             : {growth.integration_intervals:10d}")
    print(f"  Finite-width life     : {result.predicted_cycles:16.4f}  cycles")
    print(f"  Infinite-plate life   : {infinite.analytical_cycles:16.4f}  cycles")
    print("  Life reduction        : "
          f"{100 * (1 - result.predicted_cycles / infinite.analytical_cycles):16.2f}  %")
    print("  Critical-size red.    : "
          f"{100 * (1 - a_c / infinite.critical_crack_length):16.2f}  %")
    print(f"  da/dN(a0)             : {growth.initial_growth_rate:16.4e}  m/cycle")
    print(f"  da/dN(a_c)            : {growth.final_growth_rate:16.4e}  m/cycle")
    print("  Rate ratio            : "
          f"{growth.final_growth_rate / growth.initial_growth_rate:16.2f}  [-]")
    print("  The critical size moves far more than the life does: most cycles")
    print("  are spent while the crack is small and Y is still close to 1.")

    # ------------------------------------------------------------------
    _rule("FROZEN-Y DIAGNOSTIC (approximation quality, NOT a result)")
    frozen = frozen_geometry_factor_life(a0, a_c, cycle, panel, law)
    print(f"  True Y(a) integration : {result.predicted_cycles:16.4f}  cycles")
    print(f"  Y frozen at Y(a0)     : {frozen:16.4f}  cycles")
    print("  Overestimate          : "
          f"{100 * (frozen / result.predicted_cycles - 1):16.2f}  %")
    print("  Freezing Y at the initial flaw ignores the amplification that")
    print("  builds up as the crack grows, so it overstates the life. Shown")
    print("  only to quantify that error; never used as a result.")

    # ------------------------------------------------------------------
    _rule("NUMERICAL CONVERGENCE (log-grid Simpson)")
    reference = cycles_to_finite_width_fracture(
        a0, cycle, panel, law, toughness, intervals=100000
    ).predicted_cycles
    print(f"  {'n':>6}  {'life [cycles]':>18}  {'rel. vs finest':>15}")
    for n in CONVERGENCE_INTERVALS:
        life = cycles_to_finite_width_fracture(
            a0, cycle, panel, law, toughness, intervals=n
        ).predicted_cycles
        print(f"  {n:>6}  {life:>18.6f}  {abs(life - reference) / reference:>15.3e}")
    print("  Measured error ratio is ~16 per interval doubling over n = 50-400,")
    print("  i.e. fourth order, as expected for Simpson; beyond n ~ 1000 the")
    print("  result sits at the floating-point floor.")

    # ------------------------------------------------------------------
    _rule("WIDTH SENSITIVITY (a0 and cycle fixed)")
    print(f"  {'W [mm]':>8}  {'Y(a0)':>9}  {'a_c [mm]':>9}  {'2a_c/W':>8}  "
          f"{'lig.frac':>9}  {'cycles':>13}")
    for point in width_sensitivity(WIDTHS_M, a0, cycle, law, toughness):
        print(f"  {point.plate_width * MM_PER_M:8.1f}  "
              f"{point.initial_geometry_factor:9.6f}  "
              f"{point.critical_crack_length * MM_PER_M:9.4f}  "
              f"{2 * point.critical_crack_length / point.plate_width:8.4f}  "
              f"{point.ligament_fraction_at_critical:9.6f}  "
              f"{point.predicted_cycles:13.1f}")
    print(f"  Infinite plate (M2): a_c = "
          f"{infinite.critical_crack_length * MM_PER_M:.4f} mm, "
          f"life = {infinite.analytical_cycles:.1f} cycles.")
    print("  Wider panels converge on that reference from below.")

    # ------------------------------------------------------------------
    _rule("INITIAL-CRACK SENSITIVITY (W and cycle fixed)")
    print(f"  {'a0 [mm]':>8}  {'a0/W':>8}  {'Y(a0)':>9}  {'utilization':>12}  "
          f"{'cycles':>13}")
    for point in initial_crack_sensitivity_finite_width(
        INITIAL_CRACKS_M, cycle, panel, law, toughness
    ):
        print(f"  {point.parameter_value * MM_PER_M:8.2f}  "
              f"{point.initial_crack_width_ratio:8.4f}  "
              f"{point.initial_geometry_factor:9.6f}  "
              f"{point.initial_utilization:12.4f}  "
              f"{point.predicted_cycles:13.1f}")

    # ------------------------------------------------------------------
    _rule("SENSITIVITY -- sigma_max AT FIXED delta_sigma = 100 MPa")
    print(f"  {'s_max':>7}  {'R':>8}  {'dK(a0)':>9}  {'da/dN(a0)':>11}  "
          f"{'a_c [mm]':>9}  {'cycles':>13}")
    for point in finite_width_max_stress_sensitivity(
        tuple(v * 1.0e6 for v in SIGMA_MAX_FIXED_RANGE_MPA),
        100.0e6, a0, panel, law, toughness,
    ):
        s_min = point.sigma_max - point.stress_range
        print(f"  {point.sigma_max * MPA_PER_PA:7.1f}  "
              f"{s_min / point.sigma_max:8.4f}  "
              f"{point.initial_delta_k * MPA_PER_PA:9.4f}  "
              f"{point.initial_growth_rate:11.4e}  "
              f"{point.critical_crack_length * MM_PER_M:9.4f}  "
              f"{point.predicted_cycles:13.1f}")
    print("  delta_K(a0) and da/dN(a0) are IDENTICAL down the column: delta_K")
    print("  depends only on delta_sigma. The life falls purely because K_max")
    print("  rises and pulls the fracture boundary inward.")

    # ------------------------------------------------------------------
    _rule("SENSITIVITY -- sigma_max AT FIXED sigma_min = 20 MPa")
    print(f"  {'s_max':>7}  {'d_sigma':>8}  {'dK(a0)':>9}  {'da/dN(a0)':>11}  "
          f"{'a_c [mm]':>9}  {'cycles':>13}")
    for point in finite_width_stress_range_sensitivity(
        tuple(v * 1.0e6 for v in SIGMA_MAX_FIXED_MIN_MPA),
        20.0e6, a0, panel, law, toughness,
    ):
        print(f"  {point.sigma_max * MPA_PER_PA:7.1f}  "
              f"{point.stress_range * MPA_PER_PA:8.1f}  "
              f"{point.initial_delta_k * MPA_PER_PA:9.4f}  "
              f"{point.initial_growth_rate:11.4e}  "
              f"{point.critical_crack_length * MM_PER_M:9.4f}  "
              f"{point.predicted_cycles:13.1f}")
    print("  Distinct from the sweep above: here BOTH the growth driving force")
    print("  and the fracture boundary move.")

    # ------------------------------------------------------------------
    _rule("TOUGHNESS SENSITIVITY (finite width)")
    print(f"  {'K_IC':>7}  {'a_c [mm]':>9}  {'Y(a_c)':>8}  {'lig.frac':>9}  "
          f"{'cycles':>13}  {'a_c/K_IC^2':>11}")
    toughness_points = toughness_sensitivity_finite_width(
        tuple(mpa_sqrt_m_to_pa_sqrt_m(v) for v in K_IC_SWEEP_MPA),
        a0, cycle, panel, law, toughness,
    )
    base_ratio = (
        toughness_points[0].critical_crack_length / toughness_points[0].k_ic**2
    )
    for point in toughness_points:
        ratio = point.critical_crack_length / point.k_ic**2
        print(f"  {point.k_ic * MPA_PER_PA:7.1f}  "
              f"{point.critical_crack_length * MM_PER_M:9.4f}  "
              f"{point.geometry_factor_at_critical:8.4f}  "
              f"{point.ligament_fraction_at_critical:9.6f}  "
              f"{point.predicted_cycles:13.1f}  {ratio / base_ratio:11.4f}")
    print("  The last column is a_c/K_IC**2 normalised to the first row. For the")
    print("  infinite plate it would be exactly 1.0000 throughout; here it falls")
    print("  to 0.46, so the Milestone 2 law a_c ~ K_IC**2 does NOT survive")
    print("  finite width. Extra toughness buys progressively less crack length,")
    print("  because a larger a_c sits at a higher Y.")

    # ------------------------------------------------------------------
    _rule("GEOMETRY AMPLIFICATION AT W = 100 mm")
    print(f"  {'a [mm]':>8}  {'a/W':>7}  {'2a [mm]':>8}  {'lig.frac':>9}  "
          f"{'Y(a)':>8}  {'K_max':>9}  {'delta_K':>9}")
    for row in geometry_amplification_table(
        AMPLIFICATION_CRACKS_M, cycle, panel
    ):
        print(f"  {row.crack_length * MM_PER_M:8.2f}  "
              f"{row.crack_width_ratio:7.4f}  "
              f"{row.total_crack_length * MM_PER_M:8.2f}  "
              f"{row.ligament_fraction:9.6f}  "
              f"{row.geometry_factor:8.4f}  "
              f"{row.k_max * MPA_PER_PA:7.4f}  "
              f"{row.delta_k * MPA_PER_PA:9.4f}")
    print("  Y is nearly 1 for small cracks and climbs steeply past a/W ~ 0.3.")

    # ------------------------------------------------------------------
    _rule("INTERPRETATION")
    print("  * Finite width amplifies K by Y(a), and Y grows as the crack grows,")
    print("    so the driving force rises faster than the infinite-plate model.")
    print("  * Growth therefore accelerates more sharply than Milestone 2 says,")
    print("    and the toughness boundary is reached at a smaller crack.")
    print("  * The critical size falls 14 % while the life falls only 4 %: most")
    print("    of the life is spent at small cracks where Y is still near 1.")
    print("  * Y(a) must be evaluated inside the integral. Freezing it at Y(a0)")
    print("    overestimates the life by ~2 % here, and more for narrow panels.")
    print("  * delta_K(a_c)/K_IC = delta_sigma/sigma_max still holds exactly,")
    print("    because Y(a_c) multiplies both K_max and delta_K and cancels.")
    print("  * a_c ~ K_IC**2 does NOT hold for finite width.")
    print("  * This remains pure LEFM. The ligament fraction is reported as a")
    print("    DIAGNOSTIC only: there is no net-section-collapse model, no")
    print("    plastic-zone validity check, and no thickness correction. For a")
    print("    short ligament, collapse rather than fracture would realistically")
    print("    govern, and this model would not say so.")
    print("  * The finite-width critical crack length remains an LEFM screening")
    print("    boundary, NOT a certified residual-strength allowable.")
    print()


if __name__ == "__main__":
    main()
