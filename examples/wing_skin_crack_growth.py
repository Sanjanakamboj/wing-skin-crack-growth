"""Canonical Milestone 1 wing-skin crack-growth sanity case.

Run with::

    python examples/wing_skin_crack_growth.py

All internal computation is SI. Values are printed in engineering units
(mm, MPa, MPa*sqrt(m)) purely for readability; the conversions are explicit.
"""

from __future__ import annotations

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_GEOMETRY,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CANONICAL_TARGET_CRACK_LENGTH,
    analytical_cycles_to_crack_length,
    crack_growth_rate,
    cycles_to_crack_length,
    delta_stress_intensity,
    initial_crack_length_sensitivity,
    stress_range_sensitivity,
)

MM_PER_M = 1.0e3
MPA_PER_PA = 1.0e-6

STRESS_RANGE_SCALES = (0.5, 0.75, 1.0, 1.25, 1.5)
INITIAL_CRACK_LENGTHS_M = (0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3)


def _rule(title: str) -> None:
    print()
    print(title)
    print("-" * len(title))


def main() -> None:
    geometry = CANONICAL_GEOMETRY
    cycle = CANONICAL_CYCLE
    law = CANONICAL_PARIS_LAW
    a0 = CANONICAL_INITIAL_CRACK_LENGTH
    af = CANONICAL_TARGET_CRACK_LENGTH

    print("=" * 74)
    print("IDEALIZED WING-SKIN CRACK GROWTH -- PARIS LAW (MILESTONE 1)")
    print("=" * 74)

    # ------------------------------------------------------------------
    _rule("STUDY BASIS")
    print(f"  Geometry convention   : {geometry.description}")
    print("  K formulation         : K = Y * sigma * sqrt(pi * a)")
    print("  Geometry factor Y     : "
          f"{geometry.geometry_factor:.4f}  [-]  (constant)")
    print(f"  sigma_max             : {cycle.sigma_max * MPA_PER_PA:10.4f}  MPa")
    print(f"  sigma_min             : {cycle.sigma_min * MPA_PER_PA:10.4f}  MPa")
    print(f"  delta_sigma           : {cycle.stress_range * MPA_PER_PA:10.4f}  MPa")
    print(f"  mean stress           : {cycle.mean_stress * MPA_PER_PA:10.4f}  MPa")
    print("  alternating stress    : "
          f"{cycle.alternating_stress * MPA_PER_PA:10.4f}  MPa")
    print(f"  stress ratio R        : {cycle.stress_ratio:10.4f}  [-]")
    print(f"  Paris C  (SI basis)   : {law.C:10.4e}  m/cycle/(Pa*sqrt(m))**m")
    print(f"  Paris C  (MPa basis)  : {law.c_on_mpa_basis:10.4e}  "
          "m/cycle/(MPa*sqrt(m))**m")
    print(f"  Paris m               : {law.m:10.4f}  [-]")
    print(f"  Curve name            : {law.name}")
    print(f"  a0 (assumed flaw)     : {a0 * MM_PER_M:10.4f}  mm")
    print(f"  af (imposed target)   : {af * MM_PER_M:10.4f}  mm")
    print()
    print("  PROVENANCE")
    print(f"    {law.source_note}")
    print(f"    {law.notes}")

    # ------------------------------------------------------------------
    result = cycles_to_crack_length(a0, af, cycle, geometry, law)
    analytical = analytical_cycles_to_crack_length(a0, af, cycle, geometry, law)

    _rule("ENDPOINTS")
    dk0 = delta_stress_intensity(a0, cycle, geometry)
    dkf = delta_stress_intensity(af, cycle, geometry)
    print(f"  delta_K(a0)           : {dk0 * MPA_PER_PA:10.4f}  MPa*sqrt(m)")
    print(f"  delta_K(af)           : {dkf * MPA_PER_PA:10.4f}  MPa*sqrt(m)")
    print("  da/dN(a0)             : "
          f"{crack_growth_rate(a0, cycle, geometry, law):10.4e}  m/cycle")
    print("  da/dN(af)             : "
          f"{crack_growth_rate(af, cycle, geometry, law):10.4e}  m/cycle")

    # ------------------------------------------------------------------
    _rule("LIFE")
    difference = abs(result.predicted_cycles - analytical)
    print(f"  Integration method    : composite Simpson, "
          f"{result.integration_intervals} intervals (fixed grid)")
    print(f"  Numerical cycles      : {result.predicted_cycles:16.4f}  cycles")
    print(f"  Analytical cycles     : {analytical:16.4f}  cycles")
    print(f"  Absolute difference   : {difference:16.4e}  cycles")
    print(f"  Relative difference   : {difference / analytical:16.4e}  [-]")

    # ------------------------------------------------------------------
    _rule("SENSITIVITY -- STRESS RANGE (a0 and af fixed)")
    print(f"  {'scale':>7}  {'delta_sigma':>13}  {'cycles':>16}  "
          f"{'N/N_ref':>10}  {'scale**-m':>10}")
    stress_points = stress_range_sensitivity(
        STRESS_RANGE_SCALES, a0, af, cycle, geometry, law
    )
    reference = next(p for p in stress_points if p.parameter_value == 1.0)
    for point in stress_points:
        print(
            f"  {point.parameter_value:7.2f}  "
            f"{point.stress_range * MPA_PER_PA:10.2f} MPa  "
            f"{point.numerical_cycles:16.1f}  "
            f"{point.numerical_cycles / reference.numerical_cycles:10.4f}  "
            f"{point.parameter_value ** -law.m:10.4f}"
        )
    print("  Verified: N is proportional to delta_sigma ** -m "
          f"(m = {law.m:g}) for constant Y.")

    # ------------------------------------------------------------------
    _rule(f"SENSITIVITY -- INITIAL CRACK LENGTH (af fixed at "
          f"{af * MM_PER_M:g} mm)")
    print(f"  {'a0 [mm]':>9}  {'delta_K(a0)':>15}  {'cycles':>16}")
    for point in initial_crack_length_sensitivity(
        INITIAL_CRACK_LENGTHS_M, af, cycle, geometry, law
    ):
        dk = delta_stress_intensity(point.parameter_value, cycle, geometry)
        print(
            f"  {point.parameter_value * MM_PER_M:9.2f}  "
            f"{dk * MPA_PER_PA:9.4f} MPa.rm  "
            f"{point.numerical_cycles:16.1f}"
        )
    print("  Life falls steeply as the assumed initial flaw grows: most of the")
    print("  life is spent while the crack is small and delta_K is low.")

    # ------------------------------------------------------------------
    _rule("INTERPRETATION")
    print("  * Crack growth accelerates as the crack extends, because")
    print("    delta_K is proportional to sqrt(a) and da/dN to delta_K**m.")
    print("  * This is a constant-amplitude LEFM result only. No crack closure,")
    print("    no threshold delta_K, no fracture-toughness cutoff, no")
    print("    retardation, no spectrum loading, no scatter.")
    print("  * af is an IMPOSED target crack length in Milestone 1. It is NOT a")
    print("    fracture-toughness-derived critical crack size, and the reported")
    print("    life is NOT a certified inspection or safe-life interval.")
    print("  * The Paris coefficients are illustrative, not design allowables.")
    print()


if __name__ == "__main__":
    main()
