"""Milestone 5: variable-amplitude block-spectrum crack growth.

Run with::

    python examples/variable_amplitude_spectrum.py

All internal computation is SI; engineering units appear only in the printed
output. ``a`` is the HALF crack length throughout.

This is SEQUENTIAL PARIS GROWTH, not Miner's rule: no cumulative damage sum is
formed anywhere. The crack length is carried from block to block, and every
block is screened against its OWN threshold and its OWN fracture boundary.

The script runs a few hundred spectrum repeats and several sweeps, so it takes
of the order of a minute.
"""

from __future__ import annotations

import math

from crackgrowth import (
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GROWTH_THRESHOLD,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CANONICAL_PLATE_WIDTH,
    CANONICAL_SPECTRUM,
    CrackGrowthThreshold,
    LoadSpectrum,
    SpectrumBlock,
    StressCycle,
    block_boundary_table,
    block_count_sensitivity,
    cycles_to_finite_width_fracture,
    delta_k_for_geometry,
    initial_flaw_sensitivity_under_spectrum,
    mpa_sqrt_m_to_pa_sqrt_m,
    simulate_repeated_spectrum,
    threshold_sensitivity_under_spectrum,
)

MM = 1.0e3
MPA = 1.0e-6

SEVERE_COUNTS = (0, 10, 50, 100, 200)
THRESHOLDS_MPA = (None, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0)
INITIAL_FLAWS_M = (0.25e-3, 0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3)
ORDERS = ((0, 1, 2), (2, 1, 0), (1, 0, 2), (2, 0, 1))
SEQUENCE_REPEATS = 400

#: Sweeps report trends, not exact counts, so they use a looser block-solve
#: tolerance to keep the runtime reasonable. Headline numbers use the default.
#: A sweep row can therefore differ from the headline by up to ONE spectrum
#: (11,100 cycles, 0.18 % of life), because the answer is quantised to whole
#: blocks -- see the note under SENSITIVITY -- SEVERE-BLOCK COUNT.
SWEEP_TOL = 1.0e-9


def _rule(title: str) -> None:
    print()
    print(title)
    print("-" * len(title))


def _cycles(value: float) -> str:
    return "arrested (inf)" if math.isinf(value) else f"{value:,.0f}"


def _threshold(mpa: float) -> CrackGrowthThreshold:
    return CrackGrowthThreshold(
        name=f"{mpa} MPa*sqrt(m)",
        delta_k_threshold=mpa_sqrt_m_to_pa_sqrt_m(mpa),
        source_note=CANONICAL_GROWTH_THRESHOLD.source_note,
    )


def main() -> None:
    panel = CANONICAL_FINITE_WIDTH_GEOMETRY
    law = CANONICAL_PARIS_LAW
    toughness = CANONICAL_FRACTURE_TOUGHNESS
    threshold = CANONICAL_GROWTH_THRESHOLD
    spectrum = CANONICAL_SPECTRUM
    a0 = CANONICAL_INITIAL_CRACK_LENGTH

    print("=" * 78)
    print("VARIABLE-AMPLITUDE SPECTRUM CRACK GROWTH (MILESTONE 5)")
    print("=" * 78)

    # ------------------------------------------------------------------
    _rule("STUDY BASIS")
    print("  Crack convention      : a = HALF crack length; total crack = 2a")
    print("  Geometry              : centre through crack, finite-width panel")
    print(f"  Plate width W         : {CANONICAL_PLATE_WIDTH * MM:10.4f}  mm")
    print(f"  a0 (assumed flaw)     : {a0 * MM:10.4f}  mm")
    print(f"  Paris C  (SI basis)   : {law.C:10.4e}  m/cycle/(Pa*sqrt(m))**m")
    print(f"  Paris m               : {law.m:10.4f}  [-]")
    print("  K_IC                  : "
          f"{toughness.k_ic_in_mpa_sqrt_m:10.4f}  MPa*sqrt(m)")
    print("  delta_K_th            : "
          f"{threshold.delta_k_threshold_in_mpa_sqrt_m:10.4f}  MPa*sqrt(m)")
    print(f"  Cycles per spectrum   : {spectrum.cycles_per_spectrum:10,d}")
    print("  Method                : sequential Paris growth, block by block.")
    print("                          NOT Miner's rule -- no damage sum is formed.")
    print("                          Only the crack length carries between blocks.")
    print()
    for label, note in (
        ("PARIS", law.source_note),
        ("TOUGHNESS", toughness.source_note),
        ("THRESHOLD", threshold.source_note),
        ("SPECTRUM", spectrum.source_note),
    ):
        print(f"  {label} PROVENANCE")
        print(f"    {note}")
        print()

    # ------------------------------------------------------------------
    _rule("SPECTRUM (order is significant and is never sorted)")
    print(f"  {'block':>16} {'s_max':>7} {'s_min':>7} {'d_sig':>7} {'R':>7} "
          f"{'count':>8} {'dK(a0)':>9} {'state@a0':>10} {'a_th[mm]':>9} {'a_c[mm]':>9}")
    for block, row in zip(
        spectrum, block_boundary_table(spectrum, panel, toughness, threshold)
    ):
        dk = delta_k_for_geometry(a0, block.stress_cycle, panel)
        state = "ACTIVE" if threshold.is_active(dk) else "ARRESTED"
        print(f"  {block.name:>16} {block.sigma_max * MPA:>7.1f} "
              f"{block.sigma_min * MPA:>7.1f} {block.stress_range * MPA:>7.1f} "
              f"{block.stress_ratio:>7.4f} {block.cycle_count:>8,d} "
              f"{dk * MPA:>9.4f} {state:>10} "
              f"{row['threshold_crack_length'] * MM:>9.4f} "
              f"{row['critical_crack_length'] * MM:>9.4f}")
    print("  Each block has its OWN threshold size (from its delta_sigma) and its")
    print("  OWN fracture size (from its sigma_max). Neither is a spectrum average.")

    # ------------------------------------------------------------------
    result = simulate_repeated_spectrum(
        a0, spectrum, panel, law, toughness, threshold
    )

    _rule("FIRST SPECTRUM")
    print(f"  {'block':>16} {'a_start[mm]':>12} {'a_end[mm]':>11} "
          f"{'ext[um]':>9} {'cycles':>8} {'state':>10} {'util':>8}")
    for execution in result.history[: len(spectrum)]:
        a = execution.advance
        print(f"  {a.block_name:>16} {a.start_crack_length * MM:>12.6f} "
              f"{a.end_crack_length * MM:>11.6f} "
              f"{a.crack_extension * 1e6:>9.4f} {a.completed_cycles:>8,.0f} "
              f"{'ACTIVE' if a.active_at_start else 'ARRESTED':>10} "
              f"{a.start_fracture_utilization:>8.4f}")
    print("  The low block contributes nothing on this pass: its delta_K is below")
    print("  threshold at a 1 mm crack.")

    # ------------------------------------------------------------------
    _rule("LIFE RESULT")
    print(f"  Status                : {result.status.name}")
    print(f"  Full spectra completed: {result.completed_full_spectra:>16,d}")
    print(f"  Partial spectrum      : {result.partial_spectrum_cycles:>16,.0f}  cycles")
    print(f"  TOTAL CYCLES          : {result.total_cycles:>16,.1f}  cycles")
    print(f"  Final crack length    : {result.final_crack_length * MM:>16.6f}  mm")
    print(f"  Total crack extension : {result.total_crack_extension * MM:>16.6f}  mm")
    print(f"  Fracture block        : {result.fracture_block_name:>16} "
          f"(index {result.fracture_block_index}, repeat "
          f"{result.fracture_spectrum_repeat})")
    print(f"  Cycle within block    : {result.fracture_cycle_within_block:>16.4f}  "
          f"of {spectrum[result.fracture_block_index].cycle_count:,}")
    print(f"  K_max at fracture     : {result.final_k_max * MPA:>16.4f}  MPa*sqrt(m)")
    print(f"  Utilization           : {result.final_fracture_utilization:>16.6f}  [-]")
    print("  Min block a_c (envelope): "
          f"{result.minimum_block_critical_crack_length * MM:.6f} mm")
    print()
    print("  NOTE: fracture happens at cycle 0 of the severe block. The crack was")
    print("  pushed past the severe block's boundary by the two milder blocks")
    print("  earlier in the SAME pass -- neither of which fractures at that size,")
    print("  because their own boundaries are larger. The severe block then")
    print("  fractures on its very first cycle. This is exactly why fracture must")
    print("  be checked per block, in sequence, rather than against one envelope.")

    # ------------------------------------------------------------------
    _rule("GROWTH CONTRIBUTION BY BLOCK (crack extension, NOT Miner damage)")
    print(f"  {'block':>16} {'execs':>6} {'cycles':>12} {'ext[mm]':>10} "
          f"{'ext %':>7} {'1st active':>11} {'max dK':>8} {'fracture':>9}")
    for c in result.contributions:
        first = "-" if c.first_active_repeat is None else f"repeat {c.first_active_repeat}"
        print(f"  {c.block_name:>16} {c.executions:>6,d} {c.executed_cycles:>12,.0f} "
              f"{c.crack_extension * MM:>10.5f} {100 * c.extension_fraction:>7.2f} "
              f"{first:>11} {c.max_delta_k * MPA:>8.3f} "
              f"{str(c.triggered_fracture):>9}")
    print("  The low block starts ARRESTED, activates only at repeat 359, and yet")
    print("  ends up the LARGEST contributor to crack extension. It also runs 100x")
    print("  more cycles than the severe block but produces only ~4.6x the growth,")
    print("  so a cycle-weighted damage sum would badly misrank the blocks.")

    # ------------------------------------------------------------------
    _rule("COMPARISONS")
    no_threshold = simulate_repeated_spectrum(
        a0, spectrum, panel, law, toughness, None
    )
    print(f"  Threshold enabled     : {result.total_cycles:>16,.0f}  cycles, "
          f"{result.completed_full_spectra} spectra")
    print(f"  Threshold disabled    : {no_threshold.total_cycles:>16,.0f}  cycles, "
          f"{no_threshold.completed_full_spectra} spectra")
    print("  Disabling the threshold shortens the predicted life by "
          f"{100 * (1 - no_threshold.total_cycles / result.total_cycles):.1f} %,")
    print("  because the low block then grows the crack from the very first pass.")
    print()

    constant = LoadSpectrum(
        blocks=(
            SpectrumBlock(
                name="constant amplitude",
                stress_cycle=StressCycle(sigma_max=120.0e6, sigma_min=20.0e6),
                cycle_count=10_000,
            ),
        ),
        name="single-block regression",
        source_note="ILLUSTRATIVE regression spectrum - NOT FLIGHT LOAD DATA",
    )
    single = simulate_repeated_spectrum(
        a0, constant, panel, law, toughness, threshold
    )
    milestone4 = cycles_to_finite_width_fracture(
        a0, StressCycle(sigma_max=120.0e6, sigma_min=20.0e6), panel, law, toughness
    ).predicted_cycles
    print("  REGRESSION: a spectrum of one repeated block must reproduce the")
    print("  verified constant-amplitude answer.")
    print(f"    Single-block spectrum : {single.total_cycles:>16,.4f}  cycles")
    print(f"    Milestone 3/4 direct  : {milestone4:>16,.4f}  cycles")
    print("    Relative difference   : "
          f"{abs(single.total_cycles - milestone4) / milestone4:>16.3e}")

    # ------------------------------------------------------------------
    _rule(f"SEQUENCE ORDER (crack length after exactly {SEQUENCE_REPEATS} spectra)")
    print(f"  {'order':>12} {'final a [mm]':>14} {'vs first':>12}")
    baseline = None
    for order in ORDERS:
        crack = simulate_repeated_spectrum(
            a0,
            spectrum.reordered(order),
            panel,
            law,
            toughness,
            threshold,
            max_spectrum_repeats=SEQUENCE_REPEATS,
        ).final_crack_length
        if baseline is None:
            baseline = crack
        label = "->".join(spectrum.block_names[i][0] for i in order)
        print(f"  {label:>12} {crack * MM:>14.9f} "
              f"{100 * (crack / baseline - 1):>11.4f}%")
    print("  Order matters here ONLY through threshold activation: running the low")
    print("  block first catches it at the smallest crack of each pass, where it is")
    print("  most often still arrested, so it grows the crack less. Before the low")
    print("  block activates (~repeat 359) the orders agree to 1 part in 1e8.")
    print("  This is NOT overload retardation, residual stress or closure memory:")
    print("  the only state carried between blocks is the crack length.")

    # ------------------------------------------------------------------
    _rule("SENSITIVITY -- SEVERE-BLOCK COUNT")
    print(f"  {'count':>7} {'act@a0':>7} {'spectra':>8} {'cycles':>15} {'frac blk':>16}")
    for p in block_count_sensitivity(
        "C severe gust", SEVERE_COUNTS, a0, spectrum, panel, law, toughness,
        threshold, tolerance=SWEEP_TOL,
    ):
        print(f"  {int(p.parameter_value):>7,d} {p.active_blocks_at_start:>3}/"
              f"{p.total_blocks:<3} {p.completed_full_spectra:>8,d} "
              f"{_cycles(p.total_cycles):>15} {str(p.fracture_block_name):>16}")
    print("  Removing the severe block moves fracture to the manoeuvre block and")
    print("  lengthens the life: the governing boundary is then 19.4 mm, not 11.9.")
    print("  These rows use a looser solve tolerance, so the count=100 row may sit")
    print("  one spectrum (11,100 cycles, 0.18 %) away from the headline figure.")
    print("  That is the resolution of the model: the answer is quantised to whole")
    print("  blocks, so differences below one spectrum are not resolvable.")

    # ------------------------------------------------------------------
    _rule("SENSITIVITY -- THRESHOLD")
    print(f"  {'dK_th':>9} {'act@a0':>7} {'later':>6} {'spectra':>8} "
          f"{'cycles':>15} {'status':>22}")
    for p in threshold_sensitivity_under_spectrum(
        tuple(None if v is None else _threshold(v) for v in THRESHOLDS_MPA),
        a0, spectrum, panel, law, toughness, tolerance=SWEEP_TOL,
    ):
        label = (
            "disabled"
            if p.parameter_value == 0.0
            else f"{p.parameter_value * MPA:.1f}"
        )
        print(f"  {label:>9} {p.active_blocks_at_start:>3}/{p.total_blocks:<3} "
              f"{p.blocks_activating_later:>6} {p.completed_full_spectra:>8,d} "
              f"{_cycles(p.total_cycles):>15} {p.status.name:>22}")
    print("  The hard cutoff is discontinuous, so this is not a smooth trend: a")
    print("  threshold below every block's delta_K(a0) changes nothing at all,")
    print("  and a high enough one arrests the entire spectrum.")

    # ------------------------------------------------------------------
    _rule("SENSITIVITY -- INITIAL FLAW")
    print(f"  {'a0 [mm]':>8} {'act@a0':>7} {'later':>6} {'spectra':>8} "
          f"{'cycles':>15} {'status':>22}")
    for p in initial_flaw_sensitivity_under_spectrum(
        INITIAL_FLAWS_M, spectrum, panel, law, toughness, threshold,
        tolerance=SWEEP_TOL,
    ):
        print(f"  {p.parameter_value * MM:>8.2f} {p.active_blocks_at_start:>3}/"
              f"{p.total_blocks:<3} {p.blocks_activating_later:>6} "
              f"{p.completed_full_spectra:>8,d} {_cycles(p.total_cycles):>15} "
              f"{p.status.name:>22}")
    print("  A 0.25 mm flaw leaves EVERY block below threshold, so the spectrum is")
    print("  arrested outright. A 4 mm flaw activates all three immediately.")

    # ------------------------------------------------------------------
    _rule("INTERPRETATION")
    print("  * Low-amplitude cycles can be entirely arrested at first, then become")
    print("    growth-relevant later once severe cycles have grown the crack past")
    print("    their threshold size. Threshold state is recomputed on EVERY block")
    print("    execution, never cached from the initial flaw.")
    print("  * Fracture is governed by the K_max of the block actually being")
    print("    applied. A crack may safely exceed a severe block's boundary while")
    print("    milder blocks are running, and fracture the instant that severe")
    print("    block is next encountered.")
    print("  * This is sequential Paris growth, NOT Miner's rule. No damage sum is")
    print("    formed, and cycle counts do not rank blocks: the low block runs 100x")
    print("    the cycles of the severe block for ~4.6x the crack extension.")
    print("  * Sequence effects here come ONLY from evolving crack size, threshold")
    print("    activation and per-block fracture checks. No overload retardation,")
    print("    crack closure, residual stress or plasticity history is represented.")
    print("  * Blocks are piecewise constant; there is no rainflow extraction and")
    print("    no cycle-by-cycle variability within a block.")
    print("  * The spectrum, the Paris curve, K_IC and delta_K_th are ALL")
    print("    illustrative. This is not flight load data and no aircraft,")
    print("    certification spectrum or alloy is claimed.")
    print()


if __name__ == "__main__":
    main()
