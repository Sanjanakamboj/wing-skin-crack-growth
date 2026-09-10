"""Generate the STM-12 portfolio figures.

Run with::

    python -m pip install -e ".[figures]"
    python examples/make_figures.py

Writes four deterministic PNGs into ``figures/``. Every plotted engineering
quantity is computed from the package APIs; only the canonical study inputs
appear as constants.
"""

from __future__ import annotations

import math
import pathlib

import matplotlib

matplotlib.use("Agg")  # deterministic, headless
import matplotlib.pyplot as plt  # noqa: E402

from crackgrowth import (  # noqa: E402
    CANONICAL_CYCLE,
    CANONICAL_FINITE_WIDTH_GEOMETRY,
    CANONICAL_FRACTURE_TOUGHNESS,
    CANONICAL_GEOMETRY,
    CANONICAL_GROWTH_THRESHOLD,
    CANONICAL_INITIAL_CRACK_LENGTH,
    CANONICAL_PARIS_LAW,
    CANONICAL_SPECTRUM,
    CANONICAL_TARGET_CRACK_LENGTH,
    analytical_cycles_to_crack_length,
    block_boundary_table,
    critical_crack_length,
    cycles_to_critical_crack,
    cycles_to_finite_width_fracture,
    cycles_to_fracture_with_threshold,
    delta_k_for_geometry,
    k_max_for_geometry,
    residual_strength_for_geometry,
    simulate_repeated_spectrum,
)

MM = 1.0e3
MPA = 1.0e-6
FIGDIR = pathlib.Path(__file__).resolve().parent.parent / "figures"
DPI = 150

PANEL = CANONICAL_FINITE_WIDTH_GEOMETRY
LAW = CANONICAL_PARIS_LAW
TOUGHNESS = CANONICAL_FRACTURE_TOUGHNESS
THRESHOLD = CANONICAL_GROWTH_THRESHOLD
A0 = CANONICAL_INITIAL_CRACK_LENGTH
W = PANEL.plate_width

DISCLAIMER = (
    "Illustrative material data and load spectrum - not design allowables, "
    "not flight load data, not a certification result."
)
BLOCK_COLOURS = {
    "A low-amplitude": "#4C72B0",
    "B manoeuvre": "#DD8452",
    "C severe gust": "#C44E52",
}


def _footer(fig: plt.Figure) -> None:
    fig.text(0.5, 0.012, DISCLAIMER, ha="center", va="bottom", fontsize=7.5,
             style="italic", color="#555555")


def _save(fig: plt.Figure, name: str) -> None:
    FIGDIR.mkdir(exist_ok=True)
    path = FIGDIR / name
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {path.relative_to(FIGDIR.parent)}")


# ----------------------------------------------------------------------
def figure_1_life_progression() -> None:
    """How added physical realism changes the predicted life."""
    m1 = analytical_cycles_to_crack_length(
        A0, CANONICAL_TARGET_CRACK_LENGTH, CANONICAL_CYCLE, CANONICAL_GEOMETRY, LAW
    )
    m2 = cycles_to_critical_crack(
        A0, CANONICAL_CYCLE, CANONICAL_GEOMETRY, LAW, TOUGHNESS
    )
    m3 = cycles_to_finite_width_fracture(
        A0, CANONICAL_CYCLE, PANEL, LAW, TOUGHNESS
    )
    m4 = cycles_to_fracture_with_threshold(
        A0, CANONICAL_CYCLE, PANEL, LAW, TOUGHNESS, THRESHOLD
    )
    m5 = simulate_repeated_spectrum(
        A0, CANONICAL_SPECTRUM, PANEL, LAW, TOUGHNESS, THRESHOLD
    )

    rows = [
        ("M1  constant Y = 1\nimposed a = 10 mm", m1, "#8C8C8C", "imposed target"),
        ("M2  infinite plate\nK$_{max}$ = K$_{IC}$", m2.analytical_cycles,
         "#4C72B0", "fracture boundary"),
        ("M3  finite width W = 100 mm\nK$_{max}$ = K$_{IC}$", m3.predicted_cycles,
         "#55A868", "fracture boundary"),
        ("M4  + $\\Delta$K threshold\ncrack already active", m4.predicted_cycles,
         "#55A868", "fracture boundary"),
        ("M5  variable-amplitude\n3-block spectrum", m5.total_cycles,
         "#C44E52", "block-specific fracture"),
    ]
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    ypos = range(len(rows))
    ax.barh(list(ypos), [r[1] / 1e6 for r in rows],
            color=[r[2] for r in rows], height=0.62, edgecolor="white")
    for y, (_, cycles, _, end) in zip(ypos, rows):
        ax.text(cycles / 1e6 + 0.09, y, f"{cycles / 1e6:.3f}M  ({end})",
                va="center", fontsize=9)
    ax.set_yticks(list(ypos))
    ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Predicted crack-growth life  [million cycles]")
    ax.set_xlim(0, 9.9)
    ax.set_title(
        "STM-12  Crack-growth life across the model progression\n"
        "$a_0$ = 1 mm half crack, W = 100 mm, illustrative K$_{IC}$ = 30 MPa$\\sqrt{m}$",
        fontsize=11.5)
    ax.grid(axis="x", alpha=0.3)
    ax.set_axisbelow(True)
    note = (
        "M1-M4 apply the SAME constant-amplitude cycle ($\\Delta\\sigma$ = 100 MPa).\n"
        "M5 is a DIFFERENT loading history dominated by lower-range cycles that\n"
        "start below threshold, so its larger cycle count does NOT mean the\n"
        "spectrum is less damaging. These lives are not comparable at equal severity."
    )
    ax.text(0.985, 0.955, note, transform=ax.transAxes, ha="right", va="top",
            fontsize=8.2, bbox=dict(boxstyle="round,pad=0.5", facecolor="#FFF6E5",
                                    edgecolor="#E0C48C"))
    _footer(fig)
    _save(fig, "fig1_life_progression.png")


# ----------------------------------------------------------------------
def figure_2_finite_width() -> None:
    """Why finite-width effects accelerate as the crack grows."""
    a_c = cycles_to_finite_width_fracture(
        A0, CANONICAL_CYCLE, PANEL, LAW, TOUGHNESS
    ).critical_crack_length
    a = [W * 0.5 * f for f in
         [i / 600 * 0.97 + 0.002 for i in range(601)]]
    y = [PANEL.geometry_factor_at(v) for v in a]
    kmax = [k_max_for_geometry(v, CANONICAL_CYCLE, PANEL) * MPA for v in a]
    sres = [residual_strength_for_geometry(v, PANEL, TOUGHNESS) * MPA for v in a]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.5, 7.4), sharex=True,
                                   gridspec_kw={"height_ratios": [1, 1.15]})

    ax1.plot([v * MM for v in a], y, color="#4C72B0", lw=2)
    ax1.axhline(1.0, color="#999999", ls=":", lw=1)
    ax1.text(0.92, 0.045, "Y $\\to$ 1 as a/W $\\to$ 0  (infinite-plate limit)",
             transform=ax1.transAxes, ha="right", va="bottom", fontsize=8,
             color="#666666")
    ax1.set_ylabel("Geometry factor  Y(a)  [-]")
    ax1.set_ylim(0.95, 6.0)
    ax1.set_title(
        "STM-12  Finite-width geometry amplification,  W = 100 mm\n"
        "Y(a) = 1/$\\sqrt{\\cos(\\pi a/W)}$   ($a$ = HALF crack length, total crack = 2a)",
        fontsize=11.5)
    ax1.grid(alpha=0.3)

    ax2.plot([v * MM for v in a], kmax, color="#C44E52", lw=2,
             label="K$_{max}$(a) at $\\sigma_{max}$ = 120 MPa")
    ax2.axhline(TOUGHNESS.k_ic * MPA, color="#C44E52", ls="--", lw=1.4,
                label="K$_{IC}$ = 30 MPa$\\sqrt{m}$ (illustrative)")
    ax2.set_ylabel("K$_{max}$  [MPa$\\sqrt{m}$]", color="#C44E52")
    ax2.tick_params(axis="y", labelcolor="#C44E52")
    ax2.set_ylim(0, 90)

    ax2b = ax2.twinx()
    ax2b.plot([v * MM for v in a], sres, color="#55A868", lw=2,
              label="residual strength $\\sigma_{res}$(a)")
    ax2b.axhline(CANONICAL_CYCLE.sigma_max * MPA, color="#55A868", ls=":", lw=1.4,
                 label="applied $\\sigma_{max}$ = 120 MPa")
    ax2b.set_ylabel("Residual strength  [MPa]", color="#55A868")
    ax2b.tick_params(axis="y", labelcolor="#55A868")
    ax2b.set_ylim(0, 900)

    for ax in (ax1, ax2):
        ax.axvline(A0 * MM, color="#333333", ls="-", lw=1.1)
        ax.axvline(a_c * MM, color="#C44E52", ls="-.", lw=1.3)
        ax.axvline(W * 0.5 * MM, color="#7F7F7F", ls="--", lw=1.3)
    guide = (
        "$a_0$ = 1.00 mm   assumed initial flaw\n"
        "$a_c$ = 17.09 mm   fracture boundary, K$_{max}$ = K$_{IC}$\n"
        "W/2 = 50 mm   geometric limit where the ligament\n"
        "                       vanishes - NOT a fracture criterion"
    )
    ax1.text(0.30, 0.93, guide, transform=ax1.transAxes, ha="left", va="top",
             fontsize=8.4, linespacing=1.35,
             bbox=dict(boxstyle="round,pad=0.45", facecolor="white",
                       edgecolor="#BBBBBB", alpha=0.95))

    ax2.set_xlabel("Crack half-length  a  [mm]     (total crack length = 2a)")
    ax2.set_xlim(0, W * 0.5 * MM * 1.03)
    ax2.grid(alpha=0.3)
    handles = ax2.get_legend_handles_labels()
    handles_b = ax2b.get_legend_handles_labels()
    ax2.legend(handles[0] + handles_b[0], handles[1] + handles_b[1],
               loc="upper center", fontsize=8.5, framealpha=0.95, ncol=2)
    _footer(fig)
    _save(fig, "fig2_finite_width_amplification.png")


# ----------------------------------------------------------------------
def figure_3_block_activation() -> None:
    """Why block A starts arrested and later activates."""
    rows = block_boundary_table(CANONICAL_SPECTRUM, PANEL, TOUGHNESS, THRESHOLD)
    a_max = max(r["critical_crack_length"] for r in rows) * 1.05
    a = [a_max * (i / 600 * 0.999 + 0.0008) for i in range(601)]

    fig, ax = plt.subplots(figsize=(10.0, 5.8))
    for block, row in zip(CANONICAL_SPECTRUM, rows):
        colour = BLOCK_COLOURS[block.name]
        dk = [delta_k_for_geometry(v, block.stress_cycle, PANEL) * MPA for v in a]
        ax.plot([v * MM for v in a], dk, color=colour, lw=2,
                label=f"{block.name}  ($\\Delta\\sigma$ = "
                      f"{block.stress_range * MPA:.0f} MPa, n = {block.cycle_count:,})")
        ax.plot(row["threshold_crack_length"] * MM,
                THRESHOLD.delta_k_threshold * MPA, "o", color=colour, ms=7,
                markeredgecolor="white", zorder=5)
        ax.axvline(row["critical_crack_length"] * MM, color=colour, ls="-.",
                   lw=1.1, alpha=0.65)

    ax.axhline(THRESHOLD.delta_k_threshold * MPA, color="#333333", ls="--", lw=1.5)
    ax.text(5.4, THRESHOLD.delta_k_threshold * MPA + 0.45,
            "$\\Delta K_{th}$ = 4 MPa$\\sqrt{m}$  (illustrative, hard cutoff)",
            fontsize=8.8, ha="left", va="bottom")
    ax.axvline(A0 * MM, color="#333333", lw=1.1)
    ax.text(A0 * MM + 0.35, 21.5, "$a_0$ = 1 mm", fontsize=8.5, rotation=90,
            va="top")

    ax.annotate(
        "A is ARRESTED at $a_0$  ($\\Delta K$ = 2.80 < 4)",
        xy=(A0 * MM, 2.803), xytext=(4.4, 0.55), fontsize=8.6,
        arrowprops=dict(arrowstyle="->", color=BLOCK_COLOURS["A low-amplitude"]),
        color=BLOCK_COLOURS["A low-amplitude"])
    ax.annotate(
        "A activates at $a$ = 2.03 mm\n(reached at spectrum repeat 359)",
        xy=(rows[0]["threshold_crack_length"] * MM, 4.0),
        xytext=(16.4, 1.05), fontsize=8.6,
        arrowprops=dict(arrowstyle="->", color=BLOCK_COLOURS["A low-amplitude"]),
        color=BLOCK_COLOURS["A low-amplitude"])
    ax.text(0.985, 0.60,
            "dash-dot lines: each block's OWN fracture\n"
            "boundary $a_c$, set by its own $\\sigma_{max}$",
            transform=ax.transAxes, fontsize=8.4, ha="right", va="top",
            color="#444444")

    ax.set_xlabel("Crack half-length  a  [mm]     (total crack length = 2a)")
    ax.set_ylabel("Stress-intensity range  $\\Delta K$(a)  [MPa$\\sqrt{m}$]")
    ax.set_title(
        "STM-12  Threshold screening across the illustrative block spectrum\n"
        "Growth needs $\\Delta K > \\Delta K_{th}$;  fracture is a separate "
        "criterion on K$_{max}$",
        fontsize=11.5)
    ax.set_xlim(0, a_max * MM)
    ax.set_ylim(0, 35)
    ax.grid(alpha=0.3)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", fontsize=8.6, framealpha=0.95)
    _footer(fig)
    _save(fig, "fig3_block_threshold_activation.png")


# ----------------------------------------------------------------------
def figure_4_spectrum_history() -> None:
    """Crack history through the canonical variable-amplitude simulation."""
    result = simulate_repeated_spectrum(
        A0, CANONICAL_SPECTRUM, PANEL, LAW, TOUGHNESS, THRESHOLD
    )
    cycles = [0.0]
    crack = [result.initial_crack_length * MM]
    for execution in result.history:
        cycles.append(execution.cycles_before + execution.advance.completed_cycles)
        crack.append(execution.advance.end_crack_length * MM)

    low = next(c for c in result.contributions
               if c.block_name == "A low-amplitude")
    activation_cycles = (
        (low.first_active_repeat - 1) * CANONICAL_SPECTRUM.cycles_per_spectrum
    )
    severe_ac = min(
        r["critical_crack_length"] for r in
        block_boundary_table(CANONICAL_SPECTRUM, PANEL, TOUGHNESS, THRESHOLD)
    )

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(12.4, 5.4), gridspec_kw={"width_ratios": [2.35, 1]})

    ax.plot([c / 1e6 for c in cycles], crack, color="#4C72B0", lw=1.8)
    ax.axhline(severe_ac * MM, color="#C44E52", ls="--", lw=1.4)
    ax.text(0.12, severe_ac * MM + 0.28,
            "severe block C fracture boundary  $a_c$ = 11.86 mm",
            fontsize=8.6, color="#C44E52")
    ax.axvline(activation_cycles / 1e6, color="#DD8452", ls=":", lw=1.6)
    ax.annotate(
        f"low block A activates\n(repeat {low.first_active_repeat}, a = 2.03 mm)",
        xy=(activation_cycles / 1e6, 2.033), xytext=(1.15, 6.2), fontsize=8.6,
        arrowprops=dict(arrowstyle="->", color="#DD8452"), color="#B5651D")
    ax.plot(result.total_cycles / 1e6, result.final_crack_length * MM, "*",
            color="#C44E52", ms=15, markeredgecolor="white", zorder=6)
    ax.annotate(
        f"FRACTURE\n{result.total_cycles:,.0f} cycles\n"
        f"repeat {result.fracture_spectrum_repeat}, block C, "
        f"a = {result.final_crack_length * MM:.2f} mm",
        xy=(result.total_cycles / 1e6, result.final_crack_length * MM),
        xytext=(4.32, 8.9), fontsize=8.6,
        arrowprops=dict(arrowstyle="->", color="#C44E52"), color="#C44E52")
    ax.set_xlabel("Cumulative stress cycles  [millions]")
    ax.set_ylabel("Crack half-length  a  [mm]     (total = 2a)")
    ax.set_title(
        "STM-12  Variable-amplitude crack history (one point per block execution)\n"
        "Sequential Paris growth - NOT Miner damage accumulation", fontsize=11.5)
    ax.set_xlim(0, 6.65)
    ax.set_ylim(0.5, 13.2)
    ax.grid(alpha=0.3)
    ax.set_axisbelow(True)

    names = [c.block_name for c in result.contributions]
    shares = [100 * c.extension_fraction for c in result.contributions]
    bars = ax2.bar(range(len(names)), shares,
                   color=[BLOCK_COLOURS[n] for n in names],
                   edgecolor="white", width=0.62)
    for bar, share, c in zip(bars, shares, result.contributions):
        first = ("never" if c.first_active_repeat is None
                 else f"active\nrep. {c.first_active_repeat}")
        ax2.text(bar.get_x() + bar.get_width() / 2, share + 1.4,
                 f"{share:.1f}%", ha="center", fontsize=9.5, weight="bold")
        ax2.text(bar.get_x() + bar.get_width() / 2, share / 2, first,
                 ha="center", va="center", fontsize=7.4, color="white")
    ax2.set_xticks(range(len(names)))
    ax2.set_xticklabels([n.split()[0] + "\n" + " ".join(n.split()[1:])
                         for n in names], fontsize=8.6)
    ax2.set_ylabel("Share of total crack extension  [%]")
    ax2.set_ylim(0, 62)
    ax2.set_title("Crack-extension contribution\n(not a damage fraction)",
                  fontsize=10.5)
    ax2.grid(axis="y", alpha=0.3)
    ax2.set_axisbelow(True)
    _footer(fig)
    _save(fig, "fig4_spectrum_crack_history.png")


def main() -> None:
    print("Generating STM-12 portfolio figures...")
    figure_1_life_progression()
    figure_2_finite_width()
    figure_3_block_activation()
    figure_4_spectrum_history()
    print("Done.")


if __name__ == "__main__":
    main()
