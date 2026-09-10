"""Repeated variable-amplitude spectrum simulation.

The spectrum is executed block by block, in order, repeating until fracture,
until the whole spectrum is proven arrested, or until an explicit repeat guard
is hit. The only state carried forward is the crack length.

Life is reported in ACTUAL STRESS CYCLES. A spectrum repeat is a secondary
diagnostic. When fracture happens inside a block the total is

    cycles in all previously completed blocks
    + cycles completed inside the fracture block

and is never rounded up to a whole block or a whole spectrum. The fracture cycle
position is kept as a float, because the integrated Paris prediction is a
continuum estimate and generally locates fracture at a non-integer cycle.

**No Miner's rule.** Nothing here forms ``D = sum(n_i / N_i)``, and no life is
derived from a cumulative-damage sum. See :mod:`crackgrowth.spectrum`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from ._validation import require_positive
from .finite_width import FiniteWidthCenterCrack
from .finite_width_fracture import (
    DEFAULT_ROOT_TOLERANCE,
    CriticalCrackStatus,
    finite_width_critical_crack_length,
    k_max_for_geometry,
)
from .fracture import FractureToughness
from .paris import ParisLaw
from .spectrum import LoadSpectrum
from .spectrum_growth import (
    DEFAULT_BLOCK_INTERVALS,
    DEFAULT_BLOCK_MAX_ITERATIONS,
    BlockAdvanceResult,
    BlockAdvanceStatus,
    advance_crack_through_block,
)
from .threshold import CrackGrowthThreshold

__all__ = [
    "DEFAULT_MAX_SPECTRUM_REPEATS",
    "DEFAULT_MAX_HISTORY_RECORDS",
    "SpectrumStatus",
    "BlockExecution",
    "BlockContribution",
    "SpectrumCrackGrowthResult",
    "block_boundary_table",
    "minimum_block_critical_crack_length",
    "simulate_repeated_spectrum",
]

#: Default explicit guard on spectrum repeats. A numerical safety limit, NOT a
#: physical life prediction.
DEFAULT_MAX_SPECTRUM_REPEATS = 200_000

#: Default cap on stored per-block history records. Capping storage never
#: affects the computation or the growth-contribution accounting.
DEFAULT_MAX_HISTORY_RECORDS = 20_000


class SpectrumStatus(Enum):
    """How the spectrum simulation ended."""

    FRACTURE_REACHED = "fracture reached during a load block"
    SPECTRUM_ARRESTED = "every block arrested; the crack does not grow"
    MAX_REPEATS_REACHED = "repeat guard reached before fracture"
    INITIAL_FLAW_ALREADY_CRITICAL = (
        "the initial flaw is already at or beyond the first block's boundary"
    )
    NO_TENSILE_FRACTURE_BOUNDARY = (
        "no block has sigma_max > 0, so no tensile mode-I fracture boundary exists"
    )


@dataclass(frozen=True)
class BlockExecution:
    """One block execution, tagged with where it sat in the sequence."""

    spectrum_repeat: int
    """1-based index of the spectrum pass."""
    block_index: int
    """0-based index of the block within the spectrum."""
    cycles_before: float
    """Total cycles accumulated before this block started."""
    advance: BlockAdvanceResult


@dataclass(frozen=True)
class BlockContribution:
    """Growth contribution of one block, aggregated over all its executions.

    This is a CRACK-EXTENSION contribution, not Miner damage. Nothing here is a
    damage fraction and none of it is summed into a failure criterion.
    """

    block_name: str
    executions: int
    executed_cycles: float
    crack_extension: float
    extension_fraction: float
    """Share of the total crack extension [-], not a damage fraction."""
    first_active_repeat: int | None
    """1-based spectrum repeat in which this block first produced growth;
    ``None`` if it never did."""
    max_delta_k: float
    triggered_fracture: bool


@dataclass(frozen=True)
class SpectrumCrackGrowthResult:
    """Outcome of a repeated variable-amplitude spectrum simulation."""

    initial_crack_length: float
    final_crack_length: float
    total_cycles: float
    """Actual stress cycles to fracture (or applied before stopping)."""
    completed_full_spectra: int
    partial_spectrum_cycles: float
    """Cycles applied within the final, incomplete spectrum pass."""
    cycles_per_spectrum: int
    status: SpectrumStatus
    fracture_block_name: str | None
    fracture_block_index: int | None
    fracture_spectrum_repeat: int | None
    fracture_cycle_within_block: float | None
    minimum_block_critical_crack_length: float
    final_k_max: float
    """``K_max`` of the block in force at the end; ``-inf`` if undefined."""
    final_fracture_utilization: float
    history: tuple[BlockExecution, ...]
    history_truncated: bool
    contributions: tuple[BlockContribution, ...]

    @property
    def total_crack_extension(self) -> float:
        return self.final_crack_length - self.initial_crack_length

    @property
    def cycles_to_fracture(self) -> float:
        """Cycles to fracture, or ``math.inf`` when fracture is not reached."""
        if self.status is SpectrumStatus.FRACTURE_REACHED:
            return self.total_cycles
        if self.status is SpectrumStatus.INITIAL_FLAW_ALREADY_CRITICAL:
            return 0.0
        return math.inf


def minimum_block_critical_crack_length(
    spectrum: LoadSpectrum,
    geometry: FiniteWidthCenterCrack,
    toughness: FractureToughness,
) -> float:
    """Smallest ``a_c`` over the spectrum's tensile blocks [m].

    A DIAGNOSTIC ENVELOPE only. The simulation never terminates because the
    crack passed this value: fracture is checked in sequence, against the block
    actually being applied. A crack may legitimately exceed a severe block's
    ``a_c`` while that block is not being encountered.
    """
    sizes = [
        finite_width_critical_crack_length(
            block.stress_cycle, geometry, toughness
        ).critical_crack_length
        for block in spectrum
        if block.is_tensile
    ]
    return min(sizes) if sizes else math.inf


def block_boundary_table(
    spectrum: LoadSpectrum,
    geometry: FiniteWidthCenterCrack,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
) -> tuple[dict, ...]:
    """Per-block threshold and fracture boundaries, for reporting.

    Each block gets its own ``a_th`` (from its ``delta_sigma``) and its own
    ``a_c`` (from its ``sigma_max``); neither is a spectrum-wide average.
    """
    from .threshold import finite_width_threshold_crack_length

    rows = []
    for block in spectrum:
        critical = finite_width_critical_crack_length(
            block.stress_cycle, geometry, toughness
        )
        if threshold is None:
            a_th = None
        else:
            solved = finite_width_threshold_crack_length(
                block.stress_cycle, geometry, threshold
            )
            a_th = (
                solved.threshold_crack_length
                if solved.has_finite_boundary
                else math.inf
            )
        rows.append(
            {
                "name": block.name,
                "sigma_max": block.sigma_max,
                "sigma_min": block.sigma_min,
                "stress_range": block.stress_range,
                "cycle_count": block.cycle_count,
                "threshold_crack_length": a_th,
                "critical_crack_length": (
                    critical.critical_crack_length
                    if critical.status is CriticalCrackStatus.ROOT_FOUND
                    else math.inf
                ),
            }
        )
    return tuple(rows)


def simulate_repeated_spectrum(
    a_initial: float,
    spectrum: LoadSpectrum,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int = DEFAULT_BLOCK_INTERVALS,
    tolerance: float = DEFAULT_ROOT_TOLERANCE,
    max_iterations: int = DEFAULT_BLOCK_MAX_ITERATIONS,
    max_spectrum_repeats: int = DEFAULT_MAX_SPECTRUM_REPEATS,
    max_history_records: int = DEFAULT_MAX_HISTORY_RECORDS,
) -> SpectrumCrackGrowthResult:
    """Run the spectrum repeatedly until fracture, arrest, or the repeat guard.

    Parameters
    ----------
    threshold:
        ``None`` disables the threshold screen entirely.
    max_spectrum_repeats:
        Explicit numerical guard. Reaching it returns
        :attr:`SpectrumStatus.MAX_REPEATS_REACHED`; nothing is extrapolated
        beyond it, and it is not a physical life prediction.
    max_history_records:
        Cap on stored per-block records. Truncation affects reporting only --
        never the simulation or the contribution totals.
    """
    if not isinstance(spectrum, LoadSpectrum):
        raise TypeError(f"spectrum must be a LoadSpectrum, got {spectrum!r}")
    if not isinstance(geometry, FiniteWidthCenterCrack):
        raise TypeError(
            f"geometry must be a FiniteWidthCenterCrack, got {geometry!r}"
        )
    if not isinstance(paris_law, ParisLaw):
        raise TypeError(f"paris_law must be a ParisLaw, got {paris_law!r}")
    if not isinstance(toughness, FractureToughness):
        raise TypeError(f"toughness must be a FractureToughness, got {toughness!r}")
    for value, label in (
        (max_spectrum_repeats, "max_spectrum_repeats"),
        (max_history_records, "max_history_records"),
    ):
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError(f"{label} must be an int, got {value!r}")
        if value < 1:
            raise ValueError(f"{label} must be at least 1, got {value!r}")
    a0 = require_positive(a_initial, "a_initial")
    geometry.geometry_factor_at(a0)

    envelope = minimum_block_critical_crack_length(spectrum, geometry, toughness)

    # Aggregators are accumulated independently of the (capped) history.
    tallies: dict[str, dict] = {
        block.name: {
            "executions": 0,
            "cycles": 0.0,
            "extension": 0.0,
            "first_active_repeat": None,
            "max_delta_k": 0.0,
            "triggered_fracture": False,
        }
        for block in spectrum
    }

    def _finish(
        crack: float,
        total_cycles: float,
        completed_spectra: int,
        partial_cycles: float,
        status: SpectrumStatus,
        history: list[BlockExecution],
        truncated: bool,
        fracture_name: str | None = None,
        fracture_index: int | None = None,
        fracture_repeat: int | None = None,
        fracture_cycle: float | None = None,
        final_block=None,
    ) -> SpectrumCrackGrowthResult:
        total_extension = crack - a0
        contributions = tuple(
            BlockContribution(
                block_name=name,
                executions=data["executions"],
                executed_cycles=data["cycles"],
                crack_extension=data["extension"],
                extension_fraction=(
                    data["extension"] / total_extension
                    if total_extension > 0.0
                    else 0.0
                ),
                first_active_repeat=data["first_active_repeat"],
                max_delta_k=data["max_delta_k"],
                triggered_fracture=data["triggered_fracture"],
            )
            for name, data in tallies.items()
        )
        if final_block is None:
            k_max, utilization = -math.inf, 0.0
        else:
            k_max = k_max_for_geometry(crack, final_block.stress_cycle, geometry)
            utilization = 0.0 if k_max <= 0.0 else k_max / toughness.k_ic
        return SpectrumCrackGrowthResult(
            initial_crack_length=a0,
            final_crack_length=crack,
            total_cycles=total_cycles,
            completed_full_spectra=completed_spectra,
            partial_spectrum_cycles=partial_cycles,
            cycles_per_spectrum=spectrum.cycles_per_spectrum,
            status=status,
            fracture_block_name=fracture_name,
            fracture_block_index=fracture_index,
            fracture_spectrum_repeat=fracture_repeat,
            fracture_cycle_within_block=fracture_cycle,
            minimum_block_critical_crack_length=envelope,
            final_k_max=k_max,
            final_fracture_utilization=utilization,
            history=tuple(history),
            history_truncated=truncated,
            contributions=contributions,
        )

    if not spectrum.has_tensile_block:
        # No block can open the crack, so no tensile mode-I fracture boundary
        # exists anywhere in the spectrum. Algebraic delta_K from a
        # compression-only cycle must not be turned into a fracture prediction.
        return _finish(
            a0, 0.0, 0, 0.0, SpectrumStatus.NO_TENSILE_FRACTURE_BOUNDARY, [], False
        )

    crack = a0
    total_cycles = 0.0
    history: list[BlockExecution] = []
    truncated = False

    for repeat in range(1, max_spectrum_repeats + 1):
        spectrum_start_crack = crack
        cycles_at_spectrum_start = total_cycles

        for block_index, block in enumerate(spectrum):
            advance = advance_crack_through_block(
                crack,
                block,
                geometry,
                paris_law,
                toughness,
                threshold,
                intervals=intervals,
                tolerance=tolerance,
                max_iterations=max_iterations,
            )
            if len(history) < max_history_records:
                history.append(
                    BlockExecution(
                        spectrum_repeat=repeat,
                        block_index=block_index,
                        cycles_before=total_cycles,
                        advance=advance,
                    )
                )
            else:
                truncated = True

            data = tallies[block.name]
            data["executions"] += 1
            data["cycles"] += advance.completed_cycles
            data["extension"] += advance.crack_extension
            data["max_delta_k"] = max(data["max_delta_k"], advance.end_delta_k)
            if advance.crack_extension > 0.0 and data["first_active_repeat"] is None:
                data["first_active_repeat"] = repeat

            crack = advance.end_crack_length
            total_cycles += advance.completed_cycles

            if advance.fracture_occurred:
                data["triggered_fracture"] = True
                already = (
                    advance.status
                    is BlockAdvanceStatus.ALREADY_AT_OR_BEYOND_FRACTURE
                )
                status = (
                    SpectrumStatus.INITIAL_FLAW_ALREADY_CRITICAL
                    if already and total_cycles == 0.0
                    else SpectrumStatus.FRACTURE_REACHED
                )
                return _finish(
                    crack,
                    total_cycles,
                    repeat - 1,
                    total_cycles - cycles_at_spectrum_start,
                    status,
                    history,
                    truncated,
                    fracture_name=block.name,
                    fracture_index=block_index,
                    fracture_repeat=repeat,
                    fracture_cycle=advance.fracture_cycle_within_block,
                    final_block=block,
                )

            if advance.status is BlockAdvanceStatus.GEOMETRY_LIMIT_REACHED:
                return _finish(
                    crack,
                    total_cycles,
                    repeat - 1,
                    total_cycles - cycles_at_spectrum_start,
                    SpectrumStatus.MAX_REPEATS_REACHED,
                    history,
                    truncated,
                    final_block=block,
                )

        # A full pass that leaves the crack exactly where it started proves the
        # spectrum is arrested: every later identical pass does the same. There
        # is no need to loop to the repeat guard.
        if crack == spectrum_start_crack:
            return _finish(
                crack,
                total_cycles,
                repeat,
                0.0,
                SpectrumStatus.SPECTRUM_ARRESTED,
                history,
                truncated,
                final_block=spectrum[-1],
            )

    return _finish(
        crack,
        total_cycles,
        max_spectrum_repeats,
        0.0,
        SpectrumStatus.MAX_REPEATS_REACHED,
        history,
        truncated,
        final_block=spectrum[-1],
    )
