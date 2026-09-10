"""Sensitivity sweeps for repeated variable-amplitude spectrum growth.

Each sweep varies one quantity, states what is held fixed, and reports what the
model actually produces. Nothing here forms or reports a Miner damage sum.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from ._validation import require_positive
from .finite_width import FiniteWidthCenterCrack
from .finite_width_fracture import DEFAULT_ROOT_TOLERANCE
from .fracture import FractureToughness
from .loading import StressCycle
from .paris import ParisLaw
from .spectrum import LoadSpectrum, SpectrumBlock
from .spectrum_growth import DEFAULT_BLOCK_INTERVALS
from .spectrum_life import (
    DEFAULT_MAX_SPECTRUM_REPEATS,
    SpectrumStatus,
    simulate_repeated_spectrum,
)
from .threshold import CrackGrowthThreshold, finite_width_threshold_crack_length

__all__ = [
    "SpectrumPoint",
    "scaled_spectrum",
    "spectrum_with_block_count",
    "active_block_count",
    "spectrum_scale_sensitivity",
    "block_count_sensitivity",
    "threshold_sensitivity_under_spectrum",
    "toughness_sensitivity_under_spectrum",
    "width_sensitivity_under_spectrum",
    "initial_flaw_sensitivity_under_spectrum",
]


@dataclass(frozen=True)
class SpectrumPoint:
    """One swept variable-amplitude case."""

    parameter_name: str
    parameter_value: float
    initial_crack_length: float
    plate_width: float
    cycles_per_spectrum: int
    active_blocks_at_start: int
    total_blocks: int
    status: SpectrumStatus
    total_cycles: float
    completed_full_spectra: int
    final_crack_length: float
    fracture_block_name: str | None
    minimum_block_critical_crack_length: float
    blocks_activating_later: int
    """Blocks that were arrested at ``a0`` but produced growth at some point."""


def scaled_spectrum(spectrum: LoadSpectrum, scale: float) -> LoadSpectrum:
    """Scale every block's ``sigma_max`` and ``sigma_min`` about zero.

    Cycle counts and ordering are untouched, so only the stress levels change:
    ``delta_sigma``, ``K_max``, threshold activity and every block's fracture
    boundary all move together.
    """
    factor = require_positive(scale, "scale")
    return LoadSpectrum(
        blocks=tuple(
            SpectrumBlock(
                name=block.name,
                stress_cycle=StressCycle(
                    sigma_max=block.sigma_max * factor,
                    sigma_min=block.sigma_min * factor,
                ),
                cycle_count=block.cycle_count,
            )
            for block in spectrum
        ),
        name=f"{spectrum.name} (stress x{factor})",
        source_note=spectrum.source_note,
        notes=spectrum.notes,
    )


def spectrum_with_block_count(
    spectrum: LoadSpectrum, block_name: str, cycle_count: int
) -> LoadSpectrum:
    """Return the spectrum with one named block's cycle count replaced.

    A count of ``0`` REMOVES the block, since a
    :class:`~crackgrowth.spectrum.SpectrumBlock` must have at least one cycle.
    Removing the last remaining block is rejected rather than producing an
    empty spectrum.
    """
    if isinstance(cycle_count, bool) or not isinstance(cycle_count, int):
        raise TypeError(f"cycle_count must be an int, got {cycle_count!r}")
    if cycle_count < 0:
        raise ValueError(f"cycle_count must be non-negative, got {cycle_count!r}")
    names = spectrum.block_names
    if block_name not in names:
        raise ValueError(
            f"no block named {block_name!r} in this spectrum; have {names!r}"
        )
    blocks = []
    for block in spectrum:
        if block.name != block_name:
            blocks.append(block)
        elif cycle_count > 0:
            blocks.append(
                SpectrumBlock(
                    name=block.name,
                    stress_cycle=block.stress_cycle,
                    cycle_count=cycle_count,
                )
            )
    if not blocks:
        raise ValueError(
            "removing that block would leave an empty spectrum"
        )
    return LoadSpectrum(
        blocks=tuple(blocks),
        name=f"{spectrum.name} ({block_name} n={cycle_count})",
        source_note=spectrum.source_note,
        notes=spectrum.notes,
    )


def active_block_count(
    a_initial: float,
    spectrum: LoadSpectrum,
    geometry: FiniteWidthCenterCrack,
    threshold: CrackGrowthThreshold | None,
) -> int:
    """How many blocks are above threshold at a given crack length."""
    from .finite_width import delta_k_for_geometry

    if threshold is None:
        return len(spectrum)
    return sum(
        1
        for block in spectrum
        if threshold.is_active(
            delta_k_for_geometry(a_initial, block.stress_cycle, geometry)
        )
    )


def _point(
    parameter_name: str,
    parameter_value: float,
    a_initial: float,
    spectrum: LoadSpectrum,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int,
    tolerance: float,
    max_spectrum_repeats: int,
) -> SpectrumPoint:
    result = simulate_repeated_spectrum(
        a_initial,
        spectrum,
        geometry,
        paris_law,
        toughness,
        threshold,
        intervals=intervals,
        tolerance=tolerance,
        max_spectrum_repeats=max_spectrum_repeats,
    )
    from .finite_width import delta_k_for_geometry

    active_now = active_block_count(a_initial, spectrum, geometry, threshold)
    if threshold is None:
        # With the screen disabled every block is active from the outset, so
        # nothing can "activate later".
        initially_active = set(spectrum.block_names)
    else:
        initially_active = {
            block.name
            for block in spectrum
            if threshold.is_active(
                delta_k_for_geometry(a_initial, block.stress_cycle, geometry)
            )
        }
    later = sum(
        1
        for c in result.contributions
        if c.block_name not in initially_active and c.first_active_repeat is not None
    )
    return SpectrumPoint(
        parameter_name=parameter_name,
        parameter_value=float(parameter_value),
        initial_crack_length=a_initial,
        plate_width=geometry.plate_width,
        cycles_per_spectrum=spectrum.cycles_per_spectrum,
        active_blocks_at_start=active_now,
        total_blocks=len(spectrum),
        status=result.status,
        total_cycles=result.cycles_to_fracture,
        completed_full_spectra=result.completed_full_spectra,
        final_crack_length=result.final_crack_length,
        fracture_block_name=result.fracture_block_name,
        minimum_block_critical_crack_length=(
            result.minimum_block_critical_crack_length
        ),
        blocks_activating_later=later,
    )


def spectrum_scale_sensitivity(
    scales: Sequence[float],
    a_initial: float,
    spectrum: LoadSpectrum,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int = DEFAULT_BLOCK_INTERVALS,
    tolerance: float = DEFAULT_ROOT_TOLERANCE,
    max_spectrum_repeats: int = DEFAULT_MAX_SPECTRUM_REPEATS,
) -> tuple[SpectrumPoint, ...]:
    """Scale all block stresses about zero; counts and order unchanged."""
    return tuple(
        _point(
            "stress_scale",
            scale,
            a_initial,
            scaled_spectrum(spectrum, scale),
            geometry,
            paris_law,
            toughness,
            threshold,
            intervals,
            tolerance,
            max_spectrum_repeats,
        )
        for scale in scales
    )


def block_count_sensitivity(
    block_name: str,
    cycle_counts: Sequence[int],
    a_initial: float,
    spectrum: LoadSpectrum,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int = DEFAULT_BLOCK_INTERVALS,
    tolerance: float = DEFAULT_ROOT_TOLERANCE,
    max_spectrum_repeats: int = DEFAULT_MAX_SPECTRUM_REPEATS,
) -> tuple[SpectrumPoint, ...]:
    """Vary one block's cycle count; stress levels and order unchanged.

    A count of 0 removes the block. Note that adding cycles of an ARRESTED
    block adds cycles without adding growth, so total cycles to fracture can
    rise; once that block later activates the same cycles become
    growth-relevant. The trend is therefore not assumed either way.
    """
    return tuple(
        _point(
            f"{block_name} cycle_count",
            count,
            a_initial,
            spectrum_with_block_count(spectrum, block_name, count),
            geometry,
            paris_law,
            toughness,
            threshold,
            intervals,
            tolerance,
            max_spectrum_repeats,
        )
        for count in cycle_counts
    )


def threshold_sensitivity_under_spectrum(
    thresholds: Sequence[CrackGrowthThreshold | None],
    a_initial: float,
    spectrum: LoadSpectrum,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    intervals: int = DEFAULT_BLOCK_INTERVALS,
    tolerance: float = DEFAULT_ROOT_TOLERANCE,
    max_spectrum_repeats: int = DEFAULT_MAX_SPECTRUM_REPEATS,
) -> tuple[SpectrumPoint, ...]:
    """Vary ``delta_K_th``. ``None`` disables the screen for that row.

    The hard cutoff is discontinuous, so smooth behaviour is not assumed.
    """
    return tuple(
        _point(
            "delta_k_threshold",
            0.0 if threshold is None else threshold.delta_k_threshold,
            a_initial,
            spectrum,
            geometry,
            paris_law,
            toughness,
            threshold,
            intervals,
            tolerance,
            max_spectrum_repeats,
        )
        for threshold in thresholds
    )


def toughness_sensitivity_under_spectrum(
    k_ic_values: Sequence[float],
    a_initial: float,
    spectrum: LoadSpectrum,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int = DEFAULT_BLOCK_INTERVALS,
    tolerance: float = DEFAULT_ROOT_TOLERANCE,
    max_spectrum_repeats: int = DEFAULT_MAX_SPECTRUM_REPEATS,
) -> tuple[SpectrumPoint, ...]:
    """Vary ``K_IC``. Threshold activity at a fixed crack size cannot change,
    because ``K_IC`` does not enter ``delta_K``; only the per-block fracture
    boundaries move."""
    points = []
    for k_ic in k_ic_values:
        swept = FractureToughness(
            name=f"{toughness.name} (swept K_IC)",
            k_ic=k_ic,
            source_note=toughness.source_note,
            notes=toughness.notes,
        )
        points.append(
            _point(
                "k_ic",
                k_ic,
                a_initial,
                spectrum,
                geometry,
                paris_law,
                swept,
                threshold,
                intervals,
                tolerance,
                max_spectrum_repeats,
            )
        )
    return tuple(points)


def width_sensitivity_under_spectrum(
    plate_widths: Sequence[float],
    a_initial: float,
    spectrum: LoadSpectrum,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int = DEFAULT_BLOCK_INTERVALS,
    tolerance: float = DEFAULT_ROOT_TOLERANCE,
    max_spectrum_repeats: int = DEFAULT_MAX_SPECTRUM_REPEATS,
) -> tuple[SpectrumPoint, ...]:
    """Vary the panel width ``W``."""
    points = []
    for width in plate_widths:
        w = require_positive(width, "plate_width")
        if 2.0 * a_initial >= w:
            raise ValueError(
                f"initial crack is not admissible at W={w!r} m: require "
                f"2*a0 < W but 2*a0={2.0 * a_initial!r} m"
            )
        points.append(
            _point(
                "plate_width",
                w,
                a_initial,
                spectrum,
                FiniteWidthCenterCrack(plate_width=w),
                paris_law,
                toughness,
                threshold,
                intervals,
                tolerance,
                max_spectrum_repeats,
            )
        )
    return tuple(points)


def initial_flaw_sensitivity_under_spectrum(
    initial_crack_lengths: Sequence[float],
    spectrum: LoadSpectrum,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int = DEFAULT_BLOCK_INTERVALS,
    tolerance: float = DEFAULT_ROOT_TOLERANCE,
    max_spectrum_repeats: int = DEFAULT_MAX_SPECTRUM_REPEATS,
) -> tuple[SpectrumPoint, ...]:
    """Vary the assumed initial flaw. A larger flaw activates more blocks."""
    return tuple(
        _point(
            "initial_crack_length",
            a0,
            a0,
            spectrum,
            geometry,
            paris_law,
            toughness,
            threshold,
            intervals,
            tolerance,
            max_spectrum_repeats,
        )
        for a0 in initial_crack_lengths
    )


# Kept importable for callers building their own per-block diagnostics.
_ = (finite_width_threshold_crack_length, math)
