"""Advancing a crack through one constant-amplitude load block.

Why not an explicit update
--------------------------
The crude update ``a_end = a_start + n * (da/dN at a_start)`` freezes the growth
rate across a whole block and can badly mis-predict growth when ``n`` is large or
the block is severe. It is not used anywhere in this package.

Instead the block advance INVERTS the already-verified life integral. For an
active, non-fracturing block the end crack length solves

    F(a_end) = N(a_start -> a_end) - n_block = 0

where ``N`` is the unchanged Milestone 3 geometry-aware log-grid integration.
``N`` is strictly increasing in ``a_end``, so the root is unique and a bounded
bisection finds it deterministically. This reuses verified machinery and carries
no explicit time-step error.

Per-block boundaries
--------------------
Every block has its OWN ``sigma_max`` and therefore its OWN fracture boundary
``a_c,i`` solved from ``K_max(a) = K_IC``, and its own threshold verdict from its
own ``delta_sigma``. Neither is computed once for the spectrum from an average
stress. A block that is arrested at one crack length may be active at a larger
one, so the threshold state is recomputed on every execution.

State carried between blocks
----------------------------
ONLY the crack length. No overload memory, no plastic-zone size, no closure
level, no residual stress, and no prior ``K_max`` history is retained.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from ._validation import require_positive
from .finite_width import FiniteWidthCenterCrack, delta_k_for_geometry
from .finite_width_fracture import (
    DEFAULT_ROOT_TOLERANCE,
    UPPER_BOUND_MARGIN,
    CriticalCrackStatus,
    finite_width_critical_crack_length,
    k_max_for_geometry,
)
from .fracture import FractureToughness
from .integration_log import integrate_crack_growth_log_grid
from .paris import ParisLaw
from .spectrum import SpectrumBlock
from .threshold import CrackGrowthThreshold

__all__ = [
    "DEFAULT_BLOCK_INTERVALS",
    "DEFAULT_BLOCK_MAX_ITERATIONS",
    "BlockAdvanceStatus",
    "BlockAdvanceResult",
    "advance_crack_through_block",
]

#: Default even interval count for the per-block log-grid integrations.
#: A block advances the crack only a little, so the integration span is short
#: and the log-grid rule is already extremely accurate there. Measured on the
#: canonical spectrum, 50, 100, 200 and 400 intervals give an identical total
#: cycle count and an identical final crack length to nine decimal places in mm;
#: the block-solve tolerance, not the interval count, sets the accuracy.
DEFAULT_BLOCK_INTERVALS = 50

#: Default bisection iteration cap for the block-end root solve.
DEFAULT_BLOCK_MAX_ITERATIONS = 200


class BlockAdvanceStatus(Enum):
    """Outcome of applying one block to a crack."""

    ADVANCED = "block completed; crack advanced"
    ARRESTED_BELOW_THRESHOLD = "block completed; delta_K at or below threshold"
    FRACTURED_WITHIN_BLOCK = "fracture reached partway through the block"
    ALREADY_AT_OR_BEYOND_FRACTURE = (
        "crack already at or beyond this block's fracture boundary"
    )
    GEOMETRY_LIMIT_REACHED = (
        "crack would exceed the admissible panel geometry within this block"
    )


@dataclass(frozen=True)
class BlockAdvanceResult:
    """Auditable record of one block execution."""

    block_name: str
    start_crack_length: float
    end_crack_length: float
    requested_cycles: int
    completed_cycles: float
    """Cycles actually applied. Equals ``requested_cycles`` for a completed or
    arrested block; less when fracture interrupts the block. Kept as a float:
    the integrated Paris prediction locates fracture at a generally non-integer
    cycle position, and it is not silently rounded."""

    start_geometry_factor: float
    end_geometry_factor: float
    start_delta_k: float
    end_delta_k: float
    delta_k_threshold: float | None
    active_at_start: bool
    active_at_end: bool
    block_critical_crack_length: float
    """This block's own ``a_c``; ``math.inf`` when the block is non-tensile."""

    start_fracture_utilization: float
    end_fracture_utilization: float
    fracture_occurred: bool
    fracture_cycle_within_block: float | None
    solver_iterations: int
    solver_tolerance: float
    status: BlockAdvanceStatus

    @property
    def crack_extension(self) -> float:
        """``a_end - a_start`` [m]."""
        return self.end_crack_length - self.start_crack_length


def _utilization(crack_length: float, block, geometry, toughness) -> float:
    k = k_max_for_geometry(crack_length, block.stress_cycle, geometry)
    return 0.0 if k <= 0.0 else k / toughness.k_ic


def advance_crack_through_block(
    a_start: float,
    block: SpectrumBlock,
    geometry: FiniteWidthCenterCrack,
    paris_law: ParisLaw,
    toughness: FractureToughness,
    threshold: CrackGrowthThreshold | None,
    intervals: int = DEFAULT_BLOCK_INTERVALS,
    tolerance: float = DEFAULT_ROOT_TOLERANCE,
    max_iterations: int = DEFAULT_BLOCK_MAX_ITERATIONS,
) -> BlockAdvanceResult:
    """Apply one block to a crack and report the resulting state.

    Parameters
    ----------
    threshold:
        Pass ``None`` to disable the threshold screen for this block, which
        reproduces the Milestone 3 no-threshold behaviour.
    tolerance:
        Absolute tolerance on the end crack length [m]. Exposed, not hidden.

    Notes
    -----
    An arrested block still CONSUMES its cycles -- they occur, they simply
    produce no growth -- so ``completed_cycles`` equals ``requested_cycles``
    and the crack extension is exactly zero.
    """
    if not isinstance(block, SpectrumBlock):
        raise TypeError(f"block must be a SpectrumBlock, got {block!r}")
    if not isinstance(geometry, FiniteWidthCenterCrack):
        raise TypeError(
            f"geometry must be a FiniteWidthCenterCrack, got {geometry!r}"
        )
    if not isinstance(paris_law, ParisLaw):
        raise TypeError(f"paris_law must be a ParisLaw, got {paris_law!r}")
    if threshold is not None and not isinstance(threshold, CrackGrowthThreshold):
        raise TypeError(
            f"threshold must be a CrackGrowthThreshold or None, got {threshold!r}"
        )
    tol = require_positive(tolerance, "tolerance")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int):
        raise TypeError(f"max_iterations must be an int, got {max_iterations!r}")
    if max_iterations < 1:
        raise ValueError(f"max_iterations must be at least 1, got {max_iterations!r}")
    a0 = require_positive(a_start, "a_start")
    # Rejects a0 >= W/2 with the finite-width admissibility message.
    start_geometry_factor = geometry.geometry_factor_at(a0)

    cycle = block.stress_cycle
    start_delta_k = delta_k_for_geometry(a0, cycle, geometry)
    threshold_value = None if threshold is None else threshold.delta_k_threshold
    active_at_start = (
        True if threshold is None else threshold.is_active(start_delta_k)
    )

    critical = finite_width_critical_crack_length(cycle, geometry, toughness)
    if critical.status is CriticalCrackStatus.ROOT_FOUND:
        a_c = critical.critical_crack_length
    elif critical.status is CriticalCrackStatus.NO_TENSILE_BOUNDARY:
        # A non-tensile block can still grow the crack algebraically (no closure
        # is modelled), but it can never fracture it in mode I.
        a_c = math.inf
    else:
        # Already critical at the search bound: treat as at/beyond fracture.
        a_c = 0.0

    def _result(
        a_end: float,
        completed: float,
        fracture: bool,
        fracture_cycle: float | None,
        iterations: int,
        status: BlockAdvanceStatus,
    ) -> BlockAdvanceResult:
        end_delta_k = delta_k_for_geometry(a_end, cycle, geometry)
        return BlockAdvanceResult(
            block_name=block.name,
            start_crack_length=a0,
            end_crack_length=a_end,
            requested_cycles=block.cycle_count,
            completed_cycles=completed,
            start_geometry_factor=start_geometry_factor,
            end_geometry_factor=geometry.geometry_factor_at(a_end),
            start_delta_k=start_delta_k,
            end_delta_k=end_delta_k,
            delta_k_threshold=threshold_value,
            active_at_start=active_at_start,
            active_at_end=(
                True if threshold is None else threshold.is_active(end_delta_k)
            ),
            block_critical_crack_length=a_c,
            start_fracture_utilization=_utilization(a0, block, geometry, toughness),
            end_fracture_utilization=_utilization(a_end, block, geometry, toughness),
            fracture_occurred=fracture,
            fracture_cycle_within_block=fracture_cycle,
            solver_iterations=iterations,
            solver_tolerance=tol,
            status=status,
        )

    # 1. Fracture is checked FIRST, at the very start of the block. A crack grown
    #    by an earlier block may already exceed THIS block's boundary, in which
    #    case the block fractures immediately without needing any growth of its
    #    own.
    if a0 >= a_c:
        return _result(
            a0, 0.0, True, 0.0, 0, BlockAdvanceStatus.ALREADY_AT_OR_BEYOND_FRACTURE
        )

    # 2. Threshold state, recomputed for THIS block at THIS crack length.
    if not active_at_start:
        return _result(
            a0,
            float(block.cycle_count),
            False,
            None,
            0,
            BlockAdvanceStatus.ARRESTED_BELOW_THRESHOLD,
        )

    # 3. Upper bracket: this block's fracture boundary, or the geometry limit
    #    when the block is non-tensile and can never fracture.
    a_upper = (
        a_c
        if math.isfinite(a_c)
        else 0.5 * geometry.plate_width * (1.0 - UPPER_BOUND_MARGIN)
    )

    def cycles_to(a_end: float) -> float:
        # The unchanged Milestone 3 integrator, which applies NO threshold.
        # That is correct here rather than convenient: delta_K increases
        # monotonically with a for this geometry, so a block that is active at
        # a_start stays active for every larger crack within the block. The
        # threshold can never re-bind mid-block, and using the thresholded rate
        # would change nothing while duplicating the integrator.
        return integrate_crack_growth_log_grid(
            a0, a_end, cycle, geometry, paris_law, intervals=intervals
        ).predicted_cycles

    cycles_to_upper = cycles_to(a_upper)

    # 4. Does the block run out of cycles before reaching the boundary?
    if cycles_to_upper <= block.cycle_count:
        if math.isfinite(a_c):
            return _result(
                a_upper,
                cycles_to_upper,
                True,
                cycles_to_upper,
                0,
                BlockAdvanceStatus.FRACTURED_WITHIN_BLOCK,
            )
        return _result(
            a_upper,
            cycles_to_upper,
            False,
            None,
            0,
            BlockAdvanceStatus.GEOMETRY_LIMIT_REACHED,
        )

    # 5. Otherwise invert the life integral for the end crack length.
    #    F(a) = N(a0 -> a) - n is strictly increasing, F(a0) = -n < 0 and
    #    F(a_upper) > 0, so the bracket is valid.
    #
    #    Seed a tighter bracket first. Because da/dN increases with a for this
    #    geometry, an explicit Euler step UNDER-predicts the growth, so
    #    a0 + n*(da/dN at a0) is a rigorous LOWER bound on a_end. The upper end
    #    is then found by doubling that step. Every trial point stays inside
    #    [a0, a_upper]: the bracket is refined, never widened beyond the
    #    admissible interval, and the fallback is the full original bracket.
    lo, hi = a0, a_upper
    euler_step = block.cycle_count * paris_law.growth_rate(start_delta_k)
    probes = 0
    if euler_step > 0.0:
        candidate = min(a0 + euler_step, a_upper)
        if a0 < candidate < a_upper:
            lo = candidate
            step = euler_step
            while probes < 64:
                step *= 2.0
                candidate = min(a0 + step, a_upper)
                probes += 1
                if candidate >= a_upper:
                    hi = a_upper
                    break
                if cycles_to(candidate) >= block.cycle_count:
                    hi = candidate
                    break
                lo = candidate

    iterations = probes
    while iterations < max_iterations and (hi - lo) > tol:
        mid = 0.5 * (lo + hi)
        if mid <= lo or mid >= hi:  # bracket collapsed to adjacent floats
            break
        iterations += 1
        if cycles_to(mid) < block.cycle_count:
            lo = mid
        else:
            hi = mid
    a_end = 0.5 * (lo + hi)
    return _result(
        a_end,
        float(block.cycle_count),
        False,
        None,
        iterations,
        BlockAdvanceStatus.ADVANCED,
    )
