"""Variable-amplitude load-block spectrum representation.

A spectrum is an ORDERED sequence of constant-amplitude blocks. Order is
preserved exactly as supplied and is never sorted, re-ranked by severity, or
merged, because the crack length evolves as the spectrum is worked through: the
same block applied at a larger crack sees a larger ``Y(a)``, a larger
``delta_K``, a different threshold verdict and a nearer fracture boundary.

**This is not Miner's rule.** No cumulative-damage sum ``D = sum(n_i / N_i)`` is
formed anywhere in this package, and fracture is never predicted from such a
sum. Milestone 5 integrates the sequential Paris process directly, block by
block, carrying the crack length forward. A damage sum would discard exactly the
state that governs the answer.

Each block keeps its own :class:`~crackgrowth.loading.StressCycle`, so every
stress quantity comes from the verified Milestone 1 object and none of its
equations are restated here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Sequence

from .loading import StressCycle

__all__ = [
    "ILLUSTRATIVE_SPECTRUM_DISCLAIMER",
    "SpectrumBlock",
    "LoadSpectrum",
]

#: Mandatory provenance banner for any spectrum that is not measured load data.
ILLUSTRATIVE_SPECTRUM_DISCLAIMER = (
    "ILLUSTRATIVE VARIABLE-AMPLITUDE STRESS SPECTRUM - NOT FLIGHT LOAD DATA"
)


@dataclass(frozen=True)
class SpectrumBlock:
    """One constant-amplitude block within a spectrum.

    Parameters
    ----------
    name:
        Short identifier, used in reporting and in growth-contribution
        accounting. Must be non-empty.
    stress_cycle:
        The verified :class:`~crackgrowth.loading.StressCycle` for this block.
        No stress equation is duplicated here.
    cycle_count:
        Number of cycles in the block. Must be a positive ``int``; ``bool`` is
        rejected, matching the package convention for integer parameters
        (``intervals``, ``max_iterations``).
    """

    name: str
    stress_cycle: StressCycle
    cycle_count: int

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("name must be a non-empty string")
        if not isinstance(self.stress_cycle, StressCycle):
            raise TypeError(
                f"stress_cycle must be a StressCycle, got {self.stress_cycle!r}"
            )
        if isinstance(self.cycle_count, bool) or not isinstance(
            self.cycle_count, int
        ):
            raise TypeError(
                f"cycle_count must be an int, got {self.cycle_count!r}"
            )
        if self.cycle_count < 1:
            raise ValueError(
                f"cycle_count must be at least 1, got {self.cycle_count!r}"
            )

    # Convenience accessors. These DELEGATE to the stress cycle; they do not
    # re-derive anything.
    @property
    def sigma_max(self) -> float:
        """Maximum stress [Pa]."""
        return self.stress_cycle.sigma_max

    @property
    def sigma_min(self) -> float:
        """Minimum stress [Pa]."""
        return self.stress_cycle.sigma_min

    @property
    def stress_range(self) -> float:
        """``delta_sigma`` [Pa]. Algebraic range; no closure correction."""
        return self.stress_cycle.stress_range

    @property
    def stress_ratio(self) -> float:
        """``R = sigma_min / sigma_max`` [-]. Raises if ``sigma_max`` is zero."""
        return self.stress_cycle.stress_ratio

    @property
    def is_tensile(self) -> bool:
        """True when ``sigma_max > 0``, i.e. the block can open the crack."""
        return self.stress_cycle.sigma_max > 0.0


@dataclass(frozen=True)
class LoadSpectrum:
    """An ordered, repeatable sequence of load blocks.

    Parameters
    ----------
    blocks:
        One or more :class:`SpectrumBlock`, in the order they are applied. The
        order is preserved verbatim: no sorting, no severity ranking, no
        merging of like blocks.
    name:
        Short identifier for the spectrum.
    source_note:
        Provenance statement. For any spectrum that is not measured load data
        this must carry :data:`ILLUSTRATIVE_SPECTRUM_DISCLAIMER`.
    notes:
        Optional additional commentary.
    """

    blocks: tuple[SpectrumBlock, ...]
    name: str
    source_note: str
    notes: str = ""

    def __init__(
        self,
        blocks: Sequence[SpectrumBlock],
        name: str,
        source_note: str,
        notes: str = "",
    ) -> None:
        block_tuple = tuple(blocks)
        if not block_tuple:
            raise ValueError("a spectrum must contain at least one block")
        for block in block_tuple:
            if not isinstance(block, SpectrumBlock):
                raise TypeError(
                    f"every entry must be a SpectrumBlock, got {block!r}"
                )
        if not isinstance(name, str) or not name.strip():
            raise ValueError("name must be a non-empty string")
        if not isinstance(source_note, str) or not source_note.strip():
            raise ValueError("source_note must be a non-empty provenance statement")
        if not isinstance(notes, str):
            raise TypeError("notes must be a string")
        object.__setattr__(self, "blocks", block_tuple)
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "source_note", source_note)
        object.__setattr__(self, "notes", notes)

    def __len__(self) -> int:
        return len(self.blocks)

    def __iter__(self) -> Iterator[SpectrumBlock]:
        return iter(self.blocks)

    def __getitem__(self, index: int) -> SpectrumBlock:
        return self.blocks[index]

    @property
    def cycles_per_spectrum(self) -> int:
        """Total cycles in one full pass through the spectrum."""
        return sum(block.cycle_count for block in self.blocks)

    @property
    def block_names(self) -> tuple[str, ...]:
        """Block names in application order."""
        return tuple(block.name for block in self.blocks)

    @property
    def has_tensile_block(self) -> bool:
        """True when at least one block has ``sigma_max > 0``.

        When False no tensile mode-I fracture boundary exists anywhere in the
        spectrum, and no finite fracture life may be reported.
        """
        return any(block.is_tensile for block in self.blocks)

    def reordered(self, order: Sequence[int], name: str | None = None) -> "LoadSpectrum":
        """Return a new spectrum with the blocks in the given index order.

        Used only for the deliberate sequence-order study. It never mutates the
        original, and it is never applied automatically: a spectrum is otherwise
        always executed exactly as supplied.
        """
        indices = tuple(order)
        if sorted(indices) != list(range(len(self.blocks))):
            raise ValueError(
                "order must be a permutation of every block index, got "
                f"{indices!r} for {len(self.blocks)} blocks"
            )
        return LoadSpectrum(
            blocks=tuple(self.blocks[i] for i in indices),
            name=name or f"{self.name} (reordered {indices})",
            source_note=self.source_note,
            notes=self.notes,
        )
