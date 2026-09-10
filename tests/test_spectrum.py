"""A-F, BF-BH. Load-block and spectrum representation."""

from __future__ import annotations

import dataclasses

import pytest

from crackgrowth import (
    CANONICAL_CYCLE,
    CANONICAL_SPECTRUM,
    ILLUSTRATIVE_SPECTRUM_DISCLAIMER,
    LoadSpectrum,
    SpectrumBlock,
    StressCycle,
)


def _block(name="A", smax=100.0e6, smin=10.0e6, n=100):
    return SpectrumBlock(
        name=name,
        stress_cycle=StressCycle(sigma_max=smax, sigma_min=smin),
        cycle_count=n,
    )


def _spectrum(*blocks):
    return LoadSpectrum(
        blocks=blocks or (_block(),),
        name="test",
        source_note="ILLUSTRATIVE (test fixture)",
    )


# --------------------------------------------------------------------------
# A. SpectrumBlock validation
# --------------------------------------------------------------------------
def test_block_stores_its_cycle_and_count():
    block = _block(n=250)
    assert block.cycle_count == 250
    assert block.stress_cycle is not None


@pytest.mark.parametrize("bad", [0, -1, -100])
def test_non_positive_cycle_count_rejected(bad):
    with pytest.raises(ValueError, match="at least 1"):
        _block(n=bad)


@pytest.mark.parametrize("bad", [1.0, "10", None])
def test_non_integer_cycle_count_rejected(bad):
    with pytest.raises(TypeError, match="cycle_count must be an int"):
        _block(n=bad)


def test_bool_cycle_count_rejected():
    """Package convention: bool is not an acceptable int."""
    with pytest.raises(TypeError, match="cycle_count must be an int"):
        _block(n=True)


def test_empty_block_name_rejected():
    with pytest.raises(ValueError, match="name"):
        _block(name="   ")


def test_non_stress_cycle_rejected():
    with pytest.raises(TypeError, match="stress_cycle must be a StressCycle"):
        SpectrumBlock(name="x", stress_cycle="cycle", cycle_count=10)


def test_block_is_immutable():
    with pytest.raises(dataclasses.FrozenInstanceError):
        _block().cycle_count = 5


# --------------------------------------------------------------------------
# F. convenience values delegate to the verified StressCycle
# --------------------------------------------------------------------------
def test_block_convenience_values_match_the_stress_cycle():
    cycle = StressCycle(sigma_max=120.0e6, sigma_min=20.0e6)
    block = SpectrumBlock(name="b", stress_cycle=cycle, cycle_count=7)
    assert block.sigma_max == cycle.sigma_max
    assert block.sigma_min == cycle.sigma_min
    assert block.stress_range == cycle.stress_range
    assert block.stress_ratio == cycle.stress_ratio
    assert block.stress_range == pytest.approx(100.0e6)
    assert block.is_tensile


def test_non_tensile_block_is_flagged():
    block = _block(smax=-10.0e6, smin=-100.0e6)
    assert not block.is_tensile
    assert block.stress_range == pytest.approx(90.0e6)  # algebraic, no closure


# --------------------------------------------------------------------------
# B. LoadSpectrum validation
# --------------------------------------------------------------------------
def test_empty_spectrum_rejected():
    with pytest.raises(ValueError, match="at least one block"):
        LoadSpectrum(blocks=(), name="x", source_note="ILLUSTRATIVE")


def test_non_block_entry_rejected():
    with pytest.raises(TypeError, match="must be a SpectrumBlock"):
        LoadSpectrum(blocks=("nope",), name="x", source_note="ILLUSTRATIVE")


def test_empty_name_or_provenance_rejected():
    with pytest.raises(ValueError, match="name"):
        LoadSpectrum(blocks=(_block(),), name=" ", source_note="ILLUSTRATIVE")
    with pytest.raises(ValueError, match="source_note"):
        LoadSpectrum(blocks=(_block(),), name="x", source_note="")


def test_spectrum_is_immutable():
    with pytest.raises(dataclasses.FrozenInstanceError):
        _spectrum().name = "other"


def test_spectrum_accepts_a_list_and_stores_a_tuple():
    spectrum = LoadSpectrum(
        blocks=[_block("A"), _block("B")], name="x", source_note="ILLUSTRATIVE"
    )
    assert isinstance(spectrum.blocks, tuple)


# --------------------------------------------------------------------------
# C/E. order is preserved and never sorted
# --------------------------------------------------------------------------
def test_spectrum_preserves_the_given_order():
    spectrum = _spectrum(
        _block("severe", smax=200.0e6),
        _block("low", smax=60.0e6),
        _block("medium", smax=120.0e6),
    )
    assert spectrum.block_names == ("severe", "low", "medium")
    assert [b.name for b in spectrum] == ["severe", "low", "medium"]
    assert spectrum[0].name == "severe"


def test_spectrum_does_not_sort_by_severity():
    """A severity-sorted spectrum would put 'low' first. It must not."""
    spectrum = _spectrum(
        _block("severe", smax=200.0e6),
        _block("low", smax=60.0e6),
    )
    assert spectrum.block_names[0] == "severe"
    assert spectrum[0].sigma_max > spectrum[1].sigma_max


def test_spectrum_does_not_merge_like_blocks():
    spectrum = _spectrum(_block("A", n=10), _block("A-again", n=10))
    assert len(spectrum) == 2
    assert spectrum.cycles_per_spectrum == 20


# --------------------------------------------------------------------------
# D. cycle accounting
# --------------------------------------------------------------------------
def test_cycles_per_spectrum_sums_the_blocks():
    spectrum = _spectrum(_block("A", n=10_000), _block("B", n=1_000), _block("C", n=100))
    assert spectrum.cycles_per_spectrum == 11_100
    assert len(spectrum) == 3


def test_has_tensile_block():
    assert _spectrum(_block(smax=100.0e6)).has_tensile_block
    assert not _spectrum(_block(smax=-10.0e6, smin=-100.0e6)).has_tensile_block
    mixed = _spectrum(_block("t", smax=100.0e6), _block("c", smax=-10.0e6, smin=-90.0e6))
    assert mixed.has_tensile_block


# --------------------------------------------------------------------------
# reordering is explicit and non-mutating
# --------------------------------------------------------------------------
def test_reordered_returns_a_new_spectrum_and_leaves_the_original_alone():
    spectrum = _spectrum(_block("A"), _block("B"), _block("C"))
    reordered = spectrum.reordered((2, 1, 0))
    assert reordered.block_names == ("C", "B", "A")
    assert spectrum.block_names == ("A", "B", "C")
    assert reordered.cycles_per_spectrum == spectrum.cycles_per_spectrum


@pytest.mark.parametrize("bad", [(0, 1), (0, 0, 1), (0, 1, 3)])
def test_reordered_rejects_a_non_permutation(bad):
    with pytest.raises(ValueError, match="permutation"):
        _spectrum(_block("A"), _block("B"), _block("C")).reordered(bad)


# --------------------------------------------------------------------------
# canonical spectrum
# --------------------------------------------------------------------------
def test_canonical_spectrum_shape_and_provenance():
    assert CANONICAL_SPECTRUM.block_names == (
        "A low-amplitude",
        "B manoeuvre",
        "C severe gust",
    )
    counts = [b.cycle_count for b in CANONICAL_SPECTRUM]
    assert counts == [10_000, 1_000, 100]
    assert CANONICAL_SPECTRUM.cycles_per_spectrum == 11_100
    assert ILLUSTRATIVE_SPECTRUM_DISCLAIMER in CANONICAL_SPECTRUM.source_note
    assert CANONICAL_SPECTRUM.has_tensile_block


def test_canonical_spectrum_stress_levels():
    levels = [(b.sigma_max, b.sigma_min) for b in CANONICAL_SPECTRUM]
    assert levels == [(70.0e6, 20.0e6), (110.0e6, 20.0e6), (150.0e6, 20.0e6)]
    ranges = [b.stress_range for b in CANONICAL_SPECTRUM]
    assert ranges == pytest.approx([50.0e6, 90.0e6, 130.0e6])


def test_canonical_cycle_is_not_repointed_by_milestone5():
    assert CANONICAL_CYCLE.sigma_max == 120.0e6
    assert CANONICAL_CYCLE.sigma_min == 20.0e6
