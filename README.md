# wing-skin-crack-growth

Verified fatigue-crack-growth analysis for an idealized aircraft wing-skin crack
under cyclic loading, built on linear-elastic fracture mechanics (LEFM) and the
Paris law.

Python package: `crackgrowth` (src layout, no runtime dependencies).

## Objective

Predict the number of constant-amplitude load cycles required to grow a
through-thickness crack in a wing skin from an assumed initial flaw size to a
defined final crack length, with the numerical result independently verified
against a closed-form reference solution.

## Scope (Milestone 1)

Implemented:

- crack geometry with a constant geometry factor `Y`
- constant-amplitude cyclic stress representation
- mode-I stress-intensity factor `K`
- stress-intensity range `ΔK`
- Paris-law material representation with explicit unit handling
- crack-growth rate `da/dN`
- deterministic numerical integration from `a₀` to a specified `a_f`
- closed-form analytical reference for constant `Y`
- stress-range and initial-crack-size sensitivity sweeps

Not implemented — see [Limitations](#limitations).

## Units and sign conventions

All internal computation is SI. Nothing is stored in mixed units.

| Quantity | Symbol | Unit |
| --- | --- | --- |
| Crack length | `a` | m |
| Stress | `σ` | Pa |
| Stress intensity | `K`, `ΔK` | Pa·√m |
| Paris coefficient | `C` | m/cycle / (Pa·√m)^m |
| Paris exponent | `m` | – |
| Crack-growth rate | `da/dN` | m/cycle |
| Life | `N` | cycles |

Conventions:

- Tensile stress is **positive**; compressive stress is **negative**.
- Mode-I crack opening is driven by tensile stress.
- Crack length must be finite and strictly positive.
- `a_final > a_initial` is required; equal or reversed limits are rejected.
- A cycle requires `σ_min ≤ σ_max`.
- Engineering units (mm, MPa, MPa·√m) appear only in printed output, where the
  conversion is explicit.

## Governing equations

**Mode-I stress intensity** (constant geometry factor):

```
K = Y · σ · sqrt(π · a)
```

`K` is returned **signed**: positive in tension, negative in compression. No
absolute value is applied silently. A negative `K` is a bookkeeping result, not
an opening-mode driving force.

**Stress-intensity range:**

```
ΔK = Y · Δσ · sqrt(π · a),    Δσ = σ_max − σ_min
```

`Δσ` is the raw **algebraic** range. Milestone 1 does **not** model crack
closure, so a compressive `σ_min` still contributes in full to `Δσ` and
therefore inflates `ΔK`. Real aluminium behaviour would partially or wholly
remove the compressive part of the cycle. No effective-`ΔK` model is offered.

**Paris law:**

```
da/dN = C · (ΔK)^m
```

with `C > 0`, `m > 0`, both finite.

## Geometry convention

Milestone 1 supports one deliberately simple geometry: `ThroughCrackGeometry`,
an idealized through crack characterised by a single constant dimensionless
geometry factor `Y`.

The canonical value `Y = 1.0` corresponds to a central through crack of length
`2a` in an **infinite** plate. There is no finite-width correction, no surface
or corner-crack shape factor, and no crack-length-dependent `Y(a)`. `Y` is a
constant supplied by the analyst and is an explicit engineering assumption.

## Paris-law convention and the C-unit conversion

`ParisLaw` always stores `C` on the **SI basis**, i.e. with `ΔK` in Pa·√m and
`da/dN` in m/cycle. Because `m` is dimensionless, `C` carries units of
`m/cycle / (Pa·√m)^m` and its numerical value depends on `m`.

Handbook coefficients are almost always quoted on the MPa·√m basis:

```
da/dN [m/cycle] = C_MPa · (ΔK [MPa·√m])^m
```

Substituting `ΔK[MPa·√m] = ΔK[Pa·√m] / 10⁶`:

```
da/dN = C_MPa · (ΔK[Pa·√m] / 10⁶)^m
      = (C_MPa / (10⁶)^m) · (ΔK[Pa·√m])^m
```

hence the exact conversion

```
C_SI = C_MPa / (10⁶)^m
```

implemented as `paris_c_from_mpa_basis(c_mpa, m)`, with the inverse
`paris_c_to_mpa_basis` and the constructor `paris_law_from_mpa_basis`.

**The package never guesses the basis.** A `ParisLaw` value is always
interpreted as SI; conversion from the MPa basis must be requested by name.
Mixing the two bases is an error of `(10⁶)^m` — a factor of 10¹⁸ at `m = 3` —
and is covered explicitly by test.

For the canonical curve, `C_MPa = 1.0e-11` and `m = 3.0` give
`C_SI = 1.0e-11 / 10¹⁸ = 1.0e-29`.

## Paris-law data and provenance

> **ILLUSTRATIVE PARIS-LAW INPUT — NOT DESIGN ALLOWABLE**

No measured fatigue-crack-growth dataset was obtained or verified for this
milestone. The canonical curve is **illustrative**:

| Parameter | Value |
| --- | --- |
| `C` (MPa·√m basis) | 1.0 × 10⁻¹¹ m/cycle / (MPa·√m)^m |
| `C` (SI basis) | 1.0 × 10⁻²⁹ m/cycle / (Pa·√m)^m |
| `m` | 3.0 |

These coefficients were selected **after** auditing the resulting physics over
the intended crack-size range, not by tuning toward a preferred life:

- `ΔK` spans 5.60 → 17.72 MPa·√m over `a = 1 → 10 mm` at `Δσ = 100 MPa` —
  comfortably above a typical threshold and well below a typical aluminium
  toughness, so the Paris form is defensible over the whole range;
- `da/dN` spans 1.76 × 10⁻⁹ → 5.57 × 10⁻⁸ m/cycle — a plausible mid-Paris band;
- the resulting life is ≈ 7.8 × 10⁵ cycles — large enough to be meaningful,
  small enough to inspect numerically.

**No specific aerospace alloy is claimed.** These values are not measured, not
traceable to a qualification dataset, and must not be used as design allowables.

## Analytical reference

For constant `Y` and constant `Δσ`, `da/dN = C·(Y·Δσ·√π)^m · a^(m/2)`, so

```
dN/da = 1 / [ C · (Y · Δσ · √π)^m · a^(m/2) ]
```

Integrating from `a₀` to `a_f`, with `p = 1 − m/2`:

```
m ≠ 2:   N = (a_f^p − a₀^p) / [ p · C · (Y · Δσ · √π)^m ]

m = 2:   N = ln(a_f / a₀) / [ C · (Y · Δσ · √π)² ]
```

**Sign note for `m > 2`:** then `p < 0`, so `a_f^p < a₀^p` and the numerator is
negative, while the denominator carries the same negative `p`. The quotient is
positive, as a life must be. This is verified by an explicit test at
`m = 2.5, 3.0, 3.5, 4.0, 5.0`.

The `m = 2` branch is additionally checked against the `m → 2` limit of the
power branch.

## Numerical integration method

```
N = ∫[a₀ → a_f] da / (da/dN)
```

evaluated with a **fixed-grid composite Simpson rule**:

- the analyst supplies the interval count (default 1000); there is no adaptive
  refinement and no hidden tolerance;
- Simpson's rule consumes sub-intervals in pairs, so an **even** interval count
  is required — an odd count is **rejected**, not silently incremented;
- the grid spans exactly `[a₀, a_f]` and never evaluates outside it;
- the result is bit-for-bit reproducible across repeated calls.

SciPy is not required: the integrand is smooth and one-dimensional, so a
transparent hand-written Simpson rule over `math` is both sufficient and easier
to audit. The package therefore has **no runtime dependencies**.

The quadrature evaluates `da/dN` through the full modelling chain
(`ΔK` → Paris law) rather than through the analytical solution's factored
constant, so the numerical and analytical routines are genuinely independent and
their agreement is a real cross-check.

**Zero-range policy:** if `Δσ = 0` then `ΔK = 0` and `da/dN = 0`. The crack never
reaches the target, and both the numerical and analytical routines return
`math.inf` — an explicit, documented no-growth result. Nothing divides by zero.

## Verification strategy

Verification rests on four independent legs:

1. **Hand calculations.** Literal expected values, derived by hand and written
   into the tests, for `K`, `ΔK`, `R`, `da/dN`, the `C` unit conversion, and the
   analytical life on both the `m ≠ 2` and `m = 2` branches.
2. **Scaling laws.** `K ∝ σ`, `K ∝ √a`, `ΔK ∝ Δσ`, `ΔK ∝ Y`, `ΔK ∝ √a`,
   `da/dN ∝ C`, `da/dN ∝ ΔK^m`, `N ∝ Δσ^−m`, `N ∝ Y^−m`, `N ∝ C^−1`.
3. **Numerical vs analytical agreement.** Cross-checked to a relative tolerance
   of 10⁻⁸ over a 360-case sweep spanning `m ∈ [1.0, 5.0]` (including below,
   at, and above `m = 2`), three crack-size ranges, three stress ranges and four
   geometry factors. Neither routine calls the other.
4. **Numerical behaviour.** Convergence under interval refinement, the expected
   O(h⁴) Simpson error order, bit-for-bit determinism, life additivity over a
   split crack interval, and monotonic trends in every input.

Input validation (non-finite values, non-positive `a`, `C`, `m`, `Y`, reversed
crack limits, odd interval counts, wrong argument types) is tested throughout.

## Representative sanity result

Canonical idealized wing-skin case:

| Input | Value |
| --- | --- |
| `Y` | 1.0 |
| `σ_max` | 120 MPa |
| `σ_min` | 20 MPa |
| `Δσ` | 100 MPa |
| `R` | 0.1667 |
| `a₀` | 1.0 mm (assumed flaw) |
| `a_f` | 10.0 mm (**imposed** target) |
| Paris `C`, `m` | 1.0 × 10⁻²⁹ (SI), 3.0 |

| Result | Value |
| --- | --- |
| `ΔK(a₀)` | 5.6050 MPa·√m |
| `ΔK(a_f)` | 17.7245 MPa·√m |
| `da/dN(a₀)` | 1.7609 × 10⁻⁹ m/cycle |
| `da/dN(a_f)` | 5.5683 × 10⁻⁸ m/cycle |
| Numerical life (1000 intervals) | 776 634.4447 cycles |
| Analytical life | 776 634.4445 cycles |
| Absolute difference | 2.72 × 10⁻⁴ cycles |
| Relative difference | 3.50 × 10⁻¹⁰ |

Crack growth accelerates as the crack extends, because `ΔK ∝ √a` and
`da/dN ∝ ΔK^m`: the rate at the target is ~32× the rate at the initial flaw.

## Stress-range sensitivity

Both `σ_max` and `σ_min` are scaled, preserving `R` and varying only `Δσ`.
`a₀ = 1 mm`, `a_f = 10 mm`.

| Scale | `Δσ` [MPa] | Life [cycles] | `N / N_ref` | `scale^−m` |
| --- | --- | --- | --- | --- |
| 0.50 | 50 | 6 213 075.6 | 8.0000 | 8.0000 |
| 0.75 | 75 | 1 840 911.3 | 2.3704 | 2.3704 |
| 1.00 | 100 | 776 634.4 | 1.0000 | 1.0000 |
| 1.25 | 125 | 397 636.8 | 0.5120 | 0.5120 |
| 1.50 | 150 | 230 113.9 | 0.2963 | 0.2963 |

The observed ratios reproduce `N ∝ (Δσ)^−m` exactly for the constant-`Y` Paris
model, as required. A 50 % stress-range increase costs ~70 % of the life.

## Initial-crack sensitivity

`a_f` fixed at 10 mm, canonical stress cycle and Paris curve.

| `a₀` [mm] | `ΔK(a₀)` [MPa·√m] | Life [cycles] |
| --- | --- | --- |
| 0.5 | 3.9633 | 1 247 101.8 |
| 1.0 | 5.6050 | 776 634.4 |
| 2.0 | 7.9267 | 443 963.8 |
| 4.0 | 11.2100 | 208 730.1 |

Life falls steeply as the assumed initial flaw grows — an eightfold increase in
`a₀` costs ~83 % of the life. Most of the life is spent while the crack is small
and `ΔK` is low, which is why the assumed initial flaw size dominates the answer.

## Limitations

This model is a teaching-grade idealization. It is **not** a certification
method and its output is **not** an inspection or safe-life interval.

- **LEFM only** — small-scale yielding assumed throughout.
- **Constant geometry factor** — no finite-width correction, no `Y(a)`.
- **Constant-amplitude loading only.**
- **No crack closure** — a compressive `σ_min` contributes in full to `Δσ`.
- **No `ΔK` threshold** — every non-zero range produces some growth.
- **No fracture-toughness cutoff.**
- **No residual-strength model.**
- **No critical crack size** — the final crack length is **imposed**, not
  derived from fracture toughness.
- **No plastic-zone correction**; no elastic-plastic fracture mechanics.
- **No retardation or overload effects.**
- **No variable-amplitude spectrum**, rainflow counting, or damage summation.
- **No R-ratio dependence** in the Paris coefficients; no Forman/Walker/NASGRO.
- **No environmental effects** (corrosion, humidity).
- **No temperature effects.**
- **No probabilistic scatter** — the result is a single deterministic value, and
  real crack-growth data scatter by a factor of several.
- **Illustrative Paris data** — not measured, not a design allowable, no alloy
  claimed.
- **Single crack** — no multiple-site damage or crack interaction.

## Package structure

```
wing-skin-crack-growth/
├── pyproject.toml
├── README.md
├── src/crackgrowth/
│   ├── __init__.py           public API
│   ├── _validation.py        finite / positivity input guards
│   ├── geometry.py           ThroughCrackGeometry
│   ├── loading.py            StressCycle
│   ├── stress_intensity.py   K and ΔK
│   ├── paris.py              ParisLaw, unit conversion, da/dN
│   ├── integration.py        analytical reference + Simpson quadrature
│   ├── sensitivity.py        single-parameter sweeps
│   └── canonical.py          the canonical wing-skin case
├── tests/                    520 tests
└── examples/
    └── wing_skin_crack_growth.py
```

## Install, test, run

From a clean virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
```

Run the full test suite:

```bash
python -m pytest -q
```

Run the sanity case:

```bash
python examples/wing_skin_crack_growth.py
```

## Minimal usage

```python
from crackgrowth import (
    StressCycle, ThroughCrackGeometry, paris_law_from_mpa_basis,
    cycles_to_crack_length, analytical_cycles_to_crack_length,
)

geometry = ThroughCrackGeometry(geometry_factor=1.0)
cycle = StressCycle(sigma_max=120.0e6, sigma_min=20.0e6)   # Pa
law = paris_law_from_mpa_basis(
    c_mpa=1.0e-11, m=3.0,
    name="illustrative",
    source_note="ILLUSTRATIVE PARIS-LAW INPUT - NOT DESIGN ALLOWABLE",
)

result = cycles_to_crack_length(1.0e-3, 10.0e-3, cycle, geometry, law)
print(result.predicted_cycles)  # 776634.4447
```

## License

**No license has been chosen for this repository.** No `LICENSE` file is
present and no license metadata is declared, so default copyright applies and
no permissions are granted to others. A license may be added later.
