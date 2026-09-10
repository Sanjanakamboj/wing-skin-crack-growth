# wing-skin-crack-growth

Verified fatigue-crack-growth analysis for an idealized aircraft wing-skin crack
under cyclic loading, built on linear-elastic fracture mechanics (LEFM) and the
Paris law.

Python package: `crackgrowth` (src layout, no runtime dependencies).

## Objective

Predict the number of constant-amplitude load cycles required to grow a
through-thickness crack in a wing skin from an assumed initial flaw size to the
crack length at which the skin fractures, with every numerical result
independently verified against a closed-form reference solution.

The engineering question Milestone 2 answers:

> At what crack length does the maximum tensile stress drive the mode-I stress
> intensity up to the material's fracture toughness, and how many cycles remain
> from the initial flaw to that boundary?

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

## Scope (Milestone 2)

Added, without altering any Milestone 1 mechanics:

- mode-I fracture-toughness material representation
- a `K_max` fracture criterion (distinct from the `ΔK` growth driving force)
- critical crack size derived from `K_IC`
- residual strength as a function of crack size
- fracture utilization and a preliminary screening margin
- initial-flaw admissibility checks against the fracture boundary
- crack-growth life integrated to the toughness-derived critical size
- toughness, maximum-stress, geometry-factor and initial-flaw sensitivity sweeps
- a residual-strength screening table

Milestone 1's imposed 10 mm endpoint is retained unchanged for regression; it is
superseded, not replaced.

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

**Fracture criterion (Milestone 2):**

```
K_max(a) = Y · σ_max · sqrt(π · a)

fracture when   K_max = K_IC
```

**Two different driving forces.** Fatigue crack *growth* is driven by the
stress-intensity **range** `ΔK`, which depends on `Δσ`. *Fracture* is driven by
the **maximum** stress intensity in the cycle, `K_max`, which depends on `σ_max`
alone. Using `ΔK` as a fracture criterion is a real and expensive error, and is
guarded against by explicit tests. For constant `Y` the two are related by

```
ΔK(a) / K_max(a) = Δσ / σ_max
```

so at the critical crack size `ΔK(a_c) = (Δσ / σ_max) · K_IC` — for the
canonical cycle, exactly `5/6 · K_IC = 25 MPa·√m`, **not** `K_IC`.

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

## Fracture toughness, critical crack size and residual strength

`FractureToughness` is an immutable record storing `K_IC` on the **SI basis**
(Pa·√m), with a mandatory non-empty provenance note. Handbook values are quoted
in MPa·√m; the conversion `1 MPa·√m = 10⁶ Pa·√m` is the named helper
`mpa_sqrt_m_to_pa_sqrt_m` (with its inverse), never inferred.

**Critical crack size** — solving `K_IC = Y · σ_max · sqrt(π · a_c)`:

```
a_c = (1/π) · ( K_IC / (Y · σ_max) )²          (σ_max > 0)
```

Verified to scale exactly as `a_c ∝ K_IC²`, `a_c ∝ σ_max⁻²`, `a_c ∝ Y⁻²`.

**Residual strength** — the same equation solved for stress instead of size:

```
σ_residual(a) = K_IC / ( Y · sqrt(π · a) )
```

Verified to scale as `a^(−1/2)`, `∝ K_IC`, `∝ 1/Y`, and — the key consistency
check — `σ_residual(a_c) = σ_max` exactly, across a sweep of stresses,
toughnesses and geometry factors.

**Maximum-stress policy.** Mode-I fracture needs a tensile opening stress. If
`σ_max ≤ 0` the cycle never opens the crack, no tensile mode-I fracture boundary
exists in this model, and `critical_crack_length` returns `math.inf` with status
`NO_TENSILE_BOUNDARY`. **`abs(σ_max)` is never taken** — compression must not
manufacture a fictitious opening-fracture boundary.

**Utilization and margin:**

```
utilization = K_max / K_IC
MS_K        = K_IC / K_max − 1        (K_max > 0)
```

The screen **passes** when `K_max ≤ K_IC`; the boundary itself passes with zero
margin. For a non-opening cycle the criterion is non-governing: utilization is
`0.0` and the margin is `math.inf`. `MS_K` is a screening indicator only — it
carries no load factor, scatter allowance, material variability, thickness or
state-of-stress validity check, and is **not a certification margin of safety**.

## Initial-flaw admissibility

Before integrating, `a₀` is compared with `a_c`:

| Case | Status | Reported life |
| --- | --- | --- |
| `a₀ < a_c` | `BELOW_CRITICAL` | finite, positive |
| `a₀ = a_c` | `AT_CRITICAL` | `0.0` |
| `a₀ > a_c` | `ABOVE_CRITICAL` | `0.0` |
| `σ_max ≤ 0` | `NO_TENSILE_BOUNDARY` | `math.inf` |

The integrator is never asked to run backwards. For `a₀ > a_c` the zero means
"no crack-growth life remains — the flaw is already outside the LEFM fracture
boundary at the first application of `σ_max`", **not** "it grew to critical in
zero cycles". The structured status is the default reporting route; passing
`strict=True` raises `InitialFlawBeyondCriticalError` instead.

A zero stress range yields `math.inf` cycles even though a finite `a_c` exists:
the crack simply never grows to it.

## Life to the critical crack size

`cycles_to_critical_crack` composes verified pieces rather than restating any
physics. It computes `a_c`, classifies the initial flaw, and then calls the
**unchanged Milestone 1** integrator with `a_final = a_c`. No separate
"fracture-growth equation" exists or should exist: Milestone 2 only moves the
upper integration limit from an imposed value to a physics-derived one. The
analytical reference is likewise the Milestone 1 closed form evaluated at
`a_final = a_c`, which is asserted directly by test.

## Fracture-toughness data and provenance

> **ILLUSTRATIVE FRACTURE-TOUGHNESS INPUT — NOT DESIGN ALLOWABLE**

No measured fracture-toughness dataset was obtained or verified. The canonical
value is **illustrative**: `K_IC = 30 MPa·√m` (`3.0 × 10⁷ Pa·√m`).

It was selected **after** auditing the consequences across
`K_IC = 20, 25, 30, 35, 40, 50, 60 MPa·√m` for the canonical cycle, not by
tuning toward a preferred answer:

| `K_IC` [MPa·√m] | `a_c` [mm] | Life 1 mm → `a_c` [cycles] | `MS_K(a₀)` |
| --- | --- | --- | --- |
| 20 | 8.842 | 753 837 | 1.974 |
| 25 | 13.816 | 830 231 | 2.717 |
| **30** | **19.894** | **881 161** | **3.460** |
| 35 | 27.078 | 917 539 | 4.204 |
| 40 | 35.368 | 944 823 | 4.947 |
| 50 | 55.262 | 983 020 | 6.434 |
| 60 | 79.577 | 1 008 485 | 7.921 |

30 MPa·√m places `a_c` at 19.89 mm: comfortably above the 1 mm initial flaw,
within a plausible inspectable crack-size range, and deliberately **not** equal
to Milestone 1's arbitrary 10 mm target. It is treated as a single constant with
no thickness or state-of-stress validity check.

**No specific aerospace alloy is claimed.** The value is not measured, not
traceable to a qualification dataset, and must not be used as a design
allowable.

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

Input validation (non-finite values, non-positive `a`, `C`, `m`, `Y`, `K_IC`,
reversed crack limits, odd interval counts, wrong argument types) is tested
throughout.

Milestone 2 adds a fifth leg:

5. **Cross-criterion guards.** Explicit tests that fracture uses `K_max` and not
   `ΔK`: two cycles sharing a `Δσ` but differing in `σ_max` must give identical
   `ΔK` and different critical sizes; `K_max(a_c) = K_IC` exactly; and
   `ΔK(a_c)/K_IC = Δσ/σ_max` across several cycles and geometry factors.
   Milestone 1's canonical numbers, public API, conventions and sensitivity
   results are additionally locked by a dedicated regression module.

### A note on quadrature span

The Milestone 1 fixed uniform Simpson grid resolves `1/(da/dN)` less well as the
ratio `a_f/a₀` grows, because the integrand is steepest just above `a₀`. The
convergence **order** is unaffected — it stays O(h⁴) — but a wide-span case needs
more intervals to reach a given tolerance:

| `a_c/a₀` | rel. error at 1000 | at 2000 | at 8000 | at 20000 |
| --- | --- | --- | --- | --- |
| 9.9 | 3.4 × 10⁻¹⁰ | 2.1 × 10⁻¹¹ | 8.7 × 10⁻¹⁴ | 3.1 × 10⁻¹⁵ |
| 39.8 | 9.8 × 10⁻⁸ | 6.1 × 10⁻⁹ | 2.4 × 10⁻¹¹ | 6.2 × 10⁻¹³ |
| 139.9 | 1.4 × 10⁻⁵ | 9.1 × 10⁻⁷ | 3.6 × 10⁻⁹ | 9.3 × 10⁻¹¹ |
| 248.7 | 1.3 × 10⁻⁴ | 8.8 × 10⁻⁶ | 3.6 × 10⁻⁸ | 9.2 × 10⁻¹⁰ |

The canonical case has `a_c/a₀ ≈ 20` and is fully converged at the 1000-interval
default. Callers integrating to a distant fracture boundary should raise
`intervals`; this behaviour is covered by test rather than left implicit.

## Representative sanity result — Milestone 1 (imposed endpoint)

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

## Representative result — Milestone 2 (toughness-derived endpoint)

Same geometry, cycle and Paris curve; `K_IC = 30 MPa·√m`; `a₀ = 1 mm`.

| Fracture boundary | Value |
| --- | --- |
| Critical crack length `a_c` | 19.8944 mm |
| Admissibility of `a₀` | below critical |
| `K_max(a₀)` | 6.7260 MPa·√m |
| `K_max(a_c)` | 30.0000 MPa·√m (`= K_IC` exactly) |
| `ΔK(a_c)` | 25.0000 MPa·√m (`= 5/6 · K_IC`) |
| `σ_residual(a₀)` | 535.24 MPa |
| `σ_residual(a_c)` | 120.00 MPa (`= σ_max`, by construction) |
| Initial utilization | 0.2242 |
| Initial margin `MS_K` | 3.4603 (screening only) |

| Life to fracture | Value |
| --- | --- |
| Numerical (1000 intervals) | 881 160.7850 cycles |
| Analytical | 881 160.7798 cycles |
| Relative difference | 5.98 × 10⁻⁹ |
| `da/dN(a₀)` | 1.7609 × 10⁻⁹ m/cycle |
| `da/dN(a_c)` | 1.5625 × 10⁻⁷ m/cycle |
| Rate ratio | 88.7× |

**Comparison with Milestone 1.** The imposed 10 mm target gave 776 634 cycles;
the toughness-derived boundary at 19.894 mm gives 881 161 cycles — 104 526 more.
The M1 target happened to be conservative here, but only by accident: it was an
imposed number, and nothing guaranteed it fell below `a_c`.

## Residual-strength screen

`σ_max = 120 MPa`, `Y = 1.0`, `K_IC = 30 MPa·√m`.

| `a` [mm] | `K_max` [MPa·√m] | `σ_residual` [MPa] | Utilization | `MS_K` | Result |
| --- | --- | --- | --- | --- | --- |
| 0.5 | 4.7560 | 756.9 | 0.1585 | 5.3078 | PASS |
| 1.0 | 6.7260 | 535.2 | 0.2242 | 3.4603 | PASS |
| 2.0 | 9.5120 | 378.5 | 0.3171 | 2.1539 | PASS |
| 5.0 | 15.0398 | 239.4 | 0.5013 | 0.9947 | PASS |
| 10.0 | 21.2694 | 169.3 | 0.7090 | 0.4105 | PASS |
| 20.0 | 30.0795 | 119.7 | 1.0027 | −0.0026 | FAIL |

The screen crosses from PASS to FAIL at `a_c = 19.894 mm`, as it must.

## Fracture-toughness sensitivity

See the audit table under [Fracture-toughness data and
provenance](#fracture-toughness-data-and-provenance). Higher toughness gives a
larger critical crack and a longer life, as expected — but with **strongly
diminishing returns**: doubling `K_IC` from 25 to 50 MPa·√m enlarges `a_c`
fourfold (13.8 → 55.3 mm) yet adds only ~18 % life, because most cycles are
spent while the crack is small and `ΔK` is low. Toughness moves the endpoint
only; it leaves `ΔK` and `da/dN` at any given crack length untouched.

## Maximum-stress sensitivity at fixed `Δσ`

`Δσ` held at 100 MPa, `σ_min = σ_max − 100 MPa`, so `R` varies and the Paris
driving force does not. `a₀ = 1 mm`.

| `σ_max` [MPa] | `σ_min` [MPa] | `R` | `da/dN(a₀)` [m/cycle] | `a_c` [mm] | Life [cycles] |
| --- | --- | --- | --- | --- | --- |
| 100 | 0 | 0.0000 | 1.7609 × 10⁻⁹ | 28.648 | 923 602 |
| 120 | 20 | 0.1667 | 1.7609 × 10⁻⁹ | 19.894 | 881 161 |
| 140 | 40 | 0.2857 | 1.7609 × 10⁻⁹ | 14.616 | 838 719 |
| 160 | 60 | 0.3750 | 1.7609 × 10⁻⁹ | 11.191 | 796 278 |
| 180 | 80 | 0.4444 | 1.7609 × 10⁻⁹ | 8.842 | 753 837 |

This is the sharpest illustration of the `K_max` / `ΔK` split. `da/dN(a₀)` is
**identical down the whole column** — `ΔK` depends only on `Δσ` — so the entire
life reduction comes from `K_max` rising and pulling the fracture boundary
inward as `a_c ∝ σ_max⁻²`. Nothing about the growth rate changed.

## Stress-range sensitivity at fixed `σ_min`

Deliberately distinct from the sweep above: here `σ_min` is held at 20 MPa, so
`σ_max` moves with the range and **both** driving forces change at once.

| `Δσ` [MPa] | `σ_max` [MPa] | `ΔK(a₀)` [MPa·√m] | `da/dN(a₀)` [m/cycle] | `a_c` [mm] | Life [cycles] |
| --- | --- | --- | --- | --- | --- |
| 50 | 70 | 2.8025 | 2.2011 × 10⁻¹⁰ | 58.465 | 7 898 116 |
| 75 | 95 | 4.2037 | 7.4286 × 10⁻¹⁰ | 31.743 | 2 214 430 |
| 100 | 120 | 5.6050 | 1.7609 × 10⁻⁹ | 19.894 | 881 161 |
| 125 | 145 | 7.0062 | 3.4392 × 10⁻⁹ | 13.626 | 423 992 |
| 150 | 170 | 8.4075 | 5.9429 × 10⁻⁹ | 9.913 | 229 647 |

The life falls far faster here than in either single-effect sweep, because a
faster-growing crack is also chasing a nearer boundary.

## Geometry-factor sensitivity

`Y` hurts twice: the growth rate rises as `Y^m` while `a_c` falls as `Y⁻²`.

| `Y` | `da/dN(a₀)` [m/cycle] | `a_c` [mm] | Life [cycles] |
| --- | --- | --- | --- |
| 0.8 | 9.0156 × 10⁻¹⁰ | 31.085 | 1 820 489 |
| 0.9 | 1.2837 × 10⁻⁹ | 24.561 | 1 243 657 |
| 1.0 | 1.7609 × 10⁻⁹ | 19.894 | 881 161 |
| 1.1 | 2.3437 × 10⁻⁹ | 16.442 | 642 897 |
| 1.2 | 3.0428 × 10⁻⁹ | 13.816 | 480 458 |

A 50 % increase in `Y` (0.8 → 1.2) costs ~74 % of the life — more than either
effect alone would produce, which is verified by test.

## Initial-flaw sensitivity to fracture

`a_c` fixed at 19.894 mm.

| `a₀` [mm] | `a₀/a_c` | Utilization | Life to `a_c` [cycles] |
| --- | --- | --- | --- |
| 0.25 | 0.0126 | 0.1121 | 2 016 969 |
| 0.50 | 0.0251 | 0.1585 | 1 351 628 |
| 1.00 | 0.0503 | 0.2242 | 881 161 |
| 2.00 | 0.1005 | 0.3171 | 548 490 |
| 4.00 | 0.2011 | 0.4484 | 313 256 |

A 16× larger initial flaw costs ~84 % of the life. The assumed initial flaw size
dominates the answer far more strongly than the toughness does — a 3× range in
`K_IC` moves the life by ~34 %, while a 16× range in `a₀` moves it by ~84 %.
Utilization follows `sqrt(a₀/a_c)` exactly for constant `Y`.

## Milestone 1 stress-range sensitivity (imposed endpoint)

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

## Milestone 1 initial-crack sensitivity (imposed endpoint)

`a_f` fixed at the imposed 10 mm, canonical stress cycle and Paris curve.

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

> **The toughness-derived critical crack size is an LEFM screening boundary, not
> a certified residual-strength allowable.**

- **LEFM only** — small-scale yielding assumed throughout.
- **Constant geometry factor** — no finite-width correction, no `Y(a)`.
- **Constant-amplitude loading only.**
- **No crack closure** — a compressive `σ_min` contributes in full to `Δσ`.
- **No `ΔK` threshold** — every non-zero range produces some growth.
- **No fracture-toughness cutoff.**
- **No residual-strength model.**
- **Constant `K_IC`** — a single value, with no thickness or state-of-stress
  validity check and no plane-stress/plane-strain transition model. Real
  toughness depends strongly on section thickness.
- **`K_IC` is illustrative** unless genuinely sourced — not measured, not a
  design allowable, no alloy claimed.
- **No plastic-zone check or correction**; no elastic-plastic fracture
  mechanics. A `K_IC`-based critical crack size becomes questionable if
  small-scale yielding is violated — that is, if the plastic zone is not small
  relative to the crack length, the remaining ligament and the thickness. This
  model performs no such check, so it will report a critical size even where
  LEFM does not apply.
- **No net-section collapse check** — for large cracks or low toughness the net
  section may yield before the reported residual strength is reached, in which
  case collapse, not fracture, governs.
- **No ligament correction** and no finite-width correction, so the residual
  strength does not fall off as the crack approaches a boundary.
- **No proof load** and no load-factor, scatter or material-variability
  allowance in the reported margin.
- **No retardation or overload effects.**
- **No variable-amplitude spectrum**, rainflow counting, or damage summation.
- **No R-ratio dependence** in the Paris coefficients; no Forman/Walker/NASGRO.
  The maximum-stress sweep therefore changes `R` without changing `da/dN`, which
  is a property of this model, not of real material behaviour.
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

Run the Milestone 1 sanity case (imposed 10 mm endpoint):

```bash
python examples/wing_skin_crack_growth.py
```

Run the Milestone 2 study (toughness-derived critical crack size):

```bash
python examples/critical_crack_growth.py
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

Growing to the fracture boundary instead of an imposed target:

```python
from crackgrowth import (
    FractureToughness, mpa_sqrt_m_to_pa_sqrt_m, cycles_to_critical_crack,
)

toughness = FractureToughness(
    name="illustrative",
    k_ic=mpa_sqrt_m_to_pa_sqrt_m(30.0),                       # MPa√m -> Pa√m
    source_note="ILLUSTRATIVE FRACTURE-TOUGHNESS INPUT - NOT DESIGN ALLOWABLE",
)

result = cycles_to_critical_crack(1.0e-3, cycle, geometry, law, toughness)
print(result.critical_crack_length)   # 0.019894367886486918 m
print(result.admissibility)           # FlawAdmissibility.BELOW_CRITICAL
print(result.predicted_cycles)        # 881160.7850
```

## License

**No license has been chosen for this repository.** No `LICENSE` file is
present and no license metadata is declared, so default copyright applies and
no permissions are granted to others. A license may be added later.
