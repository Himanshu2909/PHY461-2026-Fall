# Uranyl Fluorescence Spectroscopy — Full Calculation & Error-Propagation Derivation

This document works through, from first principles, every calculation performed by
`code.ipynb` and derives every uncertainty formula the code uses — with the actual
numbers from your run substituted in at each step. All numbers below (including the
local-parabola fit coefficients and the weighted-regression covariance matrix) have
been confirmed against the console trace of an actual execution of the instrumented
notebook on your real data. It is organized in the same order as the pipeline:
**Hg calibration → drum-to-wavelength conversion → peak extraction → wavenumber
conversion → vibronic regression → force constant → error budget.**

**Methodology note (adopted uncertainty model):** wavelength uncertainty is
propagated using only the analytic drum/peak-fit term,
`u_lambda_D = |dλ/dD| · u_D_total`. No separate calibration-curve uncertainty term
is added, because the actual D→λ conversion is a PCHIP interpolant, which by
construction passes exactly through every Hg calibration point — so the calibration
curve is treated as exact, and all propagated uncertainty traces back to the
finite precision of reading the drum and fitting each peak's center.

---

## 1. Overview of the measurement chain

```
drum reading D (cm)  --[Hg calibration: D → λ, via PCHIP]-->  wavelength λ (nm)  --[λ → ν̄]-->  wavenumber (cm⁻¹)
```

Two independent scans were taken:

1. **Hg calibration scan** (300 V) — a mercury lamp with known emission lines is
   scanned to build the map `D → λ`.
2. **Uranyl fluorescence scan** (600 V) — the unknown emission spectrum, read only
   in drum units, is converted to wavelength using the map built in step 1.

Every downstream uncertainty (in wavelength, wavenumber, band spacing, ω_e, force
constant) traces back to **one root cause**: the finite resolution of the drum
reading and the finite precision of the local peak-center fit, mapped through the
local slope of the calibration curve.

---

## 2. Mercury calibration

### 2.1 Reference lines used

Eight Hg lines were resolved and paired with drum readings:

| Label | D (cm) | λ (nm) |
|---|---|---|
| Hg 404.7 | 12.87 | 404.7 |
| Hg 407.8 | 13.19 | 407.8 |
| Hg 435.8 | 13.51 | 435.8 |
| Hg 491.6 | 13.60 | 491.6 |
| Hg 546.1 | 14.25 | 546.1 |
| Hg 577.0 | 15.38 | 577.0 |
| Hg 579.1 | 16.29 | 579.1 |
| Hg 623.4 | 17.29 | 623.4 |

D(λ) is monotonic over this range (checked in code by requiring all `diff(L)` to have
one sign), which is required before either a polynomial or a PCHIP map can be trusted.

### 2.2 PCHIP is the actual D → λ map (polynomial is diagnostic only)

The **actual** `D → λ` conversion applied to the fluorescence scan is a **PCHIP**
(piecewise cubic Hermite, shape-preserving) interpolant through the 8 points above,
with `extrapolate=False`. Two reasons this matters:

- A global polynomial can wiggle non-physically between calibration points (visible
  in `Hg_Calibration_fit.png`, where the dashed PCHIP curve tracks the points tightly
  while the solid degree-2 curve visibly cuts corners near D ≈ 13.5–14 cm).
- PCHIP is monotone between knots by construction, so it can't invent a fake local
  extremum in λ(D) the way a polynomial can.
- **PCHIP passes exactly through every calibration point** (zero residual at the
  knots) — this is the property that justifies treating the calibration curve as
  exact in Section 3, rather than assigning it a separate uncertainty term.

`extrapolate=False` is important for error control too: any fluorescence point with
D outside [12.87, 17.29] cm is *dropped*, rather than being silently given an
extrapolated (and unbounded-error) wavelength. Your run reports this explicitly:

```
25 fluorescence points lie outside the Hg calibration range
and will not be assigned wavelengths.
```

A degree-2/degree-3 polynomial fit and its leave-one-out cross-validation (LOOCV)
RMSE are still computed in the notebook (42.50 nm for degree 2, 61.79 nm for degree
3), but purely as a **diagnostic** of how much curvature/scatter is present in the
calibration data — they are **not** combined into the wavelength uncertainty below,
since they describe a different function (a non-interpolating least-squares fit)
than the one actually used for the conversion.

---

## 3. Propagating drum uncertainty into wavelength uncertainty

### 3.1 Instrumental (least-count) uncertainties

Both readout devices are analog with a mechanical least count (LC); the standard
convention used here is **half the least count**, treated as a symmetric bound:

```
Drum:         LC = 0.010 cm        → u_D    = LC/2 = 0.00500 cm
Galvanometer: LC = 0.100 (×100 nA) → u_galv = LC/2 = 0.05000 (×100 nA)
```

### 3.2 First-order (delta-method) propagation through the calibration curve

For a smooth, locally linear map λ = f(D), a small uncertainty u_D in the input
propagates to the output as

```
u_lambda_D = |dλ/dD| · u_D
```

This is the first term of a Taylor expansion of f around the measured D — valid as
long as u_D is small compared to the curvature scale of f, which is true here (u_D =
0.005 cm is small compared to the cm-scale range of the calibration).

`dλ/dD` is obtained analytically from the PCHIP interpolant's own derivative
(`pchip.derivative()`), evaluated **at each specific point**, not from a single global
slope — this matters because the calibration is strongly nonlinear (the local slope
ranges from ~42 to ~168 nm/cm across the three fluorescence bands, see Section 4.3).
Because the calibration curve itself is treated as exact (Section 2.2), this
`u_lambda_D` term **is** the total wavelength uncertainty — nothing further is added
in quadrature.

---

## 4. Fluorescence peak (band) extraction

### 4.1 Why fit a local peak, not just take the max data point?

The raw scan is sparsely and irregularly sampled in drum-reading space, so the
single highest raw data point near a peak is not, in general, the true peak location
— it's biased by wherever the drum happened to land during the scan. The code instead
fits a local parabola to the points around each candidate peak and solves for the
vertex analytically, which recovers a sub-least-count estimate of the true peak
center and gives a **statistical uncertainty** on that estimate (rather than none at
all, which is what "just take the max" would give you).

### 4.2 The local quadratic vertex fit — derivation

Around a candidate center `D_guess`, take all points within `±0.025 cm`, re-center
them at their mean `x₀` (`u = D − x₀`, purely for numerical conditioning of the fit),
and fit

```
I(u) = a u² + b u + c
```

by nonlinear least squares (`scipy.optimize.curve_fit`), which returns the best-fit
parameters `(a, b, c)` **and their covariance matrix** `pcov`. A concave-down
parabola (`a < 0`, required and checked) has its maximum at

```
u_peak = −b / (2a)         →        D_peak = x₀ + u_peak
```

**Propagating the fit uncertainty to the vertex position (delta method):** treat
`D_peak` as a function of the fitted parameters, g(a, b, c) = −b/(2a) (the x₀
shift doesn't add variance since x₀ is fixed, not fitted). The gradient is

```
∂g/∂a = b/(2a²)      ∂g/∂b = −1/(2a)      ∂g/∂c = 0
```

and the variance of the vertex estimate follows from the standard linearized
(delta-method) propagation formula for a function of correlated fitted parameters:

```
Var(D_peak) = ∇g · Cov(a,b,c) · ∇gᵀ           (u_D_fit = sqrt of this)
```

This is exactly `grad @ pcov @ grad` in the code — a full quadratic form, not just
adding independent-variance terms, because `a`, `b`, `c` from the same least-squares
fit are correlated with each other.

### 4.2b Actual fitted parabolas (from the instrumented notebook run)

Running the instrumented code against the real data gives the following local fits
(5 points each, window ±0.025 cm, re-centered at `x0` = the guess drum position):

| v | D_guess (cm) | fitted (a, b, c) | Var(D_peak) (cm²) | u_D_fit (cm) |
|---|---|---|---|---|
| 0 | 13.64 | (−5428.571, −92.000, 98.346) | 3.913×10⁻⁶ | **0.001978** |
| 1 | 13.94 | (−8071.429, −83.000, 135.914) | 1.020×10⁻⁷ | **0.000319** |
| 2 | 14.23 | (−1428.571, −6.000, 55.626) | 7.094×10⁻⁷ | **0.000842** |

Because `b` (the linear term at the re-centered `x0`) is small relative to `a` at
each of these three peaks, essentially all of `Var(D_peak)` comes from the
`∂D_peak/∂b = −1/(2a)` term — i.e. the peak-position uncertainty is dominated by how
well the *curvature-normalized slope* `b` is pinned down by the 5 local points, not
by how well the curvature `a` itself is known. This is why the widest, shallowest
local parabola (v = 0, `|a| = 5429`) has the *largest* `u_D_fit`, while the
narrowest, steepest one (v = 1, `|a| = 8071`) has the smallest.

### 4.3 Combining the fit uncertainty with the reading uncertainty, and the analytic |dλ/dD| table

The fitted peak center has two independent uncertainty contributions: the statistical
uncertainty of the local parabola fit itself (`u_D_fit`, from 4.2) and the same
half-least-count reading uncertainty as any other drum reading (`U_DRUM = 0.005 cm`,
Section 3.1). Quadrature sum:

```
u_D_total = sqrt( U_DRUM² + u_D_fit² )
```

This total drum uncertainty is then propagated into wavelength using the **local**
PCHIP slope at the fitted peak position (Section 3.2), with no further terms added:

**Structured results table (as printed by the instrumented notebook — copy directly into the report):**

| v | D_peak (cm) | u_D_total (cm) | dλ/dD (nm/cm) | λ (nm) | u_lambda_D (nm) | ν̄ (cm⁻¹) | u_ν̄ (cm⁻¹) |
|---|---|---|---|---|---|---|---|
| 0 | 13.63153 | 0.00538 | 167.9969 | 497.128 | **0.9038** | 20115.53 | **36.57** |
| 1 | 13.93486 | 0.00501 | 67.1753 | 530.868 | **0.3366** | 18837.08 | **11.94** |
| 2 | 14.22790 | 0.00507 | 42.2102 | 545.157 | **0.2140** | 18343.35 | **7.20** |

Worked example (v = 0):

```
u_lambda_D = |dλ/dD| × u_D_total = 167.9969 × 0.00538 = 0.9038 nm
```

Compare this to the previous (now-abandoned) methodology, where a constant 42.50 nm
calibration term was added in quadrature and swamped this analytic term almost
completely (giving u_λ ≈ 42.51 nm for every band regardless of position). Removing
that term restores the analytic, position-dependent structure of the uncertainty:
band v = 0, at the steepest part of the calibration curve (dλ/dD ≈ 168 nm/cm), now
has by far the *largest* wavelength uncertainty of the three bands, while band v = 2,
at the shallowest local slope (dλ/dD ≈ 42 nm/cm), has the smallest — a sensible,
physically interpretable pattern that the calibration-dominated numbers previously
hid entirely.

---

## 5. Wavelength → wavenumber conversion and its uncertainty

```
ν̄ (cm⁻¹) = 1×10⁷ / λ (nm)
```

Propagating uncertainty through this nonlinear (inverse) relation via the delta
method:

```
dν̄/dλ = −1×10⁷ / λ²      →      u_ν̄ = (1×10⁷ / λ²) · u_lambda_D
```

Worked for each band:

```
v=0: λ=497.128 nm → ν̄ = 20115.53 cm⁻¹,  u_ν̄ = (1e7/497.128²) × 0.9038 = 36.57 cm⁻¹
v=1: λ=530.868 nm → ν̄ = 18837.08 cm⁻¹,  u_ν̄ = (1e7/530.868²) × 0.3366 = 11.94 cm⁻¹
v=2: λ=545.157 nm → ν̄ = 18343.35 cm⁻¹,  u_ν̄ = (1e7/545.157²) × 0.2140 =  7.20 cm⁻¹
```

---

## 6. Final band table (summary)

| v | D_peak (cm) | dλ/dD (nm/cm) | λ (nm) | ν̄ (cm⁻¹) |
|---|---|---|---|---|
| 0 | 13.63153 ± 0.00538 | 167.997 | 497.128 ± 0.904 | 20115.53 ± 36.57 |
| 1 | 13.93486 ± 0.00501 | 67.175 | 530.868 ± 0.337 | 18837.08 ± 11.94 |
| 2 | 14.22790 ± 0.00507 | 42.210 | 545.157 ± 0.214 | 18343.35 ± 7.20 |

Bands are sorted by increasing λ (= decreasing ν̄) and assigned v = 0, 1, 2 in that
order — the shortest-wavelength, highest-energy band is the ground-state vibrational
level v = 0, per the harmonic-progression model `ν̄_v = ν̄₀₀ − ω_e v`.

---

## 7. Adjacent vibrational spacings (independent consistency check)

Simple differences between neighboring bands, each with quadrature-combined
uncertainty (since the two ν̄ values being subtracted are statistically independent
measurements):

```
Δν̄(0→1) = ν̄₀ − ν̄₁ = 20115.53 − 18837.08 = 1278.45 cm⁻¹
u_Δ(0→1) = sqrt(36.57² + 11.94²)         = 38.47 cm⁻¹

Δν̄(1→2) = ν̄₁ − ν̄₂ = 18837.08 − 18343.35 = 493.73 cm⁻¹
u_Δ(1→2) = sqrt(11.94² + 7.20²)          = 13.94 cm⁻¹

Mean spacing        = (1278.45 + 493.73)/2                 = 886.09 cm⁻¹
SEM (spread-based)  = stdev([1278.45,493.73], ddof=1)/√2   = 392.36 cm⁻¹
```

The mean spacing and SEM are **unchanged** from before (they only depend on the
central ν̄ values, not on how the uncertainty is modeled), but the formally
propagated spacing uncertainties (38.47 and 13.94 cm⁻¹) have shrunk dramatically —
by roughly a factor of 60–150 — from the previous calibration-dominated values
(~2100–2300 cm⁻¹). Now, for the first time, the *empirical* SEM (392.36 cm⁻¹) is
**larger** than the *formally propagated* uncertainties, rather than the reverse.
This is an important flag, taken up in Section 8.

---

## 8. Weighted linear regression of the vibronic progression

### 8.1 Model

```
ν̄(v) = ν̄₀₀ − ω_e·v
```

### 8.2 Weights (now dominated by band v = 2, not roughly balanced)

Each point is weighted by `wᵢ = 1/σᵢ²`:

```
w₀ = 1/36.57² = 7.477×10⁻⁴
w₁ = 1/11.94² = 7.012×10⁻³
w₂ = 1/7.20²  = 1.929×10⁻²
```

This is a **qualitatively different weighting** from the previous (calibration-
dominated) analysis, where the three weights were all within a factor of ~1.5 of
each other. Now `w₂ / w₀ ≈ 25.8` — band v = 2 dominates the fit almost completely,
because its wavelength uncertainty is smallest (shallowest local calibration slope).

### 8.3 Fit result and covariance

```
ν̄₀₀  = 19579.98 cm⁻¹
ω_e   =   628.70 cm⁻¹

Covariance matrix (cm⁻¹)², order [ν̄₀₀, ω_e]:
            ν̄₀₀            ω_e
ν̄₀₀    424.695         230.041
ω_e    230.041         136.488

u_ν̄₀₀ = sqrt(424.695) = 20.61 cm⁻¹
u_ω_e = sqrt(136.488) = 11.68 cm⁻¹

Correlation coefficient:
ρ(ν̄₀₀, ω_e) = 230.041 / sqrt(424.695 × 136.488) = +0.955
```

**Both the central values and their uncertainties have changed substantially** from
the previous (calibration-uncertainty-included) analysis — not just the error bars.
Because band v = 2 now so heavily dominates the weighted fit, the line is pulled much
closer to passing through v = 1 and v = 2, at the cost of a large miss at v = 0. This
is the direct, mechanical consequence of the reweighting in Section 8.2, and it is
worth stating plainly: **the fitted ν̄₀₀ and ω_e are not robust to the choice of
uncertainty model** — they moved from (19950.9, 860.7) to (19580.0, 628.7) cm⁻¹ purely
because the relative weighting of the three points changed, with no new data.

### 8.4 Fit residuals — and a genuine statistical warning

```
v=0: measured=20115.53, predicted=19579.98+0=19579.98... 
```
more precisely, using the fitted line `ν̄(v) = 19579.98 − 628.70v`:

| v | measured ν̄ (cm⁻¹) | predicted ν̄ (cm⁻¹) | residual (cm⁻¹) | residual / u_ν̄ |
|---|---|---|---|---|
| 0 | 20115.53 | 19579.98 | **+535.55** | **+14.6 σ** |
| 1 | 18837.08 | 18951.28 | **−114.21** | **−9.6 σ** |
| 2 | 18343.35 | 18322.58 | **+20.76** | **+2.9 σ** |

```
Weighted chi-square = Σ (residual/u_ν̄)² = 314.2   (1 degree of freedom)
```

A reduced chi-square of ~314 (anything ≫ 1 is a bad fit) means the harmonic-line
model **does not fit within the stated uncertainties** once the calibration term is
removed — band v = 0 in particular sits 14.6σ away from the fitted line. This is the
direct, honest cost of the simplification: previously, the ~42 nm calibration term
was large enough to swallow this tension silently (every band's uncertainty was so
dominated by that one shared 42.5 nm term that a 500+ cm⁻¹ miss looked unremarkable
next to a ~1700 cm⁻¹ error bar). With that term removed, the analytic
drum/peak-fit-only uncertainty is small enough to reveal that **something beyond
drum-reading and peak-fit precision is causing v = 0 to deviate from a straight
line** — most plausibly mild anharmonicity (a real, physically expected downward
curvature in ν̄(v) that a 2-parameter straight-line model cannot capture), an
unresolved/blended shoulder on the raw v = 0 peak, or a genuine peak-assignment
issue. This is worth flagging explicitly in the report rather than silently
reporting R² alone.

```
SS_res = 535.55² + (−114.21)² + 20.76² = 289,749  (cm⁻¹)²   [unweighted, matches R² definition below]
SS_tot = Σ(ν̄ᵢ − mean(ν̄))²             = 1,672,952  (cm⁻¹)²

R� = 1 − SS_res/SS_tot = 0.8205
```

R� alone (0.82) looks like a reasonably good fit; it is the **weighted chi-square**
(314, computed against each point's own uncertainty) that reveals the real problem.
This is a useful general lesson: R² measures fit quality relative to the *scatter in
the data*, while chi-square measures it relative to the *stated measurement
uncertainties* — with only 3 points, and once the analytic uncertainties shrank far
below the actual point-to-point scatter, these two diagnostics started disagreeing
with each other, and chi-square is the one that should be trusted here.

---

## 9. The 0–0 transition (band origin) and excitation energy

### 9.1 Where ν̄₀₀ comes from

`ν̄₀₀` is **not** the raw measured wavenumber of the v = 0 band (20115.53 cm⁻¹,
Section 6). It is the **intercept of the weighted regression fitted in Section 8.3**,
`ν̄₀₀ = 19579.98 ± 20.61 cm⁻¹`. This follows directly from the model itself: since
`ν̄(v) = ν̄₀₀ − ω_e·v`, setting v = 0 gives `ν̄(0) = ν̄₀₀` by construction — so the
fitted intercept *is* the model's best estimate of the electronic origin, pooling
information from all three bands rather than relying on the v = 0 band alone. Its
uncertainty, `u_ν̄₀₀ = sqrt(pcov[0,0]) = sqrt(424.695) = 20.61 cm⁻¹`, is likewise read
directly off the diagonal of the fit covariance matrix already derived in Section 8.3
— it is not a new calculation, just the same regression output carried forward.
(This is also why `ν̄₀₀` differs so much from the raw v = 0 value: the regression is
dominated by bands v = 1 and v = 2, Section 8.2, which pull the fitted intercept well
below the measured v = 0 point — the same tension flagged in Section 8.4.)

### 9.2 Converting cm⁻¹ to eV

The photon energy corresponding to a wavenumber ν̄ is `E = hc·ν̄`. Using
`h = 6.62607015×10⁻³⁴ J·s`, `c = 2.99792458×10¹⁰ cm/s`, and `1 eV = 1.602176634×10⁻¹⁹ J`:

```
hc = 6.62607015×10⁻³⁴ J·s × 2.99792458×10¹⁰ cm/s = 1.986445857×10⁻²³ J·cm

hc (in eV·cm) = 1.986445857×10⁻²³ J·cm / 1.602176634×10⁻¹⁹ J/eV
              = 1.239841984×10⁻⁴ eV·cm
```

This constant, `1.239841984×10⁻⁴ eV·cm`, is exactly the conversion factor used below
— it is a fixed physical constant, not something fit to the data, so it carries no
uncertainty of its own and simply scales both `ν̄₀₀` and `u_ν̄₀₀` linearly:

```
ν̄₀₀ = 19579.98 ± 20.61 cm⁻¹
E₀₀ = 1.239841984×10⁻⁴ × 19579.98 = 2.427608 eV
u_E₀₀ = 1.239841984×10⁻⁴ × 20.61 = 0.002555 eV
```

The relative uncertainty on E₀₀ has shrunk to 0.11% (from 7.7% previously) — but
given the chi-square warning in Section 8.4, this small formal uncertainty should be
read as "precise under the stated drum/peak-fit-only error model," not necessarily
"accurate" — the v = 0 band's 14.6σ residual suggests the true uncertainty on ν̄₀₀ is
probably somewhat larger than 20.61 cm⁻¹ once whatever is causing that tension is
accounted for.

---

## 10. Symmetric-stretch force constant

### 10.1 Physical model

```
f = m_O · (2π c ω_e)²
```

with `m_O = 15.999 amu × 1.66053906660×10⁻²⁴ g/amu` and `c = 2.99792458×10¹⁰ cm/s`.

### 10.2 Uncertainty propagation

Since `f ∝ ω_e²`:

```
u_f/f = 2 · |u_ω_e/ω_e|
```

### 10.3 Numbers

```
m_O = 2.65670×10⁻²³ g
2π c ω_e = 2π × 2.99792458×10¹⁰ × 628.70 = 1.18425×10¹⁴ s⁻¹

f = m_O × (2π c ω_e)² = 3.72588×10⁵ dyn/cm

u_ω_e/ω_e = 11.68/628.70 = 0.01858  (1.86% relative)
u_f/f = 2 × 0.01858 = 0.03716  (3.72% relative)
u_f = 0.03716 × 3.72588×10⁵ = 1.385×10⁴ dyn/cm

f = (3.72588×10⁵ ± 1.385×10⁴) dyn/cm = (372.6 ± 13.8) N/m
```

Both the central value and the relative precision have changed substantially: f
dropped from ~698 N/m to ~373 N/m (because ω_e itself dropped, Section 8.3), and the
relative uncertainty tightened from ~258% to ~3.7%. The force constant is now
formally very well-determined — but this precision is a direct consequence of
adopting the analytic-only uncertainty model, and it inherits the same caveat as
Section 9: the underlying ω_e comes from a fit with a large chi-square, so 3.7%
almost certainly understates the *true* uncertainty on f. A more defensible
statement is that the force constant lies somewhere in the range implied by the
harmonic model given the visible tension in the fit — order 350–400 N/m — rather
than a tight ±3.7% band.

---

## 11. Complete error budget (hierarchy of contributions)

| Source | Value | Where it enters |
|---|---|---|
| Drum reading (half least count) | ±0.00500 cm | every D measurement |
| Galvanometer reading (half least count) | ±0.05000 (×100 nA) | intensity axis only |
| Local-parabola peak-fit uncertainty | 0.00032–0.00198 cm (band-dependent) | added in quadrature to drum reading |
| **u_lambda_D = \|dλ/dD\| × u_D_total** | **0.21–0.90 nm (band-dependent)** | **the sole wavelength-uncertainty term — see Section 4.3 table** |
| Polynomial LOOCV RMSE (42.50 nm) | diagnostic only | printed for reference; NOT added into u_lambda |
| Adjacent-spacing SEM (empirical) | 392.36 cm⁻¹ | internal-consistency cross-check |
| Weighted-regression uncertainty in ω_e | 11.68 cm⁻¹ (formal) | propagates into the 3.7% force-constant uncertainty |
| **Weighted chi-square of the regression fit** | **314 (1 dof) — poor fit** | flags that the formal uncertainties above likely understate the true uncertainty (Section 8.4) |

**Calibration is no longer treated as a source of numerical uncertainty** (Section
2.2) — the PCHIP curve is assumed exact at the knots. This removes the single
dominant term from the previous analysis, and the position-dependent
`dλ/dD`-weighted `u_lambda_D` (Section 4.3 table) now does all the propagation work.
The trade-off, made visible for the first time by this simplification, is the poor
chi-square in Section 8.4: with the calibration "cushion" gone, the data's internal
inconsistency (most likely mild anharmonicity or an unresolved feature near the v = 0
band) is now the dominant, unquantified source of doubt in ν̄₀₀, ω_e, and f — larger
in practice than any of the formally propagated numbers above.

**Unquantified systematics** (not combined numerically, because they weren't
independently measured in this run): mechanical backlash, slit-width/instrumental
resolution, PMT gain drift, lamp intensity drift, background/stray light, unresolved
Hg doublets, incorrect peak assignment, and interpolation/model inadequacy (including
possible real anharmonicity in the vibronic progression, which the chi-square result
above makes a live candidate rather than a purely hypothetical caveat).

---

## 12. Final results (headline numbers)

```
0-0 transition:   lambda_00 = 497.13 nm (v=0 band; regression intercept differs -- see below)
                  nu_00 (fit) = 19580.0 +/- 20.6 cm^-1   (formal, drum/peak-fit-only uncertainty)
                  E_00 = 2.4276 +/- 0.0026 eV            (formal; see Section 9 caveat)

Vibrational spacing:  omega_e = 628.7 +/- 11.7 cm^-1     (formal; see Section 8.4 -- fit chi^2 = 314, poor)

Force constant:   f = (372.6 +/- 13.8) N/m               (formal; true uncertainty likely larger, Section 10.3)
```

Compared to the previous (calibration-uncertainty-included) analysis, removing the
42.5 nm calibration term shrinks every formal error bar by one to two orders of
magnitude, but it also exposes a real weighted chi-square problem in the 3-point
harmonic-line fit that the old, calibration-dominated error bars had been masking.
The honest summary for the report: this simplified uncertainty model gives tight,
band-position-dependent wavelength uncertainties directly traceable to drum/peak-fit
precision (Section 4.3's table), but the resulting weighted regression does not fit
within those uncertainties (Section 8.4) — which is itself a meaningful experimental
finding worth reporting, not a flaw to hide.
