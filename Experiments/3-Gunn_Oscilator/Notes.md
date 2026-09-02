# Experiment 20: Study of Gunn Oscillator — Complete Theory, Setup Rationale, and Data Analysis Guide

---

## 1. Introduction and Motivation

This experiment has two intertwined goals: understand a genuinely quantum-mechanical/solid-state transport phenomenon (the Gunn effect) that produces microwave oscillations from a **DC bias alone** (no external resonant circuit needed to start the oscillation, unlike a Gunn-diode-in-a-cavity oscillator which needs the cavity only to *select* frequency, not to *create* negative resistance), and to get hands-on fluency with the standard microwave bench: isolators, PIN modulators, frequency meters, slotted lines, and square-law detectors. The five sub-experiments (a)–(e) are not independent — they build a chain: first you characterize the Gunn diode itself (a), then you use it as a source (b), then you characterize the amplitude modulator that rides on it (c), then you characterize the detector you're using to make all these measurements (d), and finally you calibrate the one "uncalibrated" component in the chain, the micrometer attenuator (e). Reading the manual out of order can be confusing for exactly this reason — each part supplies a tool needed by the next.

---

## 2. Theory of the Gunn Effect — In Depth

### 2.1 The band structure origin: the two-valley (multi-valley) model

Gallium Arsenide (GaAs) and other III-V semiconductors (InP, CdTe, GaN under some conditions) have a **conduction band with more than one local minimum ("valley") in k-space**, at different energies:

- **Γ-valley (central valley)**: This is the lowest-energy conduction band minimum, located at k = 0. Electrons here have a *small effective mass* (m*ℓ ≈ 0.067 m₀ in GaAs), which by the relation μ = eτ/m* gives them **high mobility** (μℓ ≈ 8000 cm²/V·s at room temperature).
- **L-valley (satellite/upper valley)**: Located away from k = 0 (along the ⟨111⟩ directions in GaAs), roughly **0.29–0.36 eV above** the Γ-valley. Electrons here have a much *larger effective mass* (m*h ≈ 0.55 m₀), and hence **low mobility** (μh ≈ 100–200 cm²/V·s).

This is exactly the "light mass band" and "heavy mass band" the manual refers to. At thermal equilibrium (low field), essentially all conduction electrons sit in the favorable low-energy, high-mobility Γ-valley, so the material behaves as an ordinary Ohmic conductor with a single effective mobility.

### 2.2 Why a strong electric field transfers electrons between valleys

When you apply a strong electric field, conduction electrons gain **drift energy** from the field between scattering events. Once their kinetic energy becomes comparable to the 0.3 eV Γ–L energy gap, a real-space-preserving but momentum-scattering process (electron–phonon inter-valley scattering) can kick an electron from the Γ-valley into the L-valley. This is called **inter-valley transfer** or the **Ridley–Watkins–Hilsum (RWH) mechanism**. It requires:

1. A large enough energy separation between valleys that the *lower* valley isn't thermally depopulated at equilibrium (else you'd get NDM even at zero field), but small enough (~0.3 eV, a few kT at typical operating temperatures once field-heated) that realistic fields (few kV/cm) can drive the transfer.
2. A large density-of-states in the upper valley (satisfied automatically since heavier effective mass → higher DOS, m*^(3/2)) so that once electrons start transferring, there's "room" for them there, reinforcing the transfer.

### 2.3 Deriving the negative differential mobility (NDM) — filling in the manual's algebra

The manual states the total current density as the sum of both populations:

**J = σ_eff·E = (nℓμℓ + nhμh)·eE**, with nℓ + nh = N (total electron density, constant).

Differentiating J with respect to E (this is done because the *differential* conductivity dJ/dE, not the *DC* conductivity J/E, is what determines whether the material is stable or breaks into oscillation — see §2.4):

**dJ/dE = (nℓμℓ + nhμh)e + [μℓ(dnℓ/dE) + μh(dnh/dE)]eE**

Since nℓ + nh = N = const, dnℓ/dE = −dnh/dE, so substituting:

**dJ/dE = (nℓμℓ + nhμh)e + (μh − μℓ)eE(dnh/dE)**

For **dJ/dE < 0** (the NDM condition — the differential mobility becomes negative, meaning increasing E actually *decreases* current) we need:

**(nℓμℓ + nhμh) < (μℓ − μh)E(dnh/dE)**

Physically, this inequality is satisfied when:
- **μℓ ≫ μh**: the mobility contrast between valleys must be large (in GaAs, roughly 40:1), so that transferring electrons from ℓ to h causes a large mobility *drop* per electron transferred.
- **dnh/dE is large and positive**: the field must strongly and rapidly repopulate electrons from Γ to L as E increases — i.e., you need strong inter-valley scattering right around the operating field.
- **High E**: you need enough field to actually be in the transfer regime; at low field essentially all electrons stay in Γ and dnh/dE ≈ 0, so dJ/dE > 0 (Ohmic).

**Physical picture**: As E increases past a threshold field Eth (~3.2 kV/cm in GaAs), electrons that were fast, light, high-mobility carriers get scattered into the slow, heavy, low-mobility valley faster than the extra field can push more current through them. The *increase* in field is outpaced by the *decrease* in average mobility, so net current falls even though voltage rises. This is why NDM is fundamentally a "more voltage → your average carrier gets heavier and slower faster than it gets more numerous or pushed harder → less current" story. It is entirely a bulk transport effect — no p-n junction, no depletion region physics like an ordinary diode. This is important for Question 3/4 territory: the "Gunn diode" is a two-terminal piece of *n-type* bulk GaAs (not really a diode/junction device at all — the name is a historical misnomer), so its I-V curve is symmetric in principle and has nothing to do with rectification.

### 2.4 From NDM to spontaneous oscillation: domain formation

A material with dJ/dE < 0 over some range of E, held at a DC bias inside that range, is **electrically unstable**. Here's why, step by step:

1. In a real device, doping is never perfectly uniform. Somewhere there is a small local fluctuation — a slightly higher doping density or slightly lower cross-section — that produces a small local field spike relative to the average field.
2. Because dJ/dE < 0 in this region, that small field increase locally causes a local current *decrease* relative to the surroundings.
3. Current continuity (J must be uniform along a 1-D conductor at any instant, since charge can't pile up on the timescale of interest without a compensating field change) forces the material on either side of this fluctuation to *increase* its local field to keep the current matched — this piles up a **dipole layer**: excess negative charge accumulates on the cathode-facing side of the fluctuation, excess positive (charge depletion) on the anode-facing side.
4. This dipole layer becomes a self-reinforcing **high-field domain**: the field inside the domain is much higher than outside it (the domain is now in an even deeper NDM region, or beyond it), while the field outside relaxes to keep the total voltage (∫E dx = V_applied) fixed.
5. This domain is *not static* — it is swept along by the same drift velocity vd that all carriers have inside it, so it **travels from the cathode to the anode** at vd ≈ 10⁷ cm/s.
6. When the domain reaches the anode, it is extinguished (absorbed at the contact), the field inside the device relaxes back to roughly uniform, and — because you're still holding it at a DC bias inside the NDM region — the instability condition immediately re-triggers and a **new domain nucleates at the cathode** and the cycle repeats.

Each domain transit constitutes one cycle of a current oscillation at the external terminals (since while the domain exists, the bulk of the sample sees a *different* field than the domain, so the total current — same everywhere by continuity — is set by the *domain's* velocity-field relation, and this jumps each time a domain is born/dies). This gives the **transit-time oscillation** the manual quotes:

**ν = v/d** (Eq. 4)

with v ≈ 10⁷ cm/s the domain (~carrier drift) velocity and d the cathode-to-anode length of the active GaAs slab. This is why part of your data analysis (Question 2) is to *invert* this relation: you measure the operating threshold field electrically, know d is what you're solving for, and can then estimate the frequency this predicts, or vice versa back out d from a known/assumed operating frequency.

**This transit-time mode is called the "Gunn mode" or sometimes the "transit-time dipole-domain mode."** It's worth knowing (for completeness, though not asked) that there are other modes (LSA — Limited Space-charge Accumulation, quenched-domain, delayed-domain) obtained by combining the diode with an external resonant cavity so the cavity field, not pure transit time, controls when domains nucleate/quench — these give higher frequencies and higher efficiency, which is why "real" Gunn oscillators in radar/communications systems are diode+cavity assemblies, not bare slabs. In your experiment, since you observe amplitude modulation via the PIN diode in part (b) rather than directly demodulating and counting oscillation cycles, you're using the diode as a transit-time (or resonator-tuned, since it sits inside a waveguide mount which does provide some cavity loading) microwave *source*, not resolving its internal domain dynamics directly.

### 2.5 The Gunn diode's DC I-V curve (Fig. 1) explained

At low bias, the device is Ohmic: I rises linearly with V (nearly all electrons in the Γ-valley, constant mobility). As V approaches V₀ (threshold), inter-valley transfer becomes significant, dJ/dE starts falling, and the I-V curve bends over, reaching a peak at V₀ (**this is the threshold voltage**, corresponding to the threshold field Eth = V₀/d). Beyond V₀, in the **negative differential conductance regime**, further increasing the DC-measured average voltage (which now represents the sum of a high-field domain plus a lower-field bulk region) causes the average measured current to *decrease*, because more of the sample's length is now inside the slow high-field domain. This is precisely the regime where domains nucleate and the device oscillates — hence the manual's warning **not to sit at V₀ for long periods** (the domain formation dissipates power as heat and can thermally damage the small active region) and to always return to full CCW (minimum bias) between readings.

---

## 3. Experimental Setup — Component-by-Component Rationale

Consult Fig. 2 in the manual: **Gunn Oscillator → Isolator → PIN Modulator → Frequency Meter → Attenuator → Slotted Section (shorted) → detector → indicator**, with the Gunn assembly tuning micrometer and DVM on the source end, and the uncalibrated micrometer attenuator riding on the PIN modulator stage.

### 3.1 Gunn oscillator + tuning micrometer + DVM

The Gunn diode chip is mounted inside a **waveguide cavity** (not free space) because the waveguide provides mechanical support, a defined characteristic impedance to couple power out into the rest of the waveguide train, and — critically — a resonant structure that can be mechanically tuned. The **tuning micrometer** moves a plunger/short that changes the effective cavity length, which pulls the oscillation frequency (this is why changing the micrometer setting in part (b) changes the observed frequency — you are loading the diode's intrinsic transit-time oscillation with a cavity resonance that can pull it up or down within a limited range, which is precisely the "electronic + mechanical tuning" behavior real Gunn oscillators exploit). The **DVM** reads the DC bias voltage across the diode directly (via the BNC T-connector so you can bias and read simultaneously), letting you build the I-V curve of part (a).

### 3.2 Isolator

A ferrite-based non-reciprocal device that passes microwave power in the forward direction (Gunn → PIN modulator) with very low loss (~0.5 dB) but presents very high loss (~20–30 dB) to any power reflected back toward the source. **Why it's essential here (Question 4):** the Gunn diode's oscillation frequency and amplitude are sensitive to the impedance it sees looking into the waveguide ("frequency pulling" and "load pulling"). Every other component downstream — especially the PIN diode, whose impedance is deliberately being switched between two very different states to create modulation — presents a **time-varying, mismatched load**. Without the isolator, this time-varying mismatch would reflect power back into the Gunn diode and directly modulate/perturb its own oscillation (frequency and power), corrupting the very square-wave modulation you're trying to characterize, and would also make the Gunn source's frequency dependent on whatever attenuator/load configuration is downstream — i.e., the source would not be a stable, independent oscillator. The isolator decouples the source from the downstream network so the source's characteristics (measured in part (a)) remain valid once it's connected into the full chain in parts (b)–(e).

### 3.3 PIN diode modulator

A PIN diode is a p-i-n junction (p-doped / intrinsic / n-doped) whose **microwave-frequency resistance is controlled by its DC/low-frequency bias**, not by conventional p-n rectification (at microwave frequencies, the RF period is far shorter than the minority-carrier lifetime in the intrinsic region, so the diode cannot "rectify" the microwave signal itself — instead it behaves as a bias-dependent variable resistor):
- **Forward bias**: minority carriers are injected into and stored in the intrinsic region, giving it a low resistance (RF signal passes with little attenuation — "ON"/transmit state).
- **Reverse bias (or zero bias)**: the intrinsic region is depleted of free carriers, giving it a high resistance (RF signal is strongly attenuated/reflected — "OFF"/block state).

By driving the PIN bias with a **square wave**, the through-line resistance is switched periodically between these two states, so the transmitted microwave amplitude is chopped into a square-wave-modulated (ON/OFF, i.e., amplitude-shift-keyed) envelope — exactly the "square waves on the CRO" the manual describes. Because it modulates via a **series/shunt resistance change**, not via mixing or frequency generation, it's a passive, broadband, and fast (nanosecond-scale switching) modulator, which is why it's the standard component for amplitude-modulating a fixed-frequency microwave carrier in a lab bench setting, as opposed to modulating the Gunn bias directly (Question 3 — you *can* amplitude-modulate by varying the Gunn bias itself, since output power depends on bias, but doing so simultaneously shifts the *frequency*, called AM-to-FM conversion/chirping, which is undesirable if you want clean, frequency-stable AM; the PIN diode modulates amplitude while leaving the source's frequency undisturbed, which is the whole reason it's inserted downstream of an isolator rather than just wiggling the Gunn bias).

### 3.4 Frequency meter

This is a precision **resonant cavity wavemeter** — a tunable cavity (via a micrometer dial calibrated directly in GHz or cm) coupled weakly to the waveguide. As you tune it through the operating frequency, at resonance it absorbs a small notch of power from the line (you'd see a dip on the VSWR meter/CRO), letting you read the oscillation frequency directly off its calibrated micrometer dial. It's placed after the PIN modulator (still before the attenuator/slotted section) simply so it can sample the actual propagating signal; its own insertion loss is small so it doesn't need to be switched out for the amplitude measurements downstream.

### 3.5 (Uncalibrated) micrometer attenuator

A **precision variable RF attenuator**, typically implemented as a resistive vane that can be inserted more or less deeply into the waveguide cross-section (deeper insertion near the E-field maximum → more absorption/attenuation). Its micrometer reading is a *mechanical* depth, not directly calibrated in dB — that's precisely the reason part (e) exists: to build a mechanical-reading-to-dB calibration curve by comparing it against a **known, dB-calibrated** attenuator. It's used throughout the experiment to keep detector/VSWR meter readings on-scale (protecting the detector from excessive power and keeping you within its square-law linear regime — see §4).

### 3.6 Slotted section (terminated with a short)

A waveguide section with a longitudinal slot cut into the *center of the broad wall*, through which a small probe (with a crystal detector on it) can be moved along the guide axis to sample the E-field magnitude at each point without significantly perturbing the field. **Why the slot is in the center of the broad wall, not the narrow wall (Question 5):** for the dominant TE₁₀ waveguide mode, the surface current pattern on the broad wall has a **longitudinal (axial) component that vanishes exactly along the center line** of the broad wall. A slot cut along a line where the wall current has no component crossing the slot does not interrupt the wall currents, so it doesn't radiate/leak power or perturb the field distribution — this is a "non-radiating" slot placement. If the slot were cut in the narrow wall (or off-center in the broad wall), it would cross current lines, causing radiation leakage, field perturbation, and unreliable field-sampling as well as loss of power out of the guide. This is precisely why slotted lines are always built with the slot centered on the broad wall.

Terminating the far end in a **short circuit** (rather than a matched load) is deliberate: a short reflects essentially 100% of incident power, producing a **standing wave** with well-defined, easily located voltage maxima and minima (nodes spaced λg/2 apart, with minima being sharp nulls at a short). This standing wave is exactly what you need for two purposes: (i) determining the guide wavelength λg in part (d)/method (ii) from the node spacing, and (ii) — although not explicit here — shorted-line probing is the standard technique to precisely locate field maxima for probe calibration in general microwave benches.

### 3.7 Detector (1N23/1N23B crystal diode) + indicator (VSWR meter / micro-ammeter / CRO)

The 1N23 is a **point-contact Schottky (metal-semiconductor) diode**, chosen for microwave detection because Schottky junctions have negligible minority-carrier storage and hence very low junction capacitance and very fast response — needed to rectify signals at GHz frequencies where a conventional p-n diode's junction capacitance and carrier lifetime would make it far too slow/lossy. Rectification converts the RF envelope into a DC/low-frequency current proportional (to some power n) of the incident RF voltage, which is exactly the "response law" I = kVⁿ investigated in part (d). The **VSWR meter** is simply a tuned, calibrated audio-frequency amplifier + meter designed to read the detector's output when the source is *modulated* (it's tuned to the modulation frequency, typically 1 kHz, and reads in dB directly — hence its use in the "square modulated microwave" method (i)). The **micro-ammeter** is used instead when the source is unmodulated CW (method (ii)), since there's no AC modulation envelope to lock onto — you just read the detector's DC rectified current directly. The **CRO** (oscilloscope) is used in part (c) because you need to see the actual *time-domain waveform* (the square-wave envelope shape, Vmax and Vmin levels) to compute % modulation, which a scalar meter reading cannot give you.

---

## 4. The Square-Law Detector and Response Law — Theory Behind Part (d)

### 4.1 Why I = kVⁿ, and why n = 2 is special

A Schottky diode's DC I-V relation is exponential: I = I₀(e^(qV/ηkT) − 1). For a *small* RF signal V(t) = Vrf·cos(ωt) superposed on some (usually zero or small) DC bias, Taylor-expanding this exponential in powers of V(t)/(ηkT/q) and time-averaging over one RF cycle shows that the DC (rectified) output current is proportional to the **mean square** of the input RF voltage, i.e., proportional to the *incident microwave power* (since power ∝ V²), giving n = 2 — the "square law" regime. This is valid only for **small input RF signal levels** relative to the diode's characteristic voltage kT/q (~26 mV at 300 K); as you increase the input power, higher-order terms in the exponential expansion become non-negligible and the response deviates from a pure square law (n departs from 2, typically dropping toward n = 1, "linear detection," at high power, since a strongly driven diode behaves more like a simple rectifier/peak detector rather than a small-signal square-law device). This is exactly the "obeyed for low power levels, deviates at high power" behavior the manual mentions, and is the physical reason part (d) asks you to repeat the measurement at multiple attenuator/gain settings — you're mapping out over what power range the detector is square-law, which is essential to know because *all other parts of this experiment implicitly assume you're reading the detector in its square-law region* when you interpret meter deflections as proportional to power.

### 4.2 Method (i): log-log calibration against a calibrated attenuator

Here the key trick is: the attenuator's dB reading is *already* a log measure of power (dB = 10 log₁₀(P), so a change of X dB corresponds to a known, fixed multiplicative change in power, i.e., in V², regardless of absolute calibration constants). The VSWR meter's dB scale is *also* a log measure, but of the *detector current*, not directly of power — this is fine precisely because the VSWR meter is calibrated in dB assuming a square-law detector (dB reading ∝ log(I) and you want to confirm log(I) ∝ log(P) with unit slope). So: **plotting VSWR-meter dB reading vs. attenuator dB setting should give a straight line of slope 1 if, and only if, the detector is behaving as a perfect square-law device.** Any deviation of the slope from 1 quantifies the exponent n directly (since attenuator-dB ∝ log(V²_in on the line) and VSWR-dB ∝ n·log(V_at detector after that fixed system loss), the ratio of slopes gives n/2, or equivalently you interpret the graph's slope as n directly if you've set it up as the manual suggests — work through the log algebra explicitly when you analyze your own data, since the precise slope-to-n mapping depends on exactly which axis is "attenuator dB decreasing to the right" as instructed).

### 4.3 Method (ii): direct log-log from standing-wave field profile

Here you don't need a second calibrated attenuator at all — instead you use the **known sinusoidal shape of the standing wave itself** as your independent, absolutely-known "input signal" whose relative amplitude you can compute at every probe position from Eq. 7. You measure raw detector current I at each probe position x, plot log(I) vs. log(sin βx) (i.e., log of the known relative field amplitude), and the slope directly gives n. This is elegant because it requires no absolute calibration of anything — pure geometry (probe position x, known λg) substitutes for a calibrated power reference.

### 4.4 Determining λg (guide wavelength) — necessary first step for method (ii)

Because the guide is shorted, you get a clean standing wave; **λg is twice the spacing between successive current minima (or maxima)** as you slide the probe (since the standing wave pattern in *power/current* — which goes as sin²βx for a square-law detector — has period λg/2, i.e., two minima per guide wavelength, while the *field* E itself has period λg). Be careful here: because the detector response is sin²(βx) (proportional to power, not field, if truly square-law) rather than sin(βx), the current pattern minima will be twice as closely spaced as you might naively expect from the field formula (Eq. 6/7) alone — always locate minima (nulls), which are much sharper and easier to pinpoint precisely than the broad maxima, exactly as the manual's method instructs ("select a convenient minimum ... move towards a maximum").

Note also: **λg ≠ λ₀ (free-space wavelength)**. Waveguide propagation is dispersive; the relation is
**1/λg² = 1/λ₀² − 1/λc²**, where λc is the cutoff wavelength of the waveguide (for TE₁₀ mode in standard X-band rectangular waveguide, λc = 2a where a is the broad-wall inner dimension, ≈ 2.286 cm giving λc ≈ 4.57 cm). This is worth knowing if a question asks you to cross-check your frequency-meter-measured frequency against your λg-measured value, or vice versa, via c = f·λ₀.

---

## 5. Data Analysis — Part by Part

### (a) I–V characteristic and threshold voltage — including how to handle the fluctuating/oscillatory region

**Recognize the four (not two) regimes in your own data.** Real student data (like yours) on this experiment essentially always shows four distinct regions as you sweep bias up from zero, not just a clean "rise then fall":

1. **Low-field Ohmic region** (roughly 0.1 V up to ~2.5–2.8 V in your 5 mm run): current rises smoothly and close to linearly. Nearly all electrons are still in the light-mass Γ-valley; this is ordinary Ohmic conduction.
2. **Near-threshold quasi-saturation** (~2.9–3.3 V, where your current plateaus around 500–501 mA): the rate of rise slows sharply and current appears to "cap out." This is the region right around V₀ where inter-valley transfer is starting to compete with the increasing field — the static differential conductance dI/dV is falling toward zero here.
3. **NDM onset / sudden drop**: immediately after the plateau you see a sharp, single-valued drop (in your data, roughly 501 mA → ~315 mA over a small voltage change). This is the device crossing into net negative differential mobility and beginning to sustain a traveling high-field domain.
4. **Sustained-oscillation (fluctuating) region**: this is your "oscillatory and fluctuatory" region, and it is the most important one to understand correctly. **This is not measurement error and not a sign you did something wrong.** Once a domain is continuously nucleating, transiting, and being extinguished at the anode (§2.4), the *instantaneous* terminal current is genuinely oscillating at the transit-time frequency (~GHz). Your panel meter (DPM) cannot track a GHz oscillation — it displays some time-average of it — but three things make that time-average itself unstable and non-repeatable at fixed knob position: (i) the DPM's own integration/response time only partially averages the fast oscillation, so its displayed value drifts depending on exactly which few milliseconds it's sampling; (ii) domain nucleation is not perfectly periodic — there is timing jitter, and the device can intermittently "mode-hop" between slightly different numbers of domains in transit or slightly different domain velocities, each of which corresponds to a different average current; (iii) because the diode's dynamic negative resistance interacts with the finite output impedance of your DC power supply, small relaxation oscillations can appear in the **bias voltage itself**, which is exactly why you also see the *voltage* jittering, not just the current, and why some of your recorded points even show voltage decreasing slightly between successive rows despite the sweep nominally going up. Eventually (in your data, above roughly 4.2–5 V) the device settles into a **sustained, fully periodic oscillation state** whose time-averaged current vs. bias becomes smooth again — but at a *much flatter* slope than the naive NDM picture would suggest, because the domain itself now self-regulates the internal field. This is why your data goes fluctuating → smooth-but-slowly-varying (roughly flat, then a gentle decline from ~420 mA down to ~396 mA by 8 V) rather than continuing to plunge.

**How to plot it:**
- Do **not** draw one single connected smooth line through the entire fluctuating region — that manufactures a false sense of precision. Instead:
  - Plot regions 1–3 (Ohmic rise through the sharp drop) as a normal connected line/curve — these points are single-valued and repeatable.
  - Plot region 4 (fluctuating region) as a **scatter of individual points** (no connecting line), so the spread is visible to the reader. If your raw data has multiple current readings recorded near the same nominal voltage (as yours does — e.g., several entries clustered around V ≈ 3.2–3.7 V with currents ranging from ~320 to ~411 mA), plot all of them as separate points at their own (V, I) coordinates rather than averaging them away silently.
  - Once the data becomes smooth again at higher bias (region 5, e.g. above ~4.2–5 V in your run), switch back to a connected line/curve — this part is real, repeatable, quasi-static behavior of the sustained-oscillation state and should be treated normally, including a linear fit if it looks linear (see Error Analysis below).
  - Add a shaded band or bracket on the x-axis labelled "domain-oscillation / unstable regime" over region 4 so a reader immediately understands why the data looks different there — this is expected physics, and saying so explicitly in your report turns what looks like "bad data" into a correctly-identified physical regime.

- **Locating V₀:** use the onset of the *sharp drop* (start of region 3) as your primary criterion for V₀, since this is the actual threshold field crossing. As a cross-check, also note the bias at which sustained fluctuation *begins* (start of region 4) — physically these should very nearly coincide, but a small gap between them is itself meaningful (it reflects the short delay/hysteresis between the point where the static I-V first goes non-monotonic and the point where a fully-formed domain is actually launched and sustained). Report both bias values if they differ, and briefly note the gap rather than silently picking one.
- **Uncertainty in V₀:** since the plateau/drop is not a sharp cusp, fit a smooth curve (e.g. a local polynomial) to the 5–8 points spanning the plateau-into-drop and take the inflection/maximum from the fit rather than eyeballing a single row of the table; quote V₀ ± half your voltage step size as a floor on the uncertainty even after fitting.
- Repeat this whole analysis for your second micrometer setting (8 mm) — your 8 mm data shows the same four-regime structure (Ohmic rise to ~501 mA plateau near 3.0–3.3 V, then an explicit "Fluctuations →" region you've already flagged by hand through roughly 3.4–5.4 V, then smoothing out again above ~5.5 V). Compare V₀ between the two micrometer settings; since V₀ = Eth·d is a bulk material/geometry property, not a cavity-tuning property, the two values should agree within your uncertainty — if they don't, discuss whether it's within scatter or a real systematic effect (self-heating, slightly different effective coupling/contact conditions at different micrometer settings).

### (b) Gunn as a source of microwaves — detailed procedure and analysis

**Step-by-step procedure** (the manual's version skips several practical steps you actually need on the bench):

1. Before switching anything on, confirm both the GUNN bias knob and the PIN bias knob are full CCW (minimum), exactly as the safety note under part (a) requires.
2. Set the Gunn-assembly tuning micrometer to a fixed reference position (e.g. 5 mm, matching your part (a) run).
3. Switch on the Gunn power supply, the VSWR meter, and the CRO. Set the **modulation selector switch on the Gunn power supply to MOD** (not CW). This is essential and easy to miss: a truly continuous, unmodulated microwave carrier produces no AC signal for a VSWR meter (which is an AF amplifier tuned near ~1 kHz) to lock onto — you need the PIN diode chopping the carrier into a square wave before the VSWR meter can show you anything.
4. Gradually increase the Gunn bias while watching the DPM, using your part (a) data as a guide, until you are clearly inside the **sustained-oscillation region** (region 4/5 as identified in part (a) above) — below V₀ there is no microwave emission at all, only DC current.
5. Set the PIN bias to roughly the midpoint of its range (knob "midway between full CCW and full CW," as the manual says) — this gives close to maximum, clearly-visible modulation depth to work with.
6. On the CRO, adjust the time base and vertical gain until you can clearly see a square wave. If you see nothing: check that the slotted-line probe isn't sitting at a field null, that the attenuator isn't set too high, and that the Gunn bias is genuinely above V₀ (step 4).
7. Slide the slotted-section probe along the guide to the position of **maximum** indicator deflection (a field maximum) — you want the best signal-to-noise for the frequency reading that follows.
8. With the signal maximized, tune the **modulation-frequency knob** on the Gunn power supply (this sets the audio-rate square-wave chopping frequency fed to the PIN diode — a completely separate, much lower frequency than the ~10 GHz microwave carrier) until the VSWR-meter indicator reading is maximized. This matters because the VSWR meter's internal amplifier is tuned near a specific audio frequency (often ~1 kHz); if your PIN chopping rate doesn't match it, the meter under-reads even though the microwave signal itself is fine.
9. Note the Gunn bias voltage (DVM) you're using. Now slowly rotate the **frequency meter's** micrometer dial through its range while watching the indicator — at resonance you will see a **dip** (absorption of a small notch of power). At the point of the deepest/sharpest dip, read the operating frequency directly off the frequency meter's calibrated dial.
10. Without touching the Gunn bias, change **only the Gunn-assembly tuning micrometer** to a new position and repeat step 9 — this builds one frequency-vs-micrometer curve at fixed Gunn bias.
11. Return the Gunn-assembly micrometer to its original reference and instead change the **Gunn bias** to a second fixed value (still safely inside the oscillation region), then repeat steps 8–10 to get a second curve.

**Analysis:**
- Plot frequency (y-axis, from the frequency-meter dial) vs. Gunn-assembly-tuning-micrometer reading (x-axis), as two separate curves — one per Gunn bias value used.
- Expect a roughly monotonic but *bounded* tuning curve: the cavity only "pulls" the diode's intrinsic transit-time frequency over a limited range (tens to a few hundred MHz), small compared to the ~10 GHz carrier itself — don't expect huge frequency swings across the micrometer's travel.
- Compare the two curves: a shift in the curve (not just its range) between the two Gunn-bias values shows that bias changes the operating point (and hence which part of the domain-formation dynamics sets the frequency), on top of the purely mechanical cavity-length effect of the micrometer.

### (c) Characteristics of the PIN modulator — detailed procedure and analysis

**Step-by-step procedure:**

1. Start from a working configuration from part (b) that gave a strong, clean square wave on the CRO (same Gunn bias, Gunn-micrometer setting, and probe position).
2. Tune to a suitable, convenient frequency and set the (uncalibrated) attenuator to a fixed reference value, e.g. ~3 dB (note this value — you will need it as a fixed reference so the CRO signal stays on-scale as you vary PIN bias in the following steps).
3. Set the **CRO input coupling to DC**, not AC. This is essential: you need the *absolute* Vmax and Vmin levels relative to true zero, not just the peak-to-peak swing — AC coupling strips out the DC reference and makes the individual Vmax/Vmin readings meaningless (only their difference would survive, which is not enough to compute %m).
4. Establish the true zero level: momentarily ground the CRO input channel and mark where the trace sits on the graticule — this is your reference for "0 volts" for every reading that follows. Re-check this occasionally, since scope DC offset can drift over a session.
5. Reconnect the detector output and adjust vertical gain/time base for a large, clearly readable square wave that doesn't clip off the top or bottom of the screen.
6. Starting from PIN bias at (or near) zero, step the PIN bias up through its range in convenient, regular increments (8–12 points is enough to trace the full curve). At **each** setting: read the PIN bias value (via the DVM/CRO as connected through the BNC T-connector, per the manual), and read Vmax (top of the square wave) and Vmin (bottom) directly off the CRO screen relative to the zero line from step 4.
7. Compute %m = (Vmax − Vmin)/(Vmax + Vmin) × 100 at each PIN bias value.

**Analysis:**
- Plot %m (y-axis) vs PIN bias (x-axis). Expect %m to start near 0% (PIN diode still low-resistance/"ON," barely chopping the carrier), rise steeply through an intermediate bias range as the PIN diode's resistance climbs toward its high-resistance "OFF" state, then flatten into a plateau once the diode is fully blocking and further bias no longer changes the (already very high) isolation resistance.
- This curve is, physically, just the PIN diode's own resistance-vs-bias transfer characteristic viewed through the modulation-depth formula — say so explicitly in your report, since it directly explains the S-shaped/saturating curve you should observe.

### (d) Response law of the detector — detailed procedure and analysis

**Method (i) — square-modulated microwaves, VSWR meter:**

1. Set up MOD mode with a clean, strong square wave as in parts (b)/(c).
2. Connect the crystal-detector mount (in place of the CRO) to the VSWR meter input.
3. Set the VSWR meter's gain control near maximum, but confirm the needle isn't pinned/overloaded.
4. Slide the slotted-line probe to the position of maximum field and **lock it there** for the whole of Method (i) — unlike Method (ii), you are not scanning probe position here, only the attenuator.
5. Adjust the attenuator and VSWR fine-gain together so the needle sits at a convenient full-scale reference mark (e.g. 0 dB) on the VSWR dB scale, at your starting attenuator setting. Note both readings.
6. Increase the attenuator in small, regular steps; at each step record (i) the attenuator micrometer reading and (ii) the resulting VSWR-meter dB reading. Continue until the reading approaches the bottom of the scale or becomes noise-limited.
7. Change the VSWR gain to a new fixed (lower) setting — this samples a different absolute input power range — and repeat steps 5–6. Do this at 3–4 different gain settings so you can compare the exponent n across different power levels, per the manual's instruction to determine n "at several VSWR gains in decreasing order."
8. For each gain-setting run, plot attenuator setting (x, decreasing to the right, per the manual's convention) vs. VSWR-meter dB reading (y) on a **linear** graph; a straight-line region with slope ≈ 1 confirms square-law behavior (n ≈ 2) in that power range, and you'll see the slope departing from 1 wherever the detector leaves the square-law region.

**Method (ii) — unmodulated CW microwaves, micro-ammeter, standing-wave scan:**

1. Switch the Gunn power supply's modulation selector to **CW** (not MOD).
2. Physically **disconnect the PIN bias cable** from the PIN diode — the PIN modulator plays no role in this method.
3. Connect the detector output to a **micro-ammeter** instead of the VSWR meter — with an unmodulated CW carrier there is no audio-rate envelope for the VSWR meter's tuned amplifier to lock onto, so you must read the detector's raw rectified DC current directly.
4. With the slotted section shorted at the far end (as already wired), slide the probe along the line while watching the micro-ammeter to locate a sharp current **minimum** (null) — fine-adjust position for the lowest, most well-defined reading; minima are much easier to locate precisely than the broad maxima.
5. Take this minimum's position as your reference x = 0. Record the (near-zero) current here.
6. Move the probe in small, regular steps toward the neighboring maximum, recording probe position (from the slotted line's built-in position scale) and micro-ammeter current at every step. Continue past the maximum to the following minimum if you can — a full period gives a much more robust λg than a single half-period.
7. From the spacing between successive minima, compute λg = 2 × (spacing between adjacent minima) — remember the *current/power* pattern has period λg/2, so don't mistake this spacing directly for λg without the factor of 2 (see §4.4 above).
8. For every recorded (position, current) pair, compute sin(2πx′/λg), where x′ is the distance from your chosen reference minimum. Plot log(I) vs. log|sin(2πx′/λg)| on log-log paper (or take logs numerically and use a linear graph) — the slope of the resulting straight line is n.
9. Repeat the full scan at several different attenuator settings between 0 dB and 10 dB, exactly as the manual specifies, to map how n changes with input power — and compare against your Method (i) results at similar power levels as a consistency check.

### (e) Calibrate the micrometer attenuator — detailed procedure and analysis

**Step-by-step procedure:**

1. Set up MOD mode with a clean square wave on the CRO, as in part (c), with **both** the calibrated (reference) attenuator and the uncalibrated micrometer attenuator present in the chain.
2. Set the calibrated attenuator to 0 dB.
3. Adjust the CRO vertical gain so the square wave's top level sits at a clear, easily-read graticule division, with full-scale (unclipped) deflection.
4. Record this top-level reading as your reference.
5. Increase the **calibrated** attenuator in known steps (e.g. 1, 2, 3, … dB); at each step, record the new (reduced) top-level reading of the square wave on the CRO — this is simply mapping out how the CRO trace height falls with known, exact dB increments.
6. Return the calibrated attenuator to 0 dB.
7. Now, instead, adjust the **uncalibrated** micrometer attenuator away from its reference position until the square wave's top level reproduces the *exact* reduced level you recorded in step 5 for a chosen dB value.
8. Record the micrometer reading that reproduces that level — this micrometer position is now known to correspond to that dB value.
9. Repeat steps 7–8 for each dB step recorded in step 5, building up a table of (micrometer reading, equivalent dB).

**Analysis:**
- Plot micrometer reading (x-axis) vs. attenuator dB (y-axis). Do **not** force a straight-line fit through this data: mechanical vane insertion depth versus dB attenuation is not linear in general (attenuation grows roughly exponentially with penetration into the field region for a resistive-vane attenuator), so present it as a smooth calibration curve (or simply a lookup table) and note the range over which the mapping is monotonic and usable.

### (c) PIN modulator: % modulation vs bias
- Compute %m = (Vmax − Vmin)/(Vmax + Vmin) × 100 for each PIN bias setting (using CRO-read Vmax/Vmin relative to true zero, found by grounding the input).
- Plot %m vs PIN bias; expect %m to rise from near 0% (PIN fully ON/low resistance — signal barely chopped) toward a maximum as PIN bias approaches/crosses into the regime that fully shuts the line off, then plateau once the PIN is fully "OFF" and further bias no longer changes its (already very high) blocking resistance.
- This curve is essentially probing the PIN diode's own resistance-vs-bias transfer characteristic, mapped through the modulation depth.

### (d) Detector response law
- **Method (i)**: linear graph of attenuator-dB setting (x, decreasing to the right per manual convention) vs. VSWR-meter dB reading (y); the local slope of this curve at each region gives you n at that power level. Do this at several VSWR gain settings (i.e., several absolute power levels) to map out the power range over which n = 2 holds and where it starts to depart.
- **Method (ii)**: log-log plot of raw detector current I vs. sin(2πx/λg) [i.e. sin βx, using your measured λg from §4.4 and probe displacement x from a chosen reference minimum]; slope of the log-log line = n. Repeat at multiple attenuator settings (0–10 dB) — build a small table of n vs. power level and comment on the trend.
- **Key check**: both methods should agree (within experimental uncertainty) on n ≈ 2 at low power, and both should show n departing from 2 as power increases — cross-validating the two independent methods is itself a good analysis point to include in your report.

### (e) Attenuator calibration
- Plot uncalibrated-micrometer reading (x) vs calibrated-attenuator dB reading (y) — this curve is your calibration; note it will likely be **nonlinear** (mechanical vane insertion depth vs. dB attenuation is not a linear relationship in general, since attenuation depends roughly exponentially on penetration depth into the field region), so don't force a linear fit — present it as a smooth calibration curve (or a table) and note the useful/monotonic range.
- This calibration curve is what lets your uncalibrated attenuator be used quantitatively (e.g., in Method (i) of part (d), if you ever needed to translate its readings into dB directly).

---

## 6. Error Analysis

### 6.1 Instrumental/reading errors
- **DVM/voltage readings**: ± (last-digit resolution), typically ±0.01–0.001 V depending on meter; report the manufacturer's stated accuracy if available, else use half the smallest division.
- **Current (DPM) readings**: similarly ± last-digit or a small % of reading (many digital panel meters specify accuracy as e.g. ±0.5% of reading ±1 digit) — use the same convention consistently and say so in your report.
- **Micrometer positions** (Gunn tuning, attenuator, and slotted-line probe position): mechanical micrometers typically resolve to 0.01 mm, but **backlash** (mechanical play when reversing direction of travel) is often the dominant source of positional error, larger than the nominal resolution — always approach a target reading from the *same direction* consistently to minimize this, and note it explicitly as a systematic-error source in your write-up.
- **CRO voltage/timing readings**: limited by vertical-scale graticule resolution and by your ability to judge the "true zero" level (found by grounding the input) — repeat the ground-reference check periodically, since scope DC offset can drift.

### 6.2 Systematic effects specific to this experiment
- **Thermal drift of the Gunn diode**: self-heating during a measurement run shifts the effective I-V curve slightly over time (mobility and threshold field are temperature dependent) — this is a reason to take data reasonably quickly and to avoid dwelling at V₀ (also protects the diode, as the manual stresses).
- **Detector mount coupling / probe depth**: the slotted-line probe, while designed to minimally perturb the field, is never perfectly non-invasive — deeper probe insertion increases sensitivity but also increases perturbation of the standing wave pattern (i.e., the coupling itself introduces a small effective loading). Keep probe depth fixed while scanning position for the λg determination.
- **Non-ideal isolator/mismatches**: residual reflections from imperfect isolator isolation or connector mismatches (VSWR of the connectors/adapters themselves, not just the diode under test) will imprint small ripples on what should be smooth I vs x or I vs attenuator-dB curves — if you see small periodic ripples riding on your main trend, this is the likely cause, not random noise.
- **Non-square-law departure treated as "error" rather than physics**: be careful in your report to distinguish genuine measurement error/uncertainty from the *physical* deviation of n from 2 at high power — the latter is the actual physics result of part (d), not a mistake to be corrected for.
- **Frequency meter reading vs. actual oscillation frequency**: the resonant absorption dip has some finite width (loaded Q of the wavemeter cavity), so there's an inherent resolution limit (typically a few MHz) in reading the dial at the dip center — take the reading at the point of maximum dip (minimum indicator reading), which is more precisely locatable than the dip's edges.

### 6.3a Error analysis specific to the fluctuating/oscillatory region of part (a)

This deserves separate treatment because it is a genuinely different *kind* of error from everything else in §6.1–6.2, and treating it the same way (e.g. quoting ± half a meter division) will understate your real uncertainty by an order of magnitude.

- **The dominant error source in region 4 is physical, not instrumental.** Across your data, the spread of current readings clustered near a given nominal voltage in the oscillatory region is on the order of 50–90 mA (e.g. multiple entries near V ≈ 3.2–3.7 V ranging roughly from 320 to 411 mA), whereas your DPM's own resolution/last-digit precision elsewhere in the table is ~1 mA. That is a ~50–90× larger spread than instrumental resolution alone would predict — the correct conclusion is that the dominant uncertainty here comes from the diode's own dynamical instability (domain-transit oscillation, nucleation jitter, occasional mode-hopping, and bias-circuit relaxation oscillations feeding back through the supply's output impedance — see the region-4 discussion in §5(a) above), not from you or your meter.
- **How to quantify it:** for each nominal voltage setting in the fluctuating region, take **several repeated readings over a short time window** (rather than a single instantaneous reading) if your bench time allows, and report the **mean ± half the observed range** (or the sample standard deviation, if you have ≥4–5 repeats) as your error bar — not the meter's nominal resolution. If you only have one reading per nominal voltage (as in a single sweep), you can still estimate a representative fluctuation magnitude by looking at the scatter of nearby points within a small voltage window and using that local scatter as your error-bar size for the whole region — state clearly in your report which of these two approaches you used.
- **Do not attempt a single differential-resistance (dV/dI) number for the whole oscillatory region.** Because the instantaneous operating point is not well-defined there (the system genuinely visits a range of quasi-states rather than sitting at one), a single slope number would misrepresent the physics. Instead report the region descriptively: identify its voltage extent, its mean current level and typical spread, and note that no single differential resistance is meaningful there — this is itself a legitimate, complete answer for that part of your I-V characterization.
- **The higher-bias smooth region (above where fluctuations die out) can and should be treated normally**: fit it with a simple linear (or low-order polynomial) least-squares fit, and quote a slope with its standard error exactly as you would for the low-field Ohmic region — this is real, repeatable quasi-static behavior of the sustained-oscillation state, not noise.
- **If you also observe the bias voltage itself fluctuating** (not just current) at fixed knob position in this region, note this explicitly as evidence for the bias-circuit relaxation-oscillation mechanism described in §5(a), and — if you want a cleaner practical fix for future runs — mention that inserting a small series resistance or ensuring the power supply's output impedance is low relative to the diode's negative differential resistance is the standard mitigation, though for this write-up simply documenting and quantifying the fluctuation is the expected deliverable.

### 6.3 Propagating errors into derived quantities
- For **d (thickness) estimated from Eth = V₀/d** (Question 2): since d = V₀/Eth, and you're given Eth as a fixed constant (36 kV/cm) with your own measured V₀ carrying uncertainty δV₀ (from §5(a)), propagate simply as δd/d = δV₀/V₀ (assuming Eth is treated as exact, as given in the question) — quote d with this propagated uncertainty.
- For **slope-derived quantities (n from log-log plots, calibration slopes)**: use a proper linear least-squares fit (not just "eyeballing" two points) and quote the slope's standard error from the fit — this is far more defensible than reading two points off a graph, and you should do this explicitly for every log-log or linear plot in this write-up (parts (d) and (e) especially).

---

## 7. Working Through the End-of-Experiment Questions

**Q1. Define decibel (dB). What is 3 dB attenuation? What is 80 dB noise level?**
dB is a logarithmic ratio: for power, dB = 10 log₁₀(P/Pref); for voltage/field amplitude (assuming equal impedances), dB = 20 log₁₀(V/Vref). "3 dB attenuation" means the output power is reduced to 10^(−3/10) ≈ 0.501 of the input — i.e., **almost exactly half power** (this is why 3 dB is the standard definition of a filter's "half-power point" / bandwidth edge). "80 dB noise level" refers to a noise power (or signal-to-noise figure, depending on context) expressed as 80 dB relative to some reference — e.g., if it's a noise floor 80 dB below a reference signal, that's a power ratio of 10⁻⁸, i.e., the noise power is one hundred-millionth of the reference; state your reference explicitly when you answer, since "80 dB" is only meaningful relative to a stated reference level.

**Q2. Estimate GaAs slab thickness from I-V curve and Eth = 36 kV/cm.**
Use d = V₀/Eth with your measured V₀ from part (a) (§5(a) and §6.3 above give you the method and how to propagate the uncertainty). Expect d in the tens-of-microns range (consistent with the manual's own example of d ≈ 10 μm giving ~10 GHz via v/d).

**Q3. Can Gunn bias amplitude-modulate the microwaves? If yes, how?**
Yes — since the Gunn oscillator's output power (and to some extent frequency) depends on the DC bias point within/around the NDM region, modulating the bias with a low-frequency signal will amplitude-modulate the microwave output. However (as discussed in §3.3), this simultaneously perturbs the operating point and hence the *oscillation frequency* (AM-to-FM conversion / "chirp"), which is why a dedicated, frequency-transparent element (the PIN diode, downstream of an isolator) is preferred for producing clean AM in practice — direct Gunn-bias modulation is possible but not "clean."

**Q4. Why an isolator/attenuator between Gunn and PIN?**
Covered in depth in §3.2: prevents reflected/mismatched power from the switching PIN diode (and beyond) from being fed back into the Gunn oscillator, which would perturb its frequency/amplitude stability (load pulling) and couple the modulation process back into the source itself, corrupting clean characterization of both devices independently.

**Q5. Why is the slot in the center of the broad wall? What if it were in the narrow wall?**
Covered in §3.6: the broad-wall centerline is where the TE₁₀ mode's longitudinal wall-current component is zero, so a slot there doesn't cut across current flow lines and hence doesn't radiate or perturb the internal field — allowing genuinely minimally-invasive field sampling. A slot in the narrow wall (or off-center on the broad wall) would cross non-zero current components, causing power leakage/radiation from the slot, distortion of the internal field distribution, and an unreliable, lossy probe measurement.

---

## 8. Quick-Reference Summary Table

| Component | Physical mechanism | Role in this experiment |
|---|---|---|
| Gunn diode | Inter-valley (Γ→L) electron transfer → NDM → traveling high-field domains | Microwave source via transit-time oscillation, ν = v/d |
| Isolator | Ferrite non-reciprocal (Faraday rotation based) device | Protects source from reflections/load pulling |
| PIN modulator | Bias-controlled carrier injection/depletion in i-region → variable RF resistance | Amplitude-modulates carrier without perturbing source frequency |
| Frequency meter | Tunable resonant cavity, weakly coupled | Direct frequency readout via absorption dip |
| Micrometer attenuator | Resistive vane, variable insertion depth | Power control / dynamic range management (calibrated in part e) |
| Slotted section (shorted) | Probe samples standing wave through non-radiating broad-wall-center slot | Measures λg and field/power distribution |
| 1N23 crystal detector | Schottky point-contact diode rectification | Converts RF envelope to measurable DC/AC signal; I=kVⁿ, n→2 at low power |
| VSWR meter / ammeter / CRO | Tuned AF amplifier+meter / DC meter / time-domain scope | Read detector output appropriately depending on modulated vs CW vs waveform-shape needs |

---

*This document is meant as study material to prepare for and write up Experiment 20. Cross-check all derived numerical estimates (V₀, d, n, λg) against your own actual bench data — the physics and reasoning above are general; the numbers are yours to measure.*
