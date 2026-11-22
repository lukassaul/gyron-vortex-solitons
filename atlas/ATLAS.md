# Atlas of stable vortex structures in a spin fluid

General, dimensionless; no calibration to any particular particle.

## Units and groups
- Length: the choke radius r_i = Γ/(2πc), where a free vortex of circulation Γ reaches the sound speed c (the hollow-core radius). Velocity c; density: far-field weighted density.
- A ring of winding n is the n = 1 ring rescaled by r_i, so winding enters only through the length unit (and through any other fixed length in the medium).
- Groups: P (vacuum confinement pressure / ρc²), S (surface tension / ρc²r_i), J~ = J/(ρc r_i⁴) (spin), N = C r_i²/ν_v (Carnot number: engine rate over viscous rate), Re = Γ/ν_v, F (fraction of the circulation that reaches outside the core).

## Page 1: driven (Carnot) rings
Model (`atlas/carnot_ring.py`): thin ring of radius R; cross-section with hollow core (w < r_i), free vortex out to the core edge a, return layer that cancels a fraction 1 − F of the circulation, exterior irrotational flow carrying FΓ (ring energy ½(FΓ)²R[ln(8R/w_out) − 2], impulse πR²FΓ); toroidal spin flow carrying J; weighted density from the isothermal Bernoulli relation. Engine W = C·E, viscous loss D (shear of both flows, exterior potential flow included). Ring radius in mechanical equilibrium at fixed J (fast); core edge a by the slow budget E′(a)ȧ = W − D.

Conditions: balance; Carnot-stable ((W − D)′/E′ < 0, i.e. balanced core below the passive energy minimum); subsonic; thin (R/w_out ≥ 2); edge layers stable against centrifugal (Synge, variable density) and spin-jet shear instabilities, which are damped viscously below a critical Re_c. Scan: F ∈ {0, 0.25, 0.5, 0.75, 0.9, 1}, P ∈ [10⁻³, 1], S/P ∈ {0.3, 3}, J~ ∈ [10², 10⁸].

### Results
1. **Confined circulation (F = 0) gives no stable driven ring.** The free-vortex energy between the choke and the core edge grows as π ln(a/r_i), so the passive energy rises with core size; any engine surplus inflates the core. Stable driven rings need open circulation: F ≳ 0.25 (8/70 points), F ≥ 0.5 (~45% of points).
2. **Open circulation means the ring moves.** Stable rings carry impulse πR²FΓ and self-propel at U_self ≈ 0.005–0.12 c. A driven ring at rest needs a partner that cancels its impulse (a counter-circulating second class or a coaxial ring).
3. **Weak confinement, large spin.** Windows exist for P ≲ 0.03 and J~ ≳ 10³–10⁴; they widen as P falls and J~ rises (up to a factor ~15 in core size at P = 10⁻³, J~ = 10⁸).
4. **Carnot number of order one.** N = C r_i²/ν_v ≈ 0.3–7: the engine rate must be comparable to the viscous rate across the choke radius.
5. **Viscous regime.** Edge layers stay stable only below Re_c ≈ 10²–3×10²: for F < 1 the return layer is centrifugally unstable (circulation falls outward) above Re_c; for F = 1 there is no return layer and the limit comes from the spin-jet shear. Stable driven rings are low-Reynolds structures.
6. **Profiles** (`page1_examples.png`): the weighted density dips toward the hollow core (Bernoulli), the spin flow peaks a few r_i out, and the poloidal flow falls as 1/w (F = 1) or faster (F < 1).

Figures: `page1_map.png` (stable fraction vs F; window width over P and J~ for F = 0.5 and 1), `page1_examples.png`. Data: `carnot_ring_atlas.json`.

### Caveats
Thin-ring, single-profile family; edge-stability thresholds are estimates (growth rate vs lowest viscous mode); no field sector or spin–vorticity locking yet; the engine acts on total energy (linear conversion).

### Next pages
- Page 2: impulse-free driven structures (two counter-circulating classes; coaxial ring pairs).
- Page 3: lines and helices (per-length Carnot balance; kink with closure).
- Page 4: conserved-sector (ideal) structures and the transport conditions for their persistence.
- Page 5: composites (shared cores, linked and Borromean rings) and the windings they allow.
