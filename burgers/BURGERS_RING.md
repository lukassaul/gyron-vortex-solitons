# Sink-driven Burgers ring in a two-class spin fluid: exact form and checks

Units: core radius δ = √(4ν_v/κ) (κ = conversion sink rate, ν_v = viscosity), sound speed c, far-field weighted density. Thin ring of radius R ≫ δ; cross-section radius r; s along the core. Script `burgers/ring.py`.

## Mechanism
Conversion removes gyrons from the weighted classes inside the core: ∇·u = −κ. The weighted fluid converges at u_r = −κr/2 (the Burgers inflow, with no axial outflow), and the dilatation term −ω∇·u = +κω replaces vortex stretching. Far away the intake is a 3D sink, Q/(4πr²), screened by replenishment of the weighted classes in the vacuum.

## Exact moments
| Moment | Form |
|---|---|
| Class swirl V±_θ | Γ± (1 − e^{−r²/δ²})/(2πr) (Burgers; steady with the sink) |
| Mass drift U_θ, correlation flux W_θ | Γ_U f, Γ_W f with Γ_U = (Γ₊ + Γ₋)/2, Γ_W = (Γ₊ − Γ₋)/2, f = (1 − e^{−r²/δ²})/(2πr) |
| Canonical vorticity per class | ∝ Γ± e^{−r²/δ²}/(πδ²) (Gaussian) |
| Intake U_r | −κr/2 in the core; Q/(4πr²) far away |
| Flows along the core U_s, W_s | 0 in steady state (see below) |
| Weighted density | n_w = exp(−(U² + W²)/2c²) (isothermal Bernoulli) |
| Net spin pressure | b = −U·W/η (screened overlap) |
| Charge | holding outside r_c: λ_q = −2π r_c b′(r_c), Q = 2πR λ_q ∝ Γ_UΓ_W R/(η r_c²) |

## Integrals
- Energy (= mass): E = (Γ_U² + ηΓ_W²) R [ln(8R/δ) − 1.615] + intake (small, ∝ κ²). Finite.
- Impulse: P = πR²Γ_U. Self-propulsion U_self = (Γ_U/4πR)[ln(8R/δ) − 0.558].
- Spin J = 0, magnetic moment = 0 (no flow or current along the core).
- Dissipation D = ν_v(Γ_U² + ηΓ_W²)R/δ² ~ κE/[4 ln(8R/δ)]: supplied by the intake (the engine).

## Results
| Case | Γ₊, Γ₋ (cδ) | R/δ | E | Q | Impulse | U_self/c | max speed |
|---|---|---|---|---|---|---|---|
| Static paired (n₊ = −n₋) | 3, −3 | 20 | 935 | 0 | 0 | 0 | 0.30c |
| Moving, U only | 3, 3 | 20 | 623 | 0 | 3770 | 0.054 | 0.30c |
| Moving, charged (n₊ = 2n₋) | 4, 2 | 20 | 727 | −16.4 | 3770 | 0.054 | 0.41c |
| Same, larger ring | 4, 2 | 60 | 2874 | −49.2 | 33929 | 0.022 | 0.41c |
| Near the choke | 9.5, 4 | 20 | 3939 | −101 | 8482 | 0.121 | 0.96c |

## Stability
1. **Core (Burgers attractor):** stable; a perturbed core relaxes back at rate ~κ (Burgers-vortex stability, Leibovich–Holmes, Robinson–Saffman, from memory).
2. **Centrifugal (Synge, variable density, per class):** stable in all cases; the circulation rises monotonically outward (open circulation, no return layer).
3. **Shear sheets:** none.
4. **Subsonic:** peak class speed 0.64Γ±/(2πδ) < c, so |Γ±| < 9.85 cδ (the winding cap).
5. **Circulation leakage across the axis:** ~exp(−R²/δ²), negligible for R ≳ 5δ.
6. **Ring radius:** moving rings (Γ_U ≠ 0) keep R through impulse conservation (Kelvin ring family): stable. The **static paired ring collapses**: its energy rises with R (dE/dR > 0) and nothing resists, since it has no impulse and no spin.

## Findings
- The sink sustains swirl (circulation) but **drains flow and current along the core**: the only steady along-core profile in a converging flow is E₁(r²/δ²), log-singular on the axis, i.e. it needs a line source of axial momentum; removed gyrons carry their momentum away. Spin and moment therefore decay at rate ~κ unless driven along the core.
- Sink-driven Burgers rings are therefore **moving Kelvin-type rings**: stable, finite energy, charged when both U and W swirl (Q ∝ Γ_UΓ_W R/(η r_c²), growing with ring size and dependent on where holding sets in), with no spin and no moment.
- A static ring (zero impulse) needs an outward term that the sink cannot supply: spin (a driven along-core flow) or another conserved quantity.
- Open: an along-core drive (chirality of conversion, active stress) that sustains spin and moment; the charge's dependence on R and r_c.
