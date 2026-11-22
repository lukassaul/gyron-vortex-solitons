# Soliton programme, step 1: energy functional and exact ring

Static, axisymmetric. Read with `claude/soliton-programme-brief.md` and `claude/gyron-context-brief.md`. Goal (third paper): exact solutions of the gyron transport equations for a ring vortex, from exact-solution functionals. Circulation is imposed: the state is assumed, and only its persistence is required. Solutions are needed for both handednesses and both charges.

## 1. Fields
- Bulk (Madelung): ψ = √(ρ/ρ_m) e^{iθ}, v = ζ∇θ, ζ = ħ/m_e.
- Weighted populations: 𝒦± ≥ 0 with drifts V±; net 𝒦 = 𝒦₊ − 𝒦₋, F = 𝒦₊V₊ − 𝒦₋V₋; vacuum 𝒦₊ = 𝒦₋ = 𝒦₀ (net 𝒦 = 0).
- Locking (no independent field degrees of freedom): Φ_e = g_Φ𝒦, A = g_ΦF/c², so E = −g_Φ∇𝒦, B = (g_Φ/c²)∇×F.

## 2. Functional
E = E_bulk + E_Q + E_w + E_field + E_ph + E_ent

- E_bulk = ∫ [½ρv² + 𝒰_b(ρ)], 𝒰_b = (p₀/(γ−1))[(ρ/ρ_m)^γ − 1 − γ(ρ/ρ_m − 1)], γ = 7/5, p₀ = ρ_mc²/γ.
- E_Q = ∫ K_H|∇√H|² = ∫ (ζ²/2)|∇√ρ|² (H = X₀ρ).
- E_w = Σ± ∫ [(e± + p±)γ±² − p± − e₀± − μ₀±(𝒦± − 𝒦₀)], γ± = (1 − V±²/c²)^{−1/2}. Closure e(𝒦), p(𝒦) pending. Isothermal limit (sound speed c): ½ρ_wV² + c²ρ_w[ϰ ln ϰ − ϰ + 1], ϰ = 𝒦/𝒦₀.
- E_field = ½ε*g_Φ² ∫ [|∇𝒦|² + c⁻²|∇×F|²] (via locking: a Korteweg-type stiffness of the weighted fluid; the second term is |S⁽²⁾|²/2I₂ with ε* = I₂/c²).
- E_ph = λ_ph ∫ (𝒦/𝒦₀)(p_b − p₀): cross-section excess coupled to the bulk pressure deficit.
- E_ent = ½ε_ent ρ_w ∫ ϰ₊ϰ₋|V₊ − V₋|²: conservative counterpart of inter-population drag.
- Holding: the static Coulomb condition C_F = c²∇𝒦 derives from the potential E_hold = −βc²∫h u(𝒦), which cancels the acoustic stiffness wherever h = 1. Charge sits where holding fails: ρ_e = −(ε*g_Φ/c²)∇·C_F. Proposal: h = h(n), collisional, h(1) = 1, h(0) = 0.

The functional is used as an energy–Casimir functional 𝔉 = H − Σμ±N± − Σλ_i𝒞_i − U·P: its stationary points are the steady, helicity-relaxed solutions of the transport equations. Dissipative collision parts vanish at a steady state.

## 3. Non-dimensional form
Units ξ = ζ/c, c, ρ_mc²ξ³. n = |ψ|², k± = 𝒦±/𝒦₀, v± = V±/c, k = k₊ − k₋, f = k₊v₊ − k₋v₋.

Ẽ = ∫d³x̃ [½|∇ψ|² + Ũ_b(n) + Σ±ν±𝔈(k±, v±) + ½Γ(|∇k|² + |∇×f|²) + Λk(n^γ − 1)/γ + ½ε̃k₊k₋|v₊ − v₋|²]

- Ũ_b = [n^γ − 1 − γ(n − 1)]/(γ(γ − 1)); ½|∇ψ|² holds bulk kinetic + quantum terms (ζ/ξ = c).
- Isothermal closure: 𝔈 = ½kv² + k ln k − k + 1.

## 4. Ratios
| Ratio | Definition | Status |
|---|---|---|
| Γ | ε*g_Φ²𝒦₀²/(ρ_mc²ξ²) | free (field stiffness) |
| Λ | λ_ph | free (phoretic coupling) |
| ε̃ | ε_ent ρ_w/ρ_m | free; may be set to 0 |
| ν± | ρ_w±/ρ_m | set by vacuum statistics via a = s∥²/2I⊥k_BT |
| Q | e/(ε*g_Φ𝒦₀ξ) | fixed by α |
| h̃ | unit helicity | fixed by closure (q = N_H/N_θ) plus a flux unit (§12) |
| I_s/I₂ | spin-to-weighted moment ratio (§6) | set by vacuum statistics; fixes g |
| γ, closure sound speed | 7/5, c | fixed |

- The field ratio carries 1/ξ².
- α condition: e²/(4πε*ρ_mc²ξ⁴) = ΓQ²/4π = ακ = αẼ_min, which fixes Q.
- Drag is a rate, not an energy: it sets the gradient-flow mobility and damping only.

## 5. Catalog of constraining equations

A. Integral constraints (Lagrange multipliers)
1. Charge: Q = −∮_∞∇k·dS = e/(ε*g_Φ𝒦₀ξ). Multiplier φ_q.
2. Helicity: h̃ = ∫f·∇×f = ℋ₁c²/(g_Φ²𝒦₀²ξ²). Multiplier λ_H (generalized Taylor λ). Possibly two (one per population).
3. Impulse: P̃_z = 0. Multiplier U.
4. Spin: J_z = ħ/2, handled by the spin ansatz (§6) as a boundary condition.

B. Topological constraints (discrete)
5. ∮∇θ·dl = 2π around the core circle (bulk circulation 2πζ).
6. F-line windings N_θ = 1 and the chosen N_H; closure q = N_H/N_θ on occupied tori (§12).

C. Local kinematic constraints
7. Locking (substituted).
8. ∇·(k±v±) = 0 for each population separately.
9. ∇·(n∇θ) = 0.
10. ∇·F = 0 (Lorenz condition → Coulomb gauge in the static case).
11. n ≥ 0, k± ≥ 0, |v±| < 1.

D. Boundary conditions
12. At infinity: n, k± → 1; v±, ∇θ → 0; exterior Coulomb plus dipole.
13. Axis: χ± = 0, I± = 0, ψ regular. Bulk vortex core circle: n = 0.

E. Self-consistency
14. ΓQ²/4π = αẼ_min.
15. κ = Ẽ_min (m_e = κρ_mξ³).

F. Checks
16. Von Laue: ∫Σ_iT_ii d³x = 0.
17. Derrick scaling with constraint terms.
18. g = 2m_eμ/(eJ); with the spin ansatz g = 3m_e/(eμ*I_s).

## 6. Spin ansatz
Linear response: 𝒮(x) = I_s B(x), I_s set by the vacuum statistics through a. Total spin J = I_s∫B d³x = (2/3)μ*I_s μ.
- J = ħ/2 fixes the exterior dipole coefficient of f (boundary condition, no multiplier).
- g = 3m_e/(eμ*I_s): set by medium constants, not the profile. The bulk carries no spin.
- The sign of I_s must follow the charge (§13): the spin polarization is carried by the population that sets the charge.
- Risk: I_s is a parameter unless the vacuum statistics fix it independently.

## 7. Coordinates and representation
Approach: spherical (ρ, θ) with a multipole basis for the field sector; a local frame near the O-ring for the bulk vortex core (matched asymptotics); no global toroidal coordinates (the ring is fat).

- Weighted fluxes (stream form): k±v± = ∇χ± × ∇φ + I±∇φ, so continuity holds identically.
- ∇×f = ∇I × ∇φ − (Δ*χ)∇φ, Δ* = ∇² − (2/r)∂_r.
- Helicity: h̃ = −2∫IΔ*χ/r² d³x.
- Field sector: expand in Chandrasekhar–Kendall (Beltrami) eigenfunctions, a complete set; stability separates into (ℓ, m) problems.
- Bulk core: ψ = √n e^{iθ} with unit winding around the O-ring; inner core of size ~ξ matched to the outer multipole field, giving the ln(8R/a) line tension.

## 8. Grad–Shafranov reduction
Single effective population, k locally uniform:
- δI: ΓΔ*I − (ν/k)I = 2λ_HΔ*χ
- δχ: ΓΔ*χ − (ν/k)χ + 2λ_H I = 0

Taylor limit (ν → 0): I = μχ, Δ*χ + μ²χ = 0, μ = 2λ_H/Γ. With ν: Γμ² − 2λ_Hμ + ν/k = 0, a double-Beltrami state (Mahajan–Yoshida type) with μ₁μ₂ = ν/(kΓ) = ξ_w⁻². The ν terms screen exterior fields at ξ_w; the holding term removes this for the longitudinal (Coulomb) part; the transverse (dipole) part is open.

## 9. Steady ring equations (minimal closure)
Closure: C_F± = c²h(n)∇𝒦± (holding), drag = 0 (persistence; or V₊ = V₋ where populations overlap), C_W = T′ = 0 outside the core. Two-region model: core D_c where n < n_* (h = 0), exterior h = 1.

Core: each population is a compressible Bragg–Hawthorne flow with 𝔅± = ½V±² + c² ln k± = 𝔅±(χ±), L± = rV±,φ = L±(χ±), I± = k±L±:
∇·(∇χ±/(k±r²)) = k±𝔅±′ − (k±/r²)L±L±′,  c² ln k± = 𝔅± − (|∇χ±|² + k±²L±²)/(2k±²r²).
Imposed circulation = choice of 𝔅±, L±; simplest L± = λ±χ±, 𝔅± = 𝔅₀ + β±χ±.

Exterior: ∇²k = 0; Δ*I = 0 (dipole B); Δ*χ = 0 (no toroidal B); bulk irrotational Madelung with unit winding.

Matching on ∂D_c: χ± and normal fluxes continuous; normal stress balance (weighted, bulk, quantum, field); tangential slip allowed.

## 10. Exact solution: low-Mach Beltrami core (λa = π)
Uniform k in the core (V ≪ c), β = 0, symmetric populations, U = 0. Units a = 1, amplitude A.
- Core ρ < a: χ = Aρ j₁(λρ) sin²θ, I = λχ, ∇×F = λF, B = (g_Φ/c²)λF.
- Exterior: χ = (Aa/π)(a/ρ) sin²θ, I = Aa sin²θ/ρ (harmonic; dipole B, no toroidal B).
- Matching I, ∂_ρI continuous: x j₀(x) = 0, so λa = π, D = Aa.
- O-ring (bulk vortex core): x j₀ = j₁ → x = 2.744, r_O = 0.873a.
- Topology: χ(a) = (Aa/π) sin²θ ≠ 0, so most poloidal F-lines are open (enter, wrap the O-ring, leave as a curl-free dipole flux). Closed nested tori only in χ ∈ (0.318, 0.338)A near the O-ring. Alternative branch j₁(λa) = 0 (λa = 4.493, classical spheromak): all F-lines closed, at the cost of a surface current in B_θ.
- Script: `ring_exact.py`.

## 11. Stability and the radial balance
- Exact identity: E_B = ½Γ∫|∇×f|² = (λ/2)Γh̃ = (π/2)Γh̃/a, interior plus exterior (numerical ratio 1.5708).
- Multipole matching: j_{ℓ−1}(λa) = 0 → λa = π, 4.493, 5.763 for ℓ = 1, 2, 3.
- Tilt (ℓ = 1, m = 1): neutral. Helicity is a pseudoscalar, so the equations and constraints are rotation-invariant; rotated solutions are solutions (orientation zero mode; precession with spin). Helicity breaks mirror symmetry only.
- Shape modes ℓ ≥ 2: higher λ, higher energy at fixed helicity → stable. λa = π is the ground state of the β = 0 family; λa = 4.493 relaxes to it.
- Breathing: E_B ∝ h̃/a falls with a → the bare core expands (zero circulation: it disperses, it does not dissipate).
- Radial balance with the bulk vortex ring (circulation 2πζ at r_O):
  E(a) = (π/2)Γh̃/ã + 2π²(0.873ã)[ln(8 × 0.873ã) − 7/4], E″ > 0 (stable minimum), at Γh̃ ≈ 4πã ln(7ã).
- With the flux unit h/e (Γh̃ = πκ/α): ã ≈ 1270 (about 9a₀), κ ≈ 3 × 10⁵ — too large; a ring of size ~ξ needs Γh̃ ~ 10.
- Alternative restoring force without circulation: phoretic self-trapping (a bubble), if Λ has the sign that expels the bulk.
- Not yet covered: modes coupling the core boundary to the bulk; modes with V₊ ≠ V₋.

## 12. Winding quantization: the q = 1 ring
- Helicity is fixed by the windings only when the F-lines close: ℋ = 2N_θN_HΦ_tΨ per linked tube, and closure requires the winding ratio q = N_H/N_θ to be rational on the occupied tori. Single-valuedness of the internal phase carried along F-lines around the bulk vortex ring imposes q(Ψ) = N_H/N_θ (electron: q = ±1).
- λa = π solution: q = 0.825 (nearly flat; lines cover their tori densely); ℋ/(2Φ_tΨ) = 5.25 (most helicity in open flux and exterior). Helicity is not yet quantized there.
- Swirling Hill-type family: χ = [Aρj₁(λρ) + (β/λ²)ρ²] sin²θ, I = λχ (Δ*χ + λ²χ = βr²). Matching I, ∂_ρI at a = 1 forces β = −Aλ³j₀(λa)/3: a one-parameter family in λ (λa = π ↔ β = 0).
- q along the family: λa = 2.6, 2.9, π, 3.4, 3.633, 4.0 → q = 0.65, 0.74, 0.82, 0.91, 1.000, 1.13.
- q = 1 solution: λa = 3.633, β = 2.08A, r_O = 0.847a, q flat to ±0.001 across the closed tori ((1,1) lines). Closed tori occupy ~7% of the flux, hugging the O-ring, where the bulk vortex core sits.
- Handedness: λ → −λ is exact (q = −1, helicity reversed). Scripts: `hill_q.py`, `hill_plot.py`.

## 13. Charge at O(V²/c²)
- Population split with the same profile: F± = ±(1 ± δ)F/2, so F₊ − F₋ = F. Each population obeys its own Bragg–Hawthorne equation with the same λ; its Bernoulli function scales with its amplitude squared: 𝔅± = c±²βχ (c± = (1 ± δ)/2), reference 𝔅₀± = 0 on the axis streamline.
- Net core density: k = (δ/c²)[βχ − ½|F|²]. Symmetric populations (δ = 0) give k = 0.
- Exterior: holding (h = 1), ∇²k = 0, k continuous at ρ = a, so Q = 4πa⟨k(a, θ)⟩.
- Results (A = flux amplitude): λa = π: Q = −5.46 δA²/c²; λa = 3.633 (q = 1): Q = −2.67 δA²/c². Identical for both handednesses (k depends on λ², β and |F|² only).
- Charge sign is set by δ (which population circulates faster), not by handedness: all four combinations (±helicity × ±charge) exist at this order. The core now carries an interior E (the seed's E = 0 inside is the δ = 0 limit).
- Link to the spin: with a fixed I_s the spin ansatz gives μ ∥ J for both charges, but μ·J must change sign with the charge. So the spin polarization must be carried by the population that sets the charge: I_s ∝ (the faster population's weight), sign(I_s) ∝ δ. Then charge and spin–moment sign are tied, and the helicity sign remains free: for each charge there are two chiral states (a suggestive 2 × 2 = 4-state structure).
- C symmetry: the population swap is exact only if ν₊ = ν₋; otherwise the medium breaks charge conjugation and the opposite charge is a slightly different solution. Script: `charge.py`.

## 14. Quantized ring: geometry and balance (q = 1, δ = ±1)
- Population split quantized: δ = ±1 (one population carries the flow). Charge discrete (±e, sign = carrying population), C exact if the populations are equivalent, and the carrying population carries the spin. Four solutions: (±e) × (q = ±1).
- With δ = ±1, Q = e fixes A²a (Q = q₀A²a, q₀ = 2.67), so the magnetic energy (∝ A²a) is constant at fixed charge; the radial balance is electrostatic self-energy (out) against the restoring term (in).
- Energy coefficients (a = 1, A = 1, δ = 1): magnetic m₀ = 69.4 (89% inside); c_E = 4.11 (≈100× a point charge's 1/8π: the electric energy is mostly internal polarization of the core). E(ã) = (m₀/q₀)ΓQ̃ + c_EΓQ̃²/ã + E_r(ã), with ΓQ̃² = 4πακ and κ = E_min.
- Restoring terms: picture B (thin core) E_r = 2π²r_O[ln(8r_O/ξ) − 1.615]; picture A (bubble) E_r = p̃₀(4π/3)ã³ + σ̃4πã², p̃₀ = 5/7, σ̃ ≈ 2√2/3 (order-one estimate).
- Fixed by closure (exact): λa = 3.633, r_O/a = 0.847, closed-band half-width b/a = 0.152 (midplane r ∈ [0.697a, a]).

| Γ | Picture | a/ξ | r_O/ξ | b/ξ | κ | Mach | W_B/W_E | M / E / restore |
|---|---|---|---|---|---|---|---|---|
| 1 | B | 1.4 | 1.2 | 0.22 | 154 | 0.99 | 2.4 | 0.63/0.26/0.10 |
| 1 | A | 1.3 | 1.1 | 0.19 | 191 | 1.11 | 1.9 | 0.57/0.30/0.13 |
| 10 | B | 3.0 | 2.6 | 0.46 | 967 | 0.61 | 6.4 | 0.80/0.13/0.07 |
| 10 | A | 2.1 | 1.8 | 0.32 | 1110 | 0.75 | 4.2 | 0.75/0.18/0.08 |
| 100 | B | 7.2 | 6.1 | 1.09 | 7.5×10³ | 0.37 | 17 | 0.91/0.05/0.04 |
| 100 | A | 3.8 | 3.2 | 0.58 | 8.4×10³ | 0.52 | 8.7 | 0.86/0.10/0.04 |
| 10⁴ | B | 52 | 44 | 7.9 | 6.4×10⁵ | 0.13 | 137 | 0.99/0.01/0.01 |

- Validity (corrected): the 'Mach' column is the amplitude A, not the peak speed. Peak |V| = 3.0A (on the axis, z = 0) and the deepest density dip is δk = −4.5A². Low-Mach validity (peak |V| ≲ c, |δk| ≲ 0.3) needs A ≲ 0.25: Γ ≳ 10³ in picture B (a ≈ 19ξ, κ ≈ 7×10⁴), Γ ≳ 10⁴ in picture A (a ≈ 12ξ, κ ≈ 7×10⁵). At Γ ~ 10² (a ≈ 4–7ξ) the peak drift exceeds c: that regime needs the compressible, relativistic solution. Thin core inside the closed band (b ≳ ξ) holds for Γ ≳ 10² in picture B.
- Plots: `ring_profiles.py` (Γ = 10³, picture B): weighted drift streamlines, bulk vortex flow and density, spin streamlines (S = I_sB), net spin pressure K, midplane and axis profiles.
- The restoring term is ~4% of the energy: it sets the size, barely the mass. Mass ≈ locked magnetic energy, κ ≈ 4παΓ(m₀/q₀)² ≈ 62Γ.
- κ ~ 10⁴ in the valid regime (brief assumed κ ~ 1), so ρ_m = m_e/(κξ³) is lower than 1.6 × 10⁷ kg/m³ by κ. κ ~ 1 would need Γ ~ 0.02 (supersonic; requires the compressible solution).
- Caveats: bubble surface tension order-one; weighted kinetic energy counted as field energy (reading (ii)); shape-mode stability of the q = 1 member assumed from the β = 0 family. Scripts: `balance2.py`, `solve_pics2.py`.

## 15. Zero impulse and ring size
- Circulation is fixed by topology (single-valued bulk phase: n·2πζ), not by the flow through the ring. The bulk flow through the hole is the ring's self-induced speed U_self ≈ (ζ/2R)[ln(8R/ξ) − ½] (≈ 0.1c at R ~ 15ξ).
- A static ring needs zero total impulse. The bulk ring's Lamb impulse ρ_m·2πζ·πr_O² must be cancelled by the weighted ring's impulse, which is fixed by its exterior dipole coefficient: P_w = 4πG(1)·νρ_m·A·a³ (G(1) = 0.364 for q = 1).
- Consequences: (1) the weighted and bulk circulations around the O-ring are opposite (fixes that discrete choice); (2) ν·A·(a/ξ) = πr_O²/(2G(1)) = 3.09. With the charge condition (A²a fixed by Q = e) this replaces Γ by ν (weighted inertia fraction, from vacuum statistics) as the last free ratio.
- Sizes (picture B, q = 1, δ = ±1):

| ν | a/ξ | r_O/ξ | κ | Γ (output) | A | peak \|V\|/c |
|---|---|---|---|---|---|---|
| 1 | 10.0 | 8.5 | 1.6×10⁴ | 224 | 0.31 | 0.93 |
| 0.3 | 118 | 100 | 3.7×10⁶ | 6×10⁴ | 0.087 | 0.26 |
| 0.1 | 1070 | 900 | 4.2×10⁸ | 7×10⁶ | 0.029 | 0.09 |
| 0.01 | 1.1×10⁵ | 9×10⁴ | 6.5×10¹² | 10¹¹ | 0.003 | 0.01 |

- Roughly a ∝ ν⁻²; a ring ~10ξ needs ν ~ 1. No solution found for ν = 3 or 10⁻³ in the search range.
- Stability: breathing stable (E″ > 0 at fixed charge); a radius perturbation unbalances the impulse (P_w ∝ a^{5/2}, P_b ∝ a²) and makes the ring translate; relaxation back needs the dynamic step. Tilt neutral.
- Spectrum: each (n, q = N_H/N_θ, branch) gives its own stable size and mass; in the magnetic-dominated limit the medium ratio largely cancels from mass ratios (test against m_μ/m_e ≈ 207).
- Caveats: Galilean Kelvin impulse for the bulk; the weighted impulse assumes real weighted inertia, which must be reconciled with reading (ii). Plots: `ring_full.py` (full meridional sections and 3D streamlines, bulk and weighted). Script: `impulse.py`.
- Open: the bulk is strongly disturbed near the core (|u| ~ c within a few ξ, falling below 10⁻²c by 1.5a), though it carries only ~4% of the energy. This comes from choosing a bulk vortex as the restoring force and impulse ballast, not from charge or closure; see the n = 0 alternative.

## 16. n = 0 test: no static charged ring with a magnetic moment
- Same exterior coefficient: with I = kL(χ), the magnetic dipole and the weighted impulse both come from the exterior poloidal flux χ(a) = G(1)sin²θ: dipole D = λG(1)Aa³, impulse P_w = 4πG(1)νρ_mc·Aa³ (net mass flux of the split ∝ δ). So μ ≠ 0 ⇔ open flux ⇔ P_w ≠ 0 (for δ ≠ 0).
- Fully closed weighted ring (χ(a) = 0): the axis and boundary lie on the same streamline, so I(a) = I(axis) = L(0) = 0 and the exterior dipole vanishes; with the spin ansatz J = I_s∫B, the spin vanishes too. Within the Hill family closed + matched occurs at λa ≈ 5.76, 9.09 (no moment).
- Charge ∝ δ and impulse ∝ δ, so the weighted sector alone cannot have Q ≠ 0, μ ≠ 0 and P = 0. A static charged, spinning ring needs a counter-impulse from the bulk; without circulation a static bubble's potential flow carries none. Bulk circulation (n ≠ 0) is required, consistent with "motion through the bulk is inertia".
- Alternative escape: the ring swims at U with total impulse zero (P_w = ½ρ_mV_bub U for a bubble, U/c ≈ 2.2νA); excluded for co-moving bound electrons unless U is tiny.
- Payoff: with P_w = P_b = ρ_m·2πζ·πr_O², the dipole coefficient is D = πλξr_O²/(2ν), independent of the flow amplitude A: the magnetic moment is fixed by the bulk circulation quantum, the ring size and ν. With the spin condition (J = ħ/2 ∝ I_sμ) this gives a further relation between ν, a and I_s.

## 17. Spin condition with impulse balance: the moment overshoots
- Moment from the locked field: A = g_ΦF/c² gives μ = 4πε*c g_Φ𝒦₀ D_f = 4πce D_f/(Q̃ξ). With the impulse balance, D_f = 4.09ξa/ν, so μ/μ_B = 102.9 ã/(Q̃ν). With the spin ansatz, J = ħ/2 and g = 2 together mean μ = μ_B.
- Charge (Q̃ = q₀A²ã), impulse (νAã = 3.09) and spin (μ = μ_B) give, without any energy input: A² = 38.6/ν, a/ξ = 0.50/√ν, a/ξ = 0.080A. Subluminal flow (peak |V| = 3A < c) needs ν ≳ 10³ and then a ≈ 0.016ξ: far below the healing length, outside the thin-core picture. A ring of size ≳ ξ needs A ≳ 12 (peak flow ≫ c).
- In the valid low-Mach regime (§15) the ring's moment exceeds μ_B by ~250–800 (ν = 1 ring: ≈ 400), i.e. g ≈ 2μ/μ_B ≈ 800 if J = ħ/2.
- General reason (model-independent): a soliton of size a whose mass is mostly magnetic field energy has μ ≈ 2(a/ξ)^{3/2}α^{−1/2} μ_B (23μ_B at a = ξ, 740μ_B at 10ξ). μ = μ_B with W_B ≈ mc² needs a = (α/4)^{1/3}ξ ≈ 0.12ξ; at a ~ ξ the magnetic energy must be only ~(α/4)(ξ/a)³ of the mass.
- Conclusion: with g = 2 and a ring of size ~ξ, the electron's mass cannot be locked magnetic energy (the analogue of "the spin cannot be field angular momentum"). The mass must sit in the weighted fluid's inertia and enthalpy ((e+p)/c², relativistic), in the bulk, or in the quantum stress — which reading (ii) removed. Reading (ii) is therefore in tension with g = 2.
- Next: rebuild the balance with the weighted fluid's relativistic enthalpy as the main energy store and the locked field energy at the ~α level; check whether charge + impulse + spin then admit a ring of size ~ξ with subluminal flow.

## 18. The obstruction is kinematic, not energetic
- The three conditions of §17 use no energy: charge Q̃ = q₀A²ã (Bernoulli response of the weighted density), moment μ ∝ D_f (dictionary A = g_ΦF/c²), impulse νAã = 3.09 (weighted inertia vs bulk ring ρ_m·2πζ·πr_O²). Moving the mass into the weighted enthalpy changes κ, not these relations; relativistic inertia only rescales ν (ν_eff = νhγ²), and ν cancels from the size.
- General result: ã = 0.030 · n · q₀ · A · (r_b/r_O)², with n the bulk circulation quanta, r_b the bulk ring radius, q₀ the charge per A²a (2.67 at low Mach, sound speed c). Then ν = 102.9/(q₀A²) follows.
- A ring of size ~ξ with subluminal flow (A ≈ 0.3) needs n·q₀·(r_b/r_O)² ≈ 110, i.e. ~40× the present value. Levers:
  1. Charge per flux q₀: δk = −V²/2c_w², so q₀ ∝ (c/c_w)². A weighted sound speed c_w ≈ c/6 inside the core would do it. The common-speed closure fixes c_w = c in vacuum only; the paper already expects the closure to change character near and below the healing length.
  2. More bulk circulation (n ~ 40) or a larger bulk ring: unattractive for the electron.
  3. Compression of the weighted fluid (k ≫ 1 in the core) letting F = kV exceed 𝒦₀c at V < c; this also changes the charge relation and needs the compressible solution.
- Status: with g = 2, spin ħ/2, charge e and a static (zero-impulse) ring, the low-Mach, c_w = c closure gives a ≈ 0.08Aξ ≪ ξ. The core closure (c_w, compression) is the deciding physics, not the energy bookkeeping.

## 19. What the kernel requirements decide
- Weighted sound speed in the core: the closure is the leading Chapman–Enskog order and holds wherever the scale ≫ mean free path; ℓ ~ 10⁻³² m ≪ ξ by ~20 orders, so the energy-weighted closure gives c_w = c inside the core as in vacuum (to the multi-messenger 10⁻¹⁵ level). Lever 1 of §18 (c_w ≈ c/6) is excluded. Departures appear only at scales approaching ℓ.
- "No energy terms odd in δ𝒦" (equal masses of opposite charges) forbids a phoretic coupling linear in the net 𝒦. The phoretic coupling must be even (acting on 𝒦₊ + 𝒦₋), so any core compression k_c is symmetric. Then F = k_cV, the charge per flux² falls as 1/k_c, and the size bound becomes ã = 0.030·n·q₀·v·(r_b/r_O)² in terms of the velocity amplitude v: compression does not help (lever 3 excluded). This also settles the charge-symmetry question of §13: population equivalence is a stated kernel requirement.
- Result for the proportional population split (F± = ±c±F) with a bulk-vortex impulse partner: peak speed < c (v ≲ 1/3) gives a ≲ 0.03ξ ≈ 10⁻¹⁴ m (a few classical electron radii). A bulk vortex with core ξ cannot sit on such a ring, so this family is inconsistent: the proportional split plus bulk impulse partner fails.
- Kernel-consistent escape (revises §16): the two populations need not flow in proportion. Impulse comes from the exterior poloidal mass flux ν(χ₊ + χ₋), the moment from the exterior toroidal net flux I₊ − I₋, the charge from the interior ⟨V₊² − V₋²⟩. With χ₋ = −χ₊ outside (counter-streaming poloidally) the weighted impulse vanishes while I_net = (λ₊ + λ₋)χ₊(a) ≠ 0, and charge survives if the interior profiles differ. Then no bulk impulse partner is needed (n = 0 allowed; bulk nearly undisturbed). Requirement: two closed (q = 1) interior profiles with equal exterior coefficients but different interiors, which needs a family beyond linear Hill (nonlinear L(χ)).

## 20. Swirling Hill pair: counter-streaming populations (n = 0)
- Extended family: χ = A[ρj₁(λρ) + (b/λ²)ρ²]sin²θ, I = λχ, with b free (a surface current at the core boundary is allowed; matched b = −λ³j₀/3 is one point). Closure q = 1 gives a curve λa ∈ [3.2, 5.6] (b from 1.43 to 6.99; exterior coefficient G(1) from 0.45 to 0.064; r_O from 0.96a to 0.54a); it crosses the no-surface-current point at λa ≈ 3.63. Script: `swirl_family.py`.
- Pair: population + takes member (λ₊, b₊), population − takes (λ₋, b₋), amplitudes with A₊G₊ + A₋G₋ = 0 (zero weighted impulse: poloidal mass flux has no exterior dipole). Charge Q = 4πa⟨δk₊ − δk₋⟩, δk± = A±²[b±s± − ½|f±|²]; net moment D = λ₊A₊G₊ − λ₋A₋G₋. Spin μ = μ_B gives ã = Q/(8πD) (charge + spin alone; impulse no longer needs ν or the bulk).
- Best pairs: opposite handedness, neighbouring members near λa ≈ 5.2–5.4. They counter-stream poloidally and co-rotate toroidally, so the net moment nearly cancels (D/A ≈ 0.014–0.03) while the charge survives (Q/A² ≈ 1–2). Maximum size at peak drift = c: ã ≈ 0.75 (λa = 5.4 with −5.2); ≈ 0.6–0.7 for neighbouring pairs. About 30× the proportional-split bound (0.027) and of order the healing length.
- No bulk circulation is needed: impulse cancels inside the weighted sector, so the bulk is only mildly compressed.
- The toroidal mass flux of the two populations adds, so the weighted fluid carries angular momentum ν ρ_m ∫ r(F₊φ + F₋φ) d³x: a candidate for the spin itself, in place of the I_s ansatz.
- Open: size from the energy balance (bubble picture, picture A, now consistent with n = 0); relativistic corrections (peak drift ~c at the bound); sensitivity of the small net moment to the pair choice; closure of each population's own lines vs the net lines. Scripts: `pairs.py`, `pair_plot.py`.

## 21. Spin of the pair ring: g as a shape number
- Exterior: a toroidal mass flow falling like the dipole field would give a linearly divergent angular momentum, so outside the core the two populations counter-stream completely (zero mass flux, both components) and carry the net F as a pure counterflow. Each population then has a toroidal vortex sheet at the boundary; the net flux stays continuous. This also means the vacuum dipole field costs no angular momentum.
- Spin = interior weighted angular momentum: J = νρ_m c a⁴ j̃, j̃ = ∫_in (I₊ + I₋) d³x. If the mass is the weighted kinetic energy, κ = ½νã³K̃ with K̃ = ∫(|f₊|² + |f₋|²) d³x (interior plus exterior counterflow).
- Then J = ħ/2 gives ã = K̃/(4j̃), and μ = μ_B gives ã = Q/(8πD); both scale with A. Combined: g = 4πDK̃/(Qj̃), a pure shape number, independent of amplitude, size and ν.
- Range: opposite-handed pairs g ≈ 0.69–16; same-handed pairs g ≥ 15. g = 2 lies on a curve in pair space, e.g. (λ₊, λ₋) = (5.4, −3.6): g = 2.009; (4.6, −4.2): 2.065; (5.2, −3.8): 1.93. Swapping the roles of the two members reverses the charge sign at the same g (the C pair).
- The largest-size pairs (5.4, −5.2) give g ≈ 0.69, so size and g trade off. On the g = 2 pair (5.4, −3.6): ã = 0.94A, and peak drift ≤ c gives A ≤ 0.248, ã ≤ 0.23. Mass κ = ½νã³K̃ ≈ 0.0045ν at the bound, so κ ~ 1 needs ν ~ 200.
- Open: relativistic corrections at the bound; a selection rule on the g = 2 curve (energy balance / bubble); whether the mass is the weighted kinetic energy (vs enthalpy, field and bulk parts). Script: `pair_spin.py`.

## 22. Toward Millikan: charge at the luminal limit
- Requirement: every charge an exact multiple of e, independent of structure (proton = −electron to ~10⁻²¹). In the model ξ is a medium constant (K_HX₀ fixes ζ), so universality means every charged soliton has the same Q̃ = e/(ε*g_Φ𝒦₀ξ), i.e. the same 𝒦-monopole length ℓ_q = e/(4πε*g_Φ𝒦₀). The ring charge Q̃ = q₀A²ã is continuous unless something fixes the amplitude.
- Exactly opposite e⁺/e⁻: automatic. The positron is the same pair with population roles swapped (|Q| identical), and pair creation from a neutral photon conserves charge.
- Candidate for the magnitude: luminal saturation. The carrying population streams at c (peak drift = c), so its Bernoulli deficit saturates and the amplitude is fixed. With the spin conditions (J = ħ/2, μ = μ_B), on the g ≈ 2 pairs:

| pair (λ₊, λ₋) | g | a/ξ | A at peak = c | Q̃ |
|---|---|---|---|---|
| (5.4, −3.6) | 2.009 | 0.233 | 0.248 | 0.0437 |
| (5.6, −3.6) | 1.926 | 0.242 | 0.239 | 0.0455 |
| (5.2, −3.8) | 1.931 | 0.238 | 0.258 | 0.0443 |
| (5.0, −4.0) | 1.859 | 0.246 | 0.268 | 0.0439 |
| (4.8, −4.0) | 2.125 | 0.222 | 0.278 | 0.0372 |
| (4.6, −4.2) | 2.065 | 0.248 | 0.290 | 0.0301 |

- On the main branch (λ₊ = 5.0–5.6) Q̃ = 0.044 ± 3% and a = 0.24ξ ± 3%: nearly shape-independent. The (4.6, −4.2) branch gives 0.030. Near-universal, not exact.
- If this is the mechanism: e = 0.044 ε*g_Φ𝒦₀ξ, a ≈ 0.24ξ, and α = ΓQ̃²/(4πκ).
- Open: exactness (10⁻²¹) needs an invariant, not a 3% coincidence; the luminal regime needs the relativistic solution (does the spread shrink?); the proton (Borromean triple) must reach the same Q̃; fractional quark charges. Script: `millikan.py`.

## 23. Inheritance vs dynamics, and what the lab measures
- Inheritance: 𝒦 conservation makes charge exactly conserved; every route to matter starts neutral, so each event makes exactly opposite charges. Reactions carry the unit between species: n → p + e⁻ + ν̄ gives q_p = −q_e exactly if the neutron is neutral (zero is natural by symmetry). Neutrality of matter (10⁻²¹) then needs conservation, a connected reaction network and one neutral reference, not an independent proton calculation. The proton family must admit charge e, not select it.
- Inheritance cannot make all electrons identical: a photon carries no charge scale, so a continuous family would give each pair its own |q|. Identity (Pauli, Fermi statistics) requires the electron to be an isolated solution: dynamics must fix it.
- Hybrid: dynamics fixes the charge of the lightest stable charged soliton; conservation and the reaction network fix all others as exact multiples.
- Lab frame: charge is measured with apparatus built of the same solitons, so only dimensionless combinations are observable: charge ratios (q_p/q_e, old/new electrons) and α = e²/(4πε*ħc) = ΓQ̃²/(4πκ). The subspace value of e (and of Q̃𝒦₀ξ, ε*, g_Φ separately) may drift with the medium state, unobservable like bulk velocity. The constraint on the galactic-sector "seas" is on α (10⁻⁶–10⁻⁷ per ΔΦ/c²) and on charge ratios, not on subspace e.
- Remaining sharp constraint: an old electron conserves its subspace Gauss flux exactly, while a new one forms at the luminal limit of the local medium. Their ratio is lab-measurable (same apparatus) and must be 1, so medium changes must preserve the saturated electron's subspace flux 𝒦₀ξQ̃_sat (an adiabatic-invariance condition), or old electrons must relax to the new state without changing flux.
- Motion through the medium: with the weighted sector Lorentz-covariant, a moving ring's subspace Gauss flux equals its rest value (Heaviside contraction preserves the flux), so lab charge is independent of the 370 km/s wind.

## 24. Transonic ring: choke and closure
- Weighted fluid is isothermal with sound speed c (Galilean subspace: c is the common wave speed; gyrons move faster, 1.35c). The drift limit is Mach 1: an isothermal flux kV peaks at Mach 1 at e^{−1/2}k₀c with k* = e^{−1/2}k₀ (choked flow). Charge as a choked flux would be fixed by the stagnation state alone, and the choked state is a maximum-flux (extremal) solution, so isolated.
- Compressible Bragg–Hawthorne solved numerically for one population, q = 1 shape (λa = 3.633, matched b), boundary data χ = AG(1)sin²θ, by Picard iteration (script `transonic.py`; matches the exact low-Mach solution to 6×10⁻⁴ at small A).
- Results: the throat is the hole centre on the axis (peak Mach there; density falls toward e^{−1/2}). No steady subsonic solution beyond A ≈ 0.45 (choke). Charge grows faster than the low-Mach law (Q/a = −0.65 at A = 0.43 vs −0.50 extrapolated).
- But compressibility pushes the O-ring outward (0.855a → 0.99a) and the closed q = 1 band vanishes at A ≈ 0.27 (peak Mach ≈ 0.6), before the choke. With fixed shape parameters the quantized ring cannot reach the choke: closure is lost first.
- Next: impose q = 1 and a closed band in the compressible regime by letting the shape (λ, b) move with A, and look for the choked, closed solution (candidate isolated electron). Vortex-sheet nodes at the boundary are excluded from the choke test (sheet treated as a discontinuity).

## 25. Energy balance at fixed spin: isolated solutions and α
- Method: spin is conserved, so fix J = ħ/2 (weighted angular momentum νa⁴j(A) = κ/2 in units ρ_mcξ⁴, with ħ = m_ecξ) and minimize E(A, a) = νa³e_w(A) + Γa e_f(A) + p̃₀(4π/3)a³ + 4πσ̃a² over amplitude and size, with e_w, e_f, j, Q from the compressible solutions (one carrying population, exterior counterflow, q = 1 shape). Scripts: `energy_spin.py`, `minimize_spin.py`.
- Result: for ν ≳ 10 (and ν = 1 with Γ = 100) the energy has an interior minimum below the choke and with the closed band intact (A* = 0.03–0.16 < 0.27). Energy and force balance alone give an isolated solution: one amplitude, one size, one charge Q̃ for given medium ratios, so all electrons identical. For ν ≲ 1 the minimum runs into the choke, where the band is already lost (inconsistent).
- Examples: ν = 10, Γ = 1: A* = 0.16, a = 1.47ξ, κ = 48, Q̃ = 0.15. ν = 100, Γ = 0.01: A* = 0.095, a = 0.58ξ, κ = 6.9, Q̃ = 0.022.
- Predicted α = ΓQ̃²/(4πκ) at the minimum: at most ≈ 3×10⁻⁴ (near ν ≈ 3, Γ ≈ 10), falling for larger ν; about 25× below 1/137 over the scanned range. The one-population q = 1 ring under-produces charge relative to its mass and spin.
- So the choke is not the mechanism here (it only bounds small-ν cases); energy balance at fixed spin is. The α shortfall is the sharp test: next try the counter-streaming pair (different charge/moment/spin ratios) and other closures (q, branch) to see whether α can reach 1/137.

## 26. Regime with a single stable ring size
- Map over (ν, Γ) of the energy at fixed spin ħ/2 along the compressible q = 1 family (script `size_map.py`; integrals e_w/A², e_f/A², j/A fitted with cubics to remove grid kinks, which otherwise create spurious extra minima; fit rms ≈ 10% for e_w).
- Single stable size with the closed band intact (one minimum along the one constrained degree of freedom, so stable to breathing at fixed spin) over ~60% of the scanned plane: ν ≳ 4–5 for any Γ ≲ 30; the allowed Γ widens upward as ν falls (upper edge roughly Γ ~ 10³ at ν ~ 3, ~30 at ν ~ 10³).
- Edges: too stiff a field (large Γ for given ν) → energy falls toward A → 0 (weak, large ring, no minimum). Too little weighted inertia (ν ≲ 4 at small Γ) → the minimum moves to large amplitude, past the loss of the closed band (A ≈ 0.27) and into the choke.
- Size where unique: a/ξ ≈ 0.6–1.2 at ν ~ 100; 1.5–4.5 at ν ~ 10; up to ~70 at small ν with large Γ.
- Not yet checked in this regime: shape modes and tilt for the compressible solution; two-population (counter-streaming) version; relativistic inertia.

## 27. Reformulation: one population with a w-dependent drift (U, W modes)
- There is one population of fast spinners, distributed in the weight w = (|s|² − s̄²)p(χ); vacuum moments ⟨w⟩ = 0 (net 𝒦 = 0), ⟨w²⟩ = σ_w², ⟨w³⟩ = 0 (zero skewness, which is the kernel's "no energy terms odd in δ𝒦", i.e. C symmetry). The "two populations" are the sign bins of w. What matters is that the drift depends on w (the kernel conserves the weighted flux separately from the mass flux).
- Linear drift: V(w) = U + w̃W, w̃ = w/σ_w. U is the common drift: it carries the weighted mass flux (inertia νρ_m), impulse and angular momentum, and no net 𝒦 flux. W is the correlation flux: F = σ_wW, so it carries A, B and the moment, and no mass.
- Density response (isothermal, per class): ln(f/f₀) = a(x) + w̃b(x) + w̃²c(x). Net spin pressure δ𝒦 = σ_w b; weighted mass density δρ_w ∝ a + c.
- Steady Euler, order by order in w̃:
  - w̃⁰: (U·∇)U + c²∇a = 0 (mass mode)
  - w̃¹: (U·∇)W + (W·∇)U + c²∇b = C_F/σ_w (charge mode; holding cancels c²∇b outside, so b is harmonic there)
  - w̃²: (W·∇)W + c²∇c = 0 (width mode)
  Continuity ∇·(fV) = 0 order by order.
- Charge: from the w̃¹ equation, ∇b = −∇(U·W)/c² + [U×(∇×W) + W×(∇×U)]/c²; for aligned (Beltrami-type) flows b = −U·W/c², so the charge density is the overlap of the common drift and the correlation flux: δ𝒦 ∝ −σ_wU·W. A charge needs both modes together; the antiparticle has W reversed relative to U.
- Constraints decouple:
  - zero exterior mass flux → U confined to the core (closed U-flow); then impulse = 0 automatically (no exterior U dipole);
  - spin J = ħ/2 from the interior toroidal U (angular momentum νρ_m∫ r U_φ);
  - moment from the exterior dipole of W; exterior field = W as a pure correlation flux;
  - charge from the interior overlap U·W.
  The §16 tie between moment and impulse (same exterior coefficient) disappears, since U and W have independent shapes. Closure (q = N_H/N_θ) applies to the F-lines, i.e. W.
- Revised catalog (supersedes the two-population items of §5): conservation (𝒦 and its flux, mass, momentum, energy, J, P, frozen labels); force balance (the three Euler orders, holding outside, even phoretic coupling, locking); matching (U confined, W continuous with harmonic exterior); global (P = 0 automatic, J = ħ/2, charge and moment as outputs); topology (closure of W-lines); energy (minimize at fixed J, closure, number); stability (breathing, tilt neutral, shape modes open). Medium statistics enter through ν and σ_w (with zero skewness).
- Line vortex: the straight (helical) line vortex is the natural exact testbed: the three Euler orders reduce to radial ODEs, the charge per unit length is λ_q ∝ ∫U·W 2πr dr, the field energy is logarithmic, and closure becomes pitch quantization. It is also the inner limit of a thin ring (matched asymptotics).

## 28. Exact line vortex: Rankine U with Gold–Hoyle W
- Straight line, fields depending on r only: radial components vanish by continuity; the axial and azimuthal Euler components hold identically, so any profiles U_φ, U_z, W_φ, W_z are exact steady solutions. Class densities f = f₀ exp[a + w̃b + w̃²c] with a′ = U_φ²/(c²r), b′ = 2U_φW_φ/(c²r), c′ = W_φ²/(c²r) (exact for the isothermal fluid). Selection: U confined, W current-free outside, holding (b harmonic outside), closure (pitch), energy at fixed spin.
- Line charge (Gauss with holding): λ_q = −(4πε*g_Φσ_w/c²) U_φW_φ|_{r_c} — the paper's straight-tube formula in U–W form. A smooth U_φ → 0 at the edge gives a polarized neutral tube; net charge needs U_φ finite at r_c (sheet where holding switches on).
- Families: Rankine U (U_φ = Ωr inside r_c); annular free vortex U (hollow core); Lundquist W (J₀, J₁; λr_c = 2.405 zeroes W_z at the edge); Gold–Hoyle W (uniform twist); aligned U ∥ W; helical Beltrami modes J_m(kr)e^{i(mφ+kz)}.
- Rankine–Gold–Hoyle line (W_φ = W₀x/(1+x²), W_z = W₀/(1+x²), x = r/r₀), closed forms: a = Ω²r²/2, b = ΩW₀r₀ ln(1+x²), c = W₀²x²/(2(1+x²)); W-line pitch 2πr₀ uniform (closure is one condition, not one per surface); λ_q = −4πsΩW₀r_c x_c/(1+x_c²); spin per length ℓ = νπΩr_c⁴/2; E_U = ℓ²/(νπr_c⁴); magnetic energy per length πΓs²W₀², independent of r₀ and r_c; W kinetic (ν/2)πr₀²W₀² ln(1+R²/r₀²).
- Stable core radius: at fixed spin per length the U-energy falls as r_c⁻⁴ and the field energies barely depend on r_c; the bulk confinement (p̃₀πr_c² + 2πσ̃r_c, the phoretic vacuum pressure) rises. One minimum, E″ > 0, in every case scanned; leading order r_c ≈ (2ℓ²/(νπ²p̃₀))^{1/6}. Examples: ℓ = 1, ν = 10, Γ = 1 → r_c = 0.47ξ; ν = 100 → 0.29ξ; ℓ = 10 → 1.08ξ. Without the vacuum pressure there is no stable radius (confinement is phoretic, as the kernel requires).
- Still to fix: W₀ (needs a flux or helicity quantum for the W-lines) and r₀ (closure: pitch vs period, r₀ = L/(2πN), or R/N on a thin ring); spin per length is not quantized for an infinite line (period L or ring circumference supplies it). Gold–Hoyle W carries a weak exterior current (∝ r⁻⁴). Script: `line_rgh.py`.

## 29. Thin ring: tension against backreaction
- Bent into a ring, the line's spin per length is not conserved. Conserved: the Kelvin circulation Γ_U of the U-flow around the core; J_z about the ring axis (carried by the line's axial U, now toroidal); the W-line flux Φ_W; the closure integer N (twist length r₀ = R/N).
- Energy per ring (thin-ring, leading order; electric term omitted), L = 2πR:
  - shortening (∝ R): Rankine core at fixed circulation, νΓ_U²/(16π) per length (independent of r_c); vacuum confinement p̃₀πr_c² + 2πσ̃r_c per length;
  - backreaction (rises as R falls): toroidal flow holding J, J²/(8π³νr_c²R⁴) per length (∝ R⁻³ per ring); W kinetic at fixed flux, νΦ_W²N²/(2πR² ln(1+N²)) per length (∝ R⁻¹); W magnetic, πΓs²W₀² with W₀ = Φ_WN²/(πR² ln(1+N²)) (∝ R⁻³).
  - The core radius balances the toroidal-flow term (∝ r_c⁻²) against the confinement (∝ r_c²); the ring radius balances tension against backreaction.
- With self-consistent spin J = ħ/2 = κ/2 (units ρ_mcξ⁴), one minimum in R in every case scanned (script `thin_ring.py`):

| Γ_U | Φ_W | N | ν | Γ | R/ξ | r_c/ξ | r_c/R | κ |
|---|---|---|---|---|---|---|---|---|
| 2π | 1 | 1 | 10 | 1 | 1.13 | 0.57 | 0.50* | 120 |
| 2π | 1 | 1 | 10 | 100 | 2.40 | 0.33 | 0.14 | 205 |
| 2π | 1 | 2 | 10 | 1 | 1.27 | 0.56 | 0.44 | 140 |
| 2π | 3 | 1 | 10 | 1 | 2.07 | 0.42 | 0.20 | 228 |
| 2π | 1 | 1 | 100 | 1 | 1.12 | 0.56 | 0.50* | 806 |
| π | 1 | 1 | 10 | 1 | 1.66 | 0.22 | 0.13 | 53 |

  (* r_c pinned at the fat-ring limit R/2: thin-ring formulas not valid there.) At the minimum the shortening energy is ~2–3× the backreaction energy (virial balance of terms ∝ R and ∝ R⁻¹, R⁻³).
- Open: the values of the quanta Γ_U and Φ_W (a weighted circulation quantum, and a W-flux quantum tied to the internal phase); electric energy and curvature corrections (ln(8R/r_c)); thin-ring validity needs r_c/R ≪ 1 (larger Γ or Φ_W, smaller Γ_U).

## 30. The two quanta, and quark rings
- Circulation quantum. A confined Rankine core has zero circulation outside (the edge sheet cancels it), so Γ_U is conserved but not quantized there. With a Madelung phase for the weighted mass mode (same Korteweg stress, ζ_w = ζ assumed), circulation is quantized only around a multiply connected irrotational region: the annular free-vortex core, U_φ = Γ_U/(2πr) on r_i < r < r_c, hollow inside, sheet at r_c. Then Γ_U = n·2πζ, and the hollow radius is where the drift reaches c (the choke): r_i = Γ_U/(2πc) = nξ. The choke reappears as the size of the hollow core, as in a superfluid vortex (healing length). Core energy per length νπn² ln(r_c/nξ) (units ρ_mc²ξ²).
- Flux quantum for the W-lines: candidates (a) single-valuedness of the internal clock phase along closed W-lines (closure q = N_H/N_θ already uses this; a flux unit needs one more condition); (b) an Aharonov–Bohm condition for the core's own charge encircling the hole flux (q∮A = 2πħn), which brings back the 1/α scaling that made the ring too large (§11). Open.
- Quark hypothesis (circulation shared, charge = circulation × correlation flux). Charge ∝ overlap U·W, so a ring carrying a fraction f of the circulation quantum and m units of W-flux has charge f·m (units of e):
  - electron: f = 1, m = −1 → −1;
  - quarks: f = 1/3; u: m = +2 → +2/3; d: m = −1 → −1/3;
  - proton uud: circulation 1/3 + 1/3 + 1/3 = 1 quantum, charge 2/3 + 2/3 − 1/3 = +1; neutron udd: 1 quantum, charge 0;
  - mesons q q̄: circulation +1/3 − 1/3 = 0, charge e.g. u d̄: 2/3 + 1/3 = +1 (π⁺), u ū: 0 (π⁰).
- Confinement: the total U-circulation around any closed core loop must be an integer multiple of 2πζ. A lone quark ring (1/3 quantum) has a phase that is not single-valued, so it has no steady state on its own; three rings around one shared core (the Borromean triple of the brief) or a ring–antiring pair (net 0) are allowed. The mismatch for separated quarks would be carried by a wall or string whose energy grows with separation (linear confinement), as for half-quantum vortices bound by solitonic walls in multicomponent superfluids (Volovik).
- Tests: integer charges for all free states follow automatically; the W-flux unit must make the electron's charge e (α as the check); the triple must have a stable size with the shared quantum supplying the backreaction that a lone quark ring lacks.

## 31. Beltrami route: the U question and two obstructions
- The paper settles the U question: the weighted flow has its own inertia, relativistic (e+p)/c², an energy excess, not medium mass; the bulk passes through freely, with no drag, wind or entrainment; a ring of weighted circulation self-propels unless its impulse vanishes (§8 of the draft; held changes C, I). So U is a separate weighted mode, not the bulk velocity, and persistence requires the kernel not to shear-dissipate it (new kernel requirement: conservation of the weighted momentum to the persistence level). Held change I already flags "how strongly collisions couple the two drifts" as open.
- Beltrami with a common λ (∇×U = λU, ∇×W = λW): all Lamb and cross terms vanish, Bernoulli is constant in space, and a = 𝔅₀ − U²/2c², b = −U·W/c² (exact pointwise charge density −σ_wU·W), c = −W²/2c².
- Obstruction 1 (closure): a pure Beltrami W in the sphere (b = 0, surface current allowed) has q ≈ 0.80–0.825 only (λa = 3.0–4.4); q = 1 is impossible. The simplest rational in range is q = 4/5 (near λa ≈ 4.3): lines close after 5 poloidal and 4 toroidal turns, a (5,4) torus-knot electron candidate.
- Obstruction 2 (impulse–moment tie returns): with a common λ and regular ℓ = 1 fields, U and W share the radial function j₁(λρ), so U ∝ W in shape; impulse (exterior U dipole) and moment (exterior W dipole) are tied again. A confined U shell with χ_U = 0 at r_i and a needs a second zero of j₁y₁ − y₁j₁, which does not exist for λa = π or 4.49 (needs λa ≳ 7.7: nested counter-rotating shells). Different λ for U and W leaves a cross term (λ_W − λ_U)U×W that is not a gradient, so no steady axisymmetric state without extra forcing.
- Ways forward: (i) multi-shell Beltrami (λa ≳ 7.7, U in an outer or inner shell with its own zero-flux boundaries); (ii) q = 4/5 Beltrami W with U on a different spherical harmonic; (iii) non-Beltrami U (general Bragg–Hawthorne), keeping W Beltrami.

## 32. Viscosities and dissipation channels
Persistence: atoms have existed for ~4×10¹⁷ s, so any decay rate of a particle's internal structure must be ≲ 2×10⁻¹⁸ s⁻¹ (electron charge-violating decay bounds go much further, ~10⁻³⁶ s⁻¹, but those test charge conservation, which 𝒦 conservation makes exact). Healing length ξ = 3.9×10⁻¹³ m; collision rate ν_c ~ 10⁴⁰ s⁻¹; mean free path ℓ ~ 10⁻³² m.

| Channel | Symbol | Value or bound | Effect on a ring | Status |
|---|---|---|---|---|
| Bulk shear viscosity | η/ρ_m ~ cℓ | ~3×10⁻²⁴ m²/s; diffusion across ξ in ξ²/cℓ ≈ 0.05 s | destroys any bulk vortex or Rankine sheet in the particle | ordinary; tolerable only because the n = 0 ring has no bulk vorticity (irrotational added-mass flow does not shear-dissipate) |
| Weighted mass drift U: shear viscosity | ν_U | needs ν_Uλ² ≲ 2×10⁻¹⁸ s⁻¹ → ν_U ≲ 3×10⁻⁴³ m²/s, ~19 orders below cℓ | the spin-carrying flow would unwind | new kernel requirement (weighted momentum conserved, no shear relaxation) |
| U–bulk drag (mutual friction) | ν_Ub | ≲ 2×10⁻¹⁸ s⁻¹, vs ν_c ~ 10⁴⁰ s⁻¹ (58 orders) | locks U to the bulk; U would become bulk vorticity | the held-change-I kernel question; must be effectively zero |
| Correlation flux / transverse spin flux | ν_W | ≲ 3×10⁻¹⁹ s⁻¹ (cosmic transparency) | the field and moment structure would relax | kernel requirement (near conservation of W⁽²⁾); transparency and persistence give similar bounds |
| Relaxation of the 𝒦 current | part of C_F ∝ −F | must vanish; only the holding part ∝ ∇𝒦 allowed | a static Coulomb tail would flow away | kernel requirement |
| 𝒦 production | C_𝒦 | 0 | charge non-conservation | kernel requirement (Lorenz condition) |
| Locking | κ₂ | fast (≫ dynamical rates) | not dissipative: enforces B = S⁽²⁾/I₂ | assumed |
| Axial-spin relaxation | ν∥ = εν_c | ε ≲ 10⁻⁹⁰ | gravity sector; negligible for the ring | bounded |
| Inclination, spin relaxation | ν_H, ν_J, ν_s | damped modes | inclination (neutrino-analogue) and spin waves | catalogued |
| Quantum (Korteweg) stress | K_H | conservative | sets ξ, ζ | not dissipative |

- Conclusion: the soliton must live entirely in sectors the kernel conserves almost exactly (W, 𝒦 and, newly, U). Ordinary bulk viscosity is large, which is another reason the particle cannot rely on bulk vorticity: the n = 0 picture is forced, not chosen.
- The U requirements (no shear viscosity, no drag on the bulk) are the strongest new demands on the kernel and belong in its requirements table.

## 33. Ideal spin-pressure fluid and the Hall-like vortex: first isolated solutions
- Ideal multi-fluid (Euler-like): weight classes k (weights w_k, densities n_k, drifts V_k), each inviscid with relativistic inertia h_k = (e+p)/c², isothermal pressure, force n_kw_k(E + V_k×B) (class-level minimal coupling, assumed) plus holding outside cores and even phoretic coupling; fields from locking (Φ = g_Φ𝒦, A = g_ΦF/c², 𝒦 = Σn_kw_k, F = Σn_kw_kV_k). Hamiltonian H = Σ∫[½h_kn_kV_k² + n_ku_k] + ½ε*∫(E² + c²B²) + E_hold + E_ph. Conserved: class numbers, 𝒦, momentum, energy, J, canonical circulation and helicity of each class (P_k = h_kV_k + w_kA). Kernel needs: no drag, no shear viscosity, the class coupling, fast locking.
- Hall-like vortex (two classes ±w): canonical Beltrami ∇×P± = λ±(U ± W) gives (∇× − μ₁)(∇× − μ₂)W = 0, W = W₁ + W₂, U = α₁W₁ + α₂W₂, α_i = (2Λμ_i − λ₊ − λ₋)/(λ₊ − λ₋), Λ = h + wκ (κ = 2g_Φnw/c²), μ² − [(λ₊+λ₋)(h+Λ)/(2hΛ)]μ + λ₊λ₋/(hΛ) = 0. Generalized Bernoulli constant everywhere per class (exactly steady). The second scale is a weighted skin depth ~ξ√(Γ/ν).
- Simplest version (sphere, ℓ = 1, a = 1, inertia ratio η = Λ/h): conditions (1) smooth W matching to an exterior dipole, c₁μ₁²j₀(μ₁) + c₂μ₂²j₀(μ₂) = 0; (2) zero impulse, α₁c₁j₁(μ₁) + α₂c₂j₁(μ₂) = 0; (3) q_W = 1. Isolated solutions (script `double_beltrami.py`, `db_points.py`):

| η | μ₁a | μ₂a | c₂/c₁ | α₁ | α₂ | r_O/a | W exterior dipole |
|---|---|---|---|---|---|---|---|
| 1.5 | 3.850 | 8.814 | +0.496 | −0.70 | +2.15 | 0.874 | +0.203 |
| 1.5 | 4.072 | 8.238 | +0.428 | −0.62 | +2.42 | 0.524 | +0.124 |
| 1.5 | 4.722 | 13.454 | +0.453 | −0.79 | +1.91 | 0.647 | −0.066 |
| 3.0 | 3.406 | 12.188 | −0.197 | −0.42 | +7.10 | 0.876 | +0.277 |

  For η = 10 the zero-impulse curves have q = 0.73–0.79 only (no q = 1 solution in the scan). Both single-Beltrami obstructions are lifted: U and W differ in shape (α₁ ≠ α₂), so zero impulse and a nonzero moment coexist, and q is no longer pinned at 0.80–0.825.
- The solutions are discrete in shape space for given η: the isolated-solution structure needed for electron identity. U forms nested counter-rotating cells (its toroidal part changes sign), so the net spin must be checked; W can have two O-points (two families of tori). Plot: `db_plot.py`.
- Next: spin (J ≠ 0 with the counter-rotating U cells), charge (Bernoulli with Φ), moment and g for each solution; then the size a from energy at fixed spin; and which η the medium has.

## 34. Spin, charge, moment and g of the Hall-like rings
- Charge from the per-class Bernoulli relation with the self-consistent Φ = g_Φ𝒦 inside the core: 𝒦 = −2wn U·W/η (the field coupling screens the overlap charge by 1/η, since η − 1 = 2g_Φnw²/(hc²)); holding outside, Q from the boundary monopole. Spin J = νAa⁴j_s with j_s = ∫ I_U d³x (the counter-rotating U cells do not cancel). Moment from the W exterior dipole D_s. Mass = weighted kinetic energy ½νA²a³K_s + bubble confinement; J = ħ/2 = κ/2 self-consistently. Then μ/μ_B = 8πD_s a/(Q_s A) and g = 2μ/μ_B = 16πD_s a/(Q_s A). Script: `db_props.py`.

| η | μ₁a | μ₂a | Q_s | j_s | D_s | orientation μ·J·Q | g (ν = 0.1 → 100) | a/ξ, A (ν = 10) |
|---|---|---|---|---|---|---|---|---|
| 1.5 | 3.850 | 8.814 | +0.354 | −1.96 | +1.03 | wrong | 1745 → 2046 | 2.09, 0.16 |
| 1.5 | 4.072 | 8.238 | +3.45 | −2.19 | +0.61 | wrong | 73 → 85 | 1.74, 0.19 |
| 1.5 | 4.722 | 13.454 | +3.87 | −1.69 | −0.48 | right | 59 → 69 | 2.01, 0.19 |
| 3.0 | 3.406 | 12.188 | −2.09 | −1.21 | +1.08 | right | 547 → 639 | 3.12, 0.13 |

- Spin, charge and moment are all nonzero, the size is unique for each, and two solutions have the electron's orientation (moment antiparallel to spin for negative charge, or the C image). But g is 60–2000, not 2: none of these is the electron yet. g barely depends on ν.
- Diagnosis: g = 16π(D_s/Q_s)(a/A). Energy at fixed spin gives a/A ≈ 10 (ring ~2ξ, amplitude ~0.2); with D_s/Q_s ≈ 0.12 (best case) g = 2 needs a/A ≈ 0.3, i.e. a ring ~30× smaller or a flow ~30× faster than the energy balance gives. This is the kinematic obstruction of §§17–18 in a new form: the moment is too large for the charge at the energy-preferred size.
- Levers: more charge per unit moment (less screening, η → 1; other harmonics or more classes); the field energy, which adds mass but not spin and so raises the amplitude required for J = ħ/2 (lowering a/A); relativistic inertia; a wider scan in η for solutions with larger Q_s/D_s.

## 35. Analysis pipeline and first scorecard
- Repository `lukassaul/gyron-vortex-solitons`, package `pipeline/` (`python -m pipeline.run`): targets as radial functions for W and U; scorecard of steadiness (three U–W Euler orders, and the two-class fluid with field coupling via canonical vorticity), q, impulse, spin, charge, moment, energies, rigid-rotation fraction, size at fixed spin, g, orientation, peak Mach, α. Reports in `reports/`.
- First results (ν = 10, Γ = 1): single Beltrami and Hill rings are exactly steady (residual ~10⁻⁵) but carry impulse. The Hall-like rings are steady only in the field-coupled two-class fluid (canonical residual ≤ 3×10⁻³; U–W residuals 0.5–2.7): the class-level coupling w(E + V×B) is essential to them. Rigid-U variants are not steady (cross term not a gradient). Two Hall-like rings are supersonic at the size the energy picks (peak Mach 1.06, 1.72). g ranges 48–966; best α ≈ 6×10⁻⁴ (rigid-U spheromak, not steady). No target meets the electron criteria.

## 36. Batch over the Hall-like family (eta = 1.05 to 4)

40 isolated scanner solutions, all steady in the canonical two-class test. 13 pass the pipeline's q = 1 check, 8 of those have the right orientation, and all 8 are subsonic in some medium. Best g = 31 (eta = 1.05, nu = 30, Gamma = 0.01, a/xi = 1.2, peak Mach 0.99). Along the main branch (mu1 a ~ 4.6-4.95, mu2 a ~ 13.4-13.6), g rises with eta. alpha is below 1e-6. The obstruction g = 16 pi (D/Q)(a/A) stands: the needed a/A is ~25x smaller than the energy balance gives. Full report: batch/REPORT.md.

## 37. Transport coefficients in the pipeline
- `pipeline/transport.py` now tests every channel: bulk shear and volume viscosity, U shear viscosity, W diffusivity, spin-pressure (Fick) diffusion, heat conduction, U–bulk drag, inter-class drag, 𝒦-current relaxation, transverse spin, inclination and axial-spin relaxation; locking as a rate against the rotation. Each is evaluated to first order on the sized ring and gives a rate, a lifetime at the ordinary collisional value, the value allowed by persistence, and the suppression needed. `batch/transport_report.py` writes batch/TRANSPORT.md.
- Result for the 8 candidates: the diffusive channels (U viscosity, W diffusivity, Fick) must be suppressed by 10¹⁸–10²¹ below cℓ, the relaxation channels (U–bulk drag, inter-class drag, current relaxation, spin relaxation) by ~5×10⁵⁷ below ν_c. Bulk viscosity, heat conduction, inclination and axial-spin relaxation do not act at this order (n = 0, isothermal, no inclination structure). At ordinary collisional values the rings live 3×10⁻⁴–0.06 s (diffusive) or 10⁻⁴⁰ s (relaxation).
- Reading: an ideal ring persists only if the kernel conserves U, W and 𝒦 flux to these factors. The ideal solutions are the paper's starting point; the next step is the dissipative problem itself (time-dependent relaxation with all coefficients; structures selected by decay, or driven by the bulk flow through the ring).
