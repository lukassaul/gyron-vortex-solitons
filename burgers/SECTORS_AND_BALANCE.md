# Vortex sources by sector, and the ring's force balance with spin pressure

## 1. The three sectors (Paper I definitions)
| Sector | Moment | Couples to |
|---|---|---|
| Spin pressure | ε_s = n⟨|s|² − s̄²⟩/2I (excess spin energy; unweighted), its flux, pressure p_s | pressure balance; heated or polarized by rotation; blind to inclination |
| Inclination | H = ⟨cos²χ⟩, flux J_H, axial label s∥ = |s|cos χ | gravity (axial-spin conversion and counterflow), Bohm/Korteweg stress K_H|∇√H|², inclination waves |
| Electromagnetism | 𝒦 = ⟨(|s|² − s̄²)p(χ)⟩ (covariance), flux F | charge and field; nonzero only where spin magnitude and inclination are correlated |

## 2. What a vortex feeds into each
**Spin pressure.**
- Polarization: spin–vorticity exchange locks the mean spin to ω/2; ordered spin energy χ_s ω²/8 per volume (χ_s = nI).
- Heating: rotational-viscosity loss 4η_r(Ω_s − ω/2)² (transients, edge layers) and the shear loss D, deposited in the core.
- Cooling: conversion moves tumbling (transverse-spin) energy into axial spin, which exerts no pressure: −q∥ = −εκ₁δe per volume.
- Steady core excess: Δp_s ≈ (heating − conversion)/ν_ex, with ν_ex the exchange rate to the translational bath.

**Inclination.** Conversion raises the axial fraction in the core: δH > 0, the gravity source, and an inclination-gradient (Bohm) stress from ∇√H; inclination waves on any change.

**𝒦 (charge).** 𝒦 ≠ 0 only if the core's change in spin magnitude correlates with its change in inclination: whether the faster (or slower) spinners in the core are more axial or more tumbling. Sign of the charge = sign of that correlation. Spin heating alone, or inclination change alone, is neutral.

**Neutrality of matter.** The counterflow (axial spinners out, tumblers in) must carry zero net 𝒦: the outgoing and incoming gyrons must have equal spin-magnitude–inclination correlation (a kernel condition), or every massive body would carry charge ∝ its conversion rate.

## 3. Net core heating and the Carnot efficiency
In the engine picture the shear loss D equals the work W, and the conversion removes Q_c = W/(T₀/T∥ − 1). If both act in the core, the net core heating is

  D − Q_c = Q_c (T₀/T∥ − 2),

so **the core is over-pressured (Δp_s > 0) iff T₀/T∥ > 2**, i.e. Carnot efficiency above ½; otherwise the conversion cools it below the vacuum (a low-pressure core). (Assumes the dissipated heat is deposited in the core's spins before exchange; conduction moves some outward.)

## 4. Ring dynamics: two cases
**Rings with mass circulation (Γ_U ≠ 0).** The core obeys Kutta–Joukowski: a force per length f on the core gives a relative velocity f × Γ̂/(ρΓ). A radial (hoop) force changes the translation speed, not R; R changes only through axial forces or impulse changes. The sink's inward pull on the core (u_R ≈ −(q/4πR)[ln(8R/s) − 1] for line sink strength q and width s) is countered only by an axial thrust f_z = ρΓ|u_R|. An isotropic exhaust gives none; a chiral or asymmetric one could. A scalar spin over-pressure cannot hold R here.

**Static paired rings (Γ_U = 0, W circulation only).** No Magnus force on the mass; ordinary force balance holds. With W line tension T_W ∝ ηΓ_W²[ln(8R/a) + const] and core over-pressure Δp_s over cross-section πa²,

  E(R) = 2πR [T_W ln(8R/a) − πa²Δp_s] + const,
  dE/dR = 0 → R* = (a/8) exp(πa²Δp_s/T_W − 1),  d²E/dR² = 2πT_W/R > 0 (stable),

which exists when πa²Δp_s ≳ T_W. The ring size depends exponentially on the over-pressure-to-tension ratio, so sizes span a wide range. With the §3 result, a static ring can be held open only by an engine with T₀/T∥ > 2.

## 5. What to compute next
1. Δp_s(r) from the heating and conversion profiles with exchange and conduction (core over-pressure vs radius).
2. The static paired ring: R* from the balance, including the sink's advection of a zero-circulation core, then a time-evolution check like §49.
3. 𝒦 from the correlated changes in spin pressure and inclination: charge sign and size.
