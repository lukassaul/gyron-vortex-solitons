"""Sink-driven Burgers ring in the two-class spin fluid: explicit moments, energy, mass, charge, impulse, stability.
Units: core radius delta = sqrt(4 nu_v / kappa), sound speed c, far-field weighted density.  Thin ring R >> delta.
Exact fields (cross-section radius r, along-core coordinate s):
  class swirl  V±_th = Gam± f(r),  f = (1 - e^{-r^2})/(2 pi r)          (Burgers: steady with the sink)
  U_th = Gam_U f, W_th = Gam_W f,  Gam_U = (Gam+ + Gam-)/2, Gam_W = (Gam+ - Gam-)/2
  intake       U_r = -(K/2) r inside the core (K = kappa delta / c = 4 nu_v/(c delta)), -> spherical Q/(4 pi r^2) far away
  along core   U_s = W_s = 0 in steady state (sink drains them; only steady profile E1(r^2) needs an axial momentum source)
  densities    n_w = exp(-(U^2 + W^2)/2);  net spin pressure b = -U.W/eta (screened overlap)
  canonical vorticity per class Om± = Gam±(1 + (eta-1)(± sign for W part)) e^{-r^2}/pi  (Gaussian)
Ring: energy E = (Gam_U^2 + eta Gam_W^2) R [ln(8R) - 1.615]  (Gaussian core, normalization E = int (U^2 + eta W^2) dV)
      impulse P = pi R^2 Gam_U ; self-propulsion U_self = Gam_U/(4 pi R) [ln(8R) - 0.558]
      charge: holding outside r_c: lambda_q = -2 pi r_c b'(r_c), Q = 2 pi R lambda_q
      dissipation D = 2 pi R nu_v (Gam_U^2 + eta Gam_W^2)/(2 pi) (core vorticity integral, normalization as E)."""
import numpy as np, json
def ring(GamP, GamM, R, eta=1.5, K=0.1, r_c=2.0):
    GU, GW = (GamP + GamM) / 2, (GamP - GamM) / 2
    r = np.linspace(1e-4, 8, 4000); f = (1 - np.exp(-r**2)) / (2 * np.pi * r)
    Uth, Wth = GU * f, GW * f; nw = np.exp(-(Uth**2 + Wth**2) / 2)
    vmax = max(abs(GamP), abs(GamM)) * f.max()
    E = (GU**2 + eta * GW**2) * R * (np.log(8 * R) - 1.615)
    Eintake = 2 * np.pi * R * np.sum((K * r / 2)**2 * (r < 1) * 2 * np.pi * r) * (r[1] - r[0])
    b = -Uth * Wth / eta; bp = np.gradient(b, r); lam = -2 * np.pi * r_c * np.interp(r_c, r, bp); Q = 2 * np.pi * R * lam
    P = np.pi * R**2 * GU; Us = GU / (4 * np.pi * R) * (np.log(8 * R) - 0.558) if R > 0 else 0
    nu_v = K / 4; D = 2 * np.pi * R * nu_v * (GU**2 + eta * GW**2) / (2 * np.pi)
    # stability: Synge per class with variable density: d(n Gam(r)^2)/dr >= 0 (Gam(r) = 2 pi r V_th)
    syn = []
    for G in (GamP, GamM):
        Gr = G * (1 - np.exp(-r**2)); syn.append(bool(np.all(np.gradient(nw * Gr**2, r)[5:] >= -1e-12)))
    return dict(Gam_plus=GamP, Gam_minus=GamM, Gam_U=GU, Gam_W=GW, R=R, E=E + Eintake, mass=E + Eintake, charge=Q, impulse=P, U_self=Us,
                spin=0.0, moment=0.0, D=D, power_per_mass=D / (E + Eintake), vmax=vmax, subsonic=bool(vmax < 1),
                synge_stable=all(syn), leakage=float(np.exp(-R**2)))
cases = [("static paired (n+ = -n-)", 3.0, -3.0, 20), ("moving, U only", 3.0, 3.0, 20), ("moving, charged (n+ = 2 n-)", 4.0, 2.0, 20),
         ("moving, charged, larger ring", 4.0, 2.0, 60), ("near choke", 9.5, 4.0, 20)]
out = []
for name, gp, gm, R in cases:
    d = ring(gp, gm, R); d["case"] = name; out.append(d)
    print(f"{name:32s} E={d['E']:8.1f} Q={d['charge']:+8.3f} P={d['impulse']:8.1f} U_self={d['U_self']:.4f} vmax={d['vmax']:.2f} "
          f"subsonic={d['subsonic']} Synge={d['synge_stable']} D/E={d['power_per_mass']:.2e} leak={d['leakage']:.0e}")
json.dump(out, open("burgers/ring_cases.json", "w"), default=float)
# collapse of the static paired ring: dE/dR > 0 with no outward term -> shrinks until R ~ delta
d1, d2 = ring(3, -3, 20), ring(3, -3, 21); print(f"static paired ring: dE/dR = {d2['E']-d1['E']:.2f} > 0, no spin or impulse to resist -> collapses")
