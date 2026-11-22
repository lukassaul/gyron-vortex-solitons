"""Atlas page 1: driven (Carnot) vortex rings in a two-class spin fluid, in dimensionless groups.
Units: length r_i = Gamma/(2 pi c) (the choke radius: free-vortex speed reaches the sound speed c), velocity c,
density = far-field weighted density. Winding enters only through r_i (a ring of winding n is the n = 1 ring scaled).
Cross-section (radius w, ring radius R >> w): hollow core w < 1; free vortex 1 < w < a, U_th = 1/w; return layer
U_th = exp(-((w-a)/d)^2)/w (d = a/2) out to w_out = a + 3d; toroidal spin flow U_s = U0 sin^2(pi (w-1)/(w_out-1));
weighted density n_w = exp(-(U_th^2 + U_s^2)/2) (isothermal Bernoulli).
Groups: P = vacuum confinement pressure / (rho c^2); S = surface tension / (rho c^2 r_i); J~ = J/(rho c r_i^4);
Carnot number N = C r_i^2 / nu_v (engine rate C over viscous rate); Re = Gamma/nu_v.
Ring: E = 2 pi R [K + P pi w_out^2 + 2 pi S w_out], D = 2 pi R nu_v Dhat, W = C E, J = 2 pi R^2 L. R in mechanical
equilibrium at fixed J (fast). Slow budget E'(a) a' = W - D at fixed J and C.
Conditions: balance (N = Dhat/E per unit), Carnot-stable ((W-D)'/E' < 0), subsonic (max |U| < 1 outside the hollow core),
thin (R/w_out >= 2), and edge stability vs Re: centrifugal (Synge, return layer) and spin-jet shear, damped viscously."""
import json, itertools, numpy as np
from scipy.optimize import minimize_scalar
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
NW = 500
F = 0.0   # exterior circulation fraction (0 = confined, 1 = open hollow ring)
def section(a, U0, dfrac=0.5):
    d = dfrac * a; wo = a + 3 * d; w = np.linspace(1.0, wo, NW); dw = w[1] - w[0]
    Uth = np.where(w <= a, 1 / w, (F + (1 - F) * np.exp(-((w - a) / d)**2)) / w); Us = U0 * np.sin(np.pi * (w - 1) / (wo - 1))**2
    nw = np.exp(-(Uth**2 + Us**2) / 2); dA = 2 * np.pi * w * dw
    return w, Uth, Us, nw, dA, wo, d
def ring_at(a, J, P, S):
    """Mechanical equilibrium at fixed J: returns dict or None."""
    def parts(U0):
        w, Uth, Us, nw, dA, wo, d = section(a, U0)
        K = np.sum(0.5 * (Uth**2 + Us**2) * nw * dA); L = np.sum(nw * Us * dA)
        Ec = P * np.pi * wo**2 + 2 * np.pi * S * wo
        Dh = np.sum(nw * ((w * np.gradient(Uth / w, w))**2 + np.gradient(Us, w)**2) * dA)
        R = np.sqrt(J / (2 * np.pi * L)); return R, K, Ec, Dh, wo
    G = 2 * np.pi * F
    Eext = lambda R, wo: 0.5 * G**2 * R * max(np.log(8 * R / wo) - 2, 0.0)
    Dext = lambda R, wo: 2 * np.pi * R * G**2 / (np.pi * wo**2)
    def E(lu):
        R, K, Ec, Dh, wo = parts(np.exp(lu)); return 2 * np.pi * R * (K + Ec) + Eext(R, wo)
    r = minimize_scalar(E, bounds=(np.log(1e-4), np.log(0.98)), method="bounded", options=dict(xatol=1e-6))
    U0 = float(np.exp(r.x)); R, K, Ec, Dh, wo = parts(U0)
    Uself = (G / (4 * np.pi * R)) * (np.log(8 * R / wo) - 0.5) if F > 0 else 0.0
    return dict(a=a, U0=U0, R=R, E=2 * np.pi * R * (K + Ec) + Eext(R, wo), Dh=2 * np.pi * R * Dh + Dext(R, wo), wo=wo, edgeU0=U0 > 0.97,
                impulse=float(np.pi * R**2 * G), Uself=float(Uself))
def edge_stability(a, U0):
    """Critical Reynolds numbers above which the edge layers go unstable (estimates)."""
    w, Uth, Us, nw, dA, wo, d = section(a, U0)
    Phi = np.gradient(nw * (w * Uth)**2, w) / (nw * w**3)                  # Synge (variable density) discriminant
    sig_c = np.sqrt(max(-Phi.min(), 0)); Re_c1 = np.inf if sig_c == 0 else 2 * np.pi * np.pi**2 / (sig_c * d**2)
    sh = np.max(np.abs(np.gradient(Us, w))); wj = wo - 1
    Re_c2 = np.inf if sh == 0 else 2 * np.pi * (2 * np.pi / wj)**2 / (0.2 * sh)
    return Re_c1, Re_c2
def page(P, S, J, agrid):
    rows = [ring_at(a, J, P, S) for a in agrid]
    E = np.array([r["E"] for r in rows]); Dh = np.array([r["Dh"] for r in rows]); N = Dh / E
    dE = np.gradient(E, agrid); dD = np.gradient(Dh, agrid)
    stab = (N * dE - dD) / dE < 0                                          # (W - D)'/E' < 0 with C = N nu_v
    out = []
    for r, Ni, st in zip(rows, N, stab):
        sub = (not r["edgeU0"]) and r["U0"] < 1 / np.sqrt(2)            # combined |U| < 1 where U_th ~ 1/w <= 1
        thin = r["R"] / r["wo"] >= 2
        Re1, Re2 = edge_stability(r["a"], r["U0"])
        out.append(dict(**r, N=float(Ni), carnot_stable=bool(st), subsonic=bool(sub), thin=bool(thin), Re_c1=float(Re1), Re_c2=float(Re2),
                        ok=bool(st and sub and thin)))
    return out
if __name__ == "__main__":
    import carnot_ring as M
    agrid = np.geomspace(1.05, 60, 40)
    Fs = [0.0, 0.25, 0.5, 0.75, 0.9, 1.0]; Ps = list(np.geomspace(1e-3, 1, 5)); Srat = [0.3, 3.0]; Js = list(np.geomspace(1e2, 1e8, 7))
    atlas = []
    for f in Fs:
        M.F = f; globals()["F"] = f
        for P, sr, J in itertools.product(Ps, Srat, Js):
            rows = page(P, sr * P, J, agrid); ok = [r for r in rows if r["ok"]]
            aE = float(agrid[np.argmin([r["E"] for r in rows])])
            rec = dict(F=f, P=P, S_over_P=sr, J=J, aE=aE, n_ok=len(ok))
            if ok:
                ex = ok[len(ok) // 2]
                rec.update(a_min=min(r["a"] for r in ok), a_max=max(r["a"] for r in ok), N_min=min(r["N"] for r in ok), N_max=max(r["N"] for r in ok),
                           Re_c=max(min(r["Re_c1"], r["Re_c2"]) for r in ok), example=ex)
            atlas.append(rec)
        n = [x for x in atlas if x["F"] == f]
        print(f"F={f}: {sum(1 for x in n if x['n_ok'])}/{len(n)} points with a stable window; a_E median {np.median([x['aE'] for x in n]):.2f}", flush=True)
    json.dump(atlas, open("carnot_ring_atlas.json", "w"), default=float)
