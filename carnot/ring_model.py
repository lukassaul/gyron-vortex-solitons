"""Explicit Carnot ring with Barnett-locked weighted spin (route 2). Thin ring, units xi, c, rho_m.
Cross-section (radius w around the core centre, ring radius R):
  hollow core   w < r_i = n            (choke: |U| = c at the free surface)
  free vortex   r_i < w < a            U_th = n/w               (quantized circulation 2 pi n zeta)
  return layer  w > a                  U_th = (n/w) exp(-((w-a)/d)^2)
  spin flow     U_s(w) = U0 sin^2(pi (w - r_i)/(w_out - r_i)) on r_i < w < w_out = a + 3d (toroidal, carries J)
  weighted density (isothermal Bernoulli) n_w = exp(-(U_th^2 + U_s^2)/2), zero in the hollow core
  locked intrinsic spin s = chi(w) omega/2, chi = chi0 n_w; weighted spin S2 = p_w s (moment), J_int = -pi^2 int chi' w^2 U_s dw
Energy E = 2 pi R [ nu int (U_th^2 + U_s^2)/2 n_w dA + int chi omega^2/8 dA + p0 pi w_out^2 + 2 pi sig w_out ].
Dissipation D = 2 pi R nu nu_v int n_w [ (w (U_th/w)')^2 + U_s'^2 ] dA.  Engine W = C E.
Fast: R from dE/dR = 0 at fixed J.  Self-consistent spin: J_orb + J_int = kappa/2, kappa = E.  Slow: core edge a by E'(a) a' = W - D."""
import numpy as np, json
from scipy.optimize import brentq, minimize_scalar
P0, SIG = 5 / 7, 2 * np.sqrt(2) / 3

class Ring:
    def __init__(self, n=1, nu=100.0, nu_v=1e-3, chi0=1e-6, p_w=0.5, dfrac=0.5, N=1500):
        self.n, self.nu, self.nu_v, self.chi0, self.p_w, self.dfrac, self.N = n, nu, nu_v, chi0, p_w, dfrac, N
    def profiles(self, a, U0):
        n = self.n; ri = n; d = self.dfrac * a; wo = a + 3 * d
        w = np.linspace(ri, wo, self.N); Uth = np.where(w <= a, n / w, n / w * np.exp(-((w - a) / d)**2))
        Us = U0 * np.sin(np.pi * (w - ri) / (wo - ri))**2
        for _ in range(3): pass
        nw = np.exp(-(Uth**2 + Us**2) / 2)
        return w, Uth, Us, nw, wo
    def per_length(self, a, U0):
        w, Uth, Us, nw, wo = self.profiles(a, U0); dA = 2 * np.pi * w * np.gradient(w)
        I = lambda f: float(np.sum(f * dA))
        om_s = np.gradient(w * Uth, w) / w; om_t = -np.gradient(Us, w); chi = self.chi0 * nw
        K = self.nu * I(0.5 * (Uth**2 + Us**2) * nw); Esp = I(chi * (om_s**2 + om_t**2) / 8)
        Econf = P0 * np.pi * wo**2 + 2 * np.pi * SIG * wo
        Dl = self.nu * self.nu_v * I(nw * ((w * np.gradient(Uth / w, w))**2 + np.gradient(Us, w)**2))
        L = self.nu * I(nw * Us)                                        # toroidal momentum per length
        jint = -np.pi**2 * float(np.sum(np.gradient(chi, w) * w**2 * Us * np.gradient(w)))   # J_int / (2 pi R)?  see below
        return dict(K=K, Esp=Esp, Econf=Econf, D=Dl, L=L, jint=jint)
    def at(self, a, kappa):
        """Given core edge a and mass kappa (J = kappa/2): solve U0 and R from J and dE/dR = 0."""
        J = kappa / 2
        def resid(U0):
            p = self.per_length(a, U0)
            # E(R) = 2 pi R (K + Esp + Econf); with J_orb = 2 pi R^2 L fixed: L ~ U0 -> K_tor ~ U0^2 ~ 1/R^4
            R = np.sqrt(max(J - p["jint"], 1e-12) / (2 * np.pi * max(p["L"], 1e-30)))
            return R, p
        # mechanical equilibrium: minimize E over U0 (equivalently R) at fixed J
        def Eof(lu):
            R, p = resid(np.exp(lu)); return 2 * np.pi * R * (p["K"] + p["Esp"] + p["Econf"])
        r = minimize_scalar(Eof, bounds=(np.log(1e-4), np.log(0.95)), method="bounded")
        U0 = float(np.exp(r.x)); R, p = resid(U0)
        return dict(a=a, U0=U0, R=R, E=Eof(r.x), D=2 * np.pi * R * p["D"], J_orb=2 * np.pi * R**2 * p["L"], J_int=p["jint"],
                    Esp_frac=p["Esp"] / (p["K"] + p["Esp"] + p["Econf"]), wout=a + 3 * self.dfrac * a)
    def solve(self, C, a_grid=None):
        a_grid = a_grid if a_grid is not None else np.geomspace(self.n * 1.02, self.n * 30, 120)
        kappa = 200.0
        for it in range(60):
            rows = [self.at(a, kappa) for a in a_grid]; E = np.array([r["E"] for r in rows]); D = np.array([r["D"] for r in rows])
            f = C * E - D; idx = np.nonzero(np.diff(np.sign(f)))[0]
            if not len(idx): return None, rows
            i = idx[0]; x = a_grid[i] - f[i] * (a_grid[i + 1] - a_grid[i]) / (f[i + 1] - f[i]); s = self.at(x, kappa)
            if abs(s["E"] - kappa) < 1e-6 * kappa: break
            kappa = 0.5 * kappa + 0.5 * s["E"]
        dE = np.gradient(E, a_grid); df = np.gradient(f, a_grid)
        s.update(kappa=kappa, stable=bool(np.interp(x, a_grid, df) / np.interp(x, a_grid, dE) < 0),
                 aE=float(a_grid[np.argmin(E)]), mach_spin=s["U0"], C=C)
        s["moment_integral"] = self.p_w * s["J_int"]          # int S2 dV (z): m = (3/2) g_e^2 * this
        return s, rows
if __name__ == "__main__":
    out = {}
    for C in (0.003, 0.01, 0.03, 0.1):
        rg = Ring(n=1, nu=100.0, nu_v=1e-3, chi0=1e-6); s, _ = rg.solve(C)
        if s is None: print(f"C={C}: no balance"); continue
        out[C] = s
        print(f"C={C:<6} a*={s['a']:.3f} (choke 1, a_E={s['aE']:.2f}) R*={s['R']:.2f} R/a={s['R']/s['a']:.1f} kappa={s['kappa']:.1f} "
              f"U0={s['U0']:.3f} J_orb={s['J_orb']:.2f} J_int={s['J_int']:+.2e} stable={s['stable']} spin-energy frac={s['Esp_frac']:.1e}")
    json.dump({str(k): v for k, v in out.items()}, open("carnot/ring_model.json", "w"), default=float)
