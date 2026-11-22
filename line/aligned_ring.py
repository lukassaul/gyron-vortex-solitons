"""Optimal ring closure of the aligned Beltrami line (U = eps W, W Lundquist with twist k r_c = kh).
Closure: q(r) = r W_z/(R W_phi) falls outward, so kink safety (q >= 1) plus an electron surface (q = 1) put q = 1 at
the core edge: R/r_c = J0(kh)/J1(kh). Thin-ring energy at fixed spin J = hbar/2 = kappa/2 (units rho_m c xi^4):
  J = 2 pi nu R^2 eps W0 r_c^2 Phi,   Q = 2 pi R eps W0^2 lam/eta,   toroidal current I_t = kh^2 W0 Phi,
  g = 4 pi I_t R^2 / Q,   alpha = Gam Q^2/(4 pi kappa)."""
import sys, json, numpy as np
sys.path.insert(0, ".")
from scipy.special import jv
from scipy.optimize import minimize_scalar
from pipeline.line import aligned, analyze as lanalyze, stability, P0, SIG

def shape(kh, n=4000):
    r = np.linspace(1e-4, 1, n); dA = 2 * np.pi * r * (r[1] - r[0])
    Wp, Wz = jv(1, kh * r), jv(0, kh * r); bp = 2 * Wp**2 / r          # eps W0^2/eta scale
    return dict(Phi=float(np.sum(Wz * dA)), K=float(np.sum((Wp**2 + Wz**2) * dA)), lam=float(-2 * np.pi * bp[-1]),
                KE=float(np.sum(bp**2 * dA)), Wpe=float(Wp[-1]), peak=float(np.max(np.hypot(Wp, Wz))), Rh=float(jv(0, kh) / jv(1, kh)))

def ring(kh, eps, nu, Gam, eta=1.0):
    s = shape(kh); Rh = s["Rh"]
    if Rh <= 1.0: return None
    L = np.log(8 * Rh)
    def E_of(rc, kap):
        R = Rh * rc; W0 = kap / (4 * np.pi * nu * eps * R**2 * rc**2 * s["Phi"])
        It = kh**2 * W0 * s["Phi"]; Cw = 2 * np.pi * rc * W0 * s["Wpe"]; lam = eps * W0**2 * s["lam"] / eta
        per = (0.5 * nu * W0**2 * rc**2 * (1 + eps**2) * s["K"] + 0.5 * Gam * W0**2 * kh**2 * s["K"]
               + P0 * np.pi * rc**2 + 2 * np.pi * SIG * rc + 0.5 * nu * Cw**2 / (2 * np.pi) * L
               + 0.5 * Gam * (eps**2 * W0**4 / eta**2 * s["KE"] + lam**2 / (2 * np.pi) * L))
        Eind = 0.5 * Gam * It**2 * R * max(L - 2, 0)
        return 2 * np.pi * R * per + Eind, dict(R=R, W0=W0, It=It, lam=lam)
    kap = 50.0
    for _ in range(400):
        xs = np.exp(np.linspace(np.log(1e-2), np.log(1e3), 500)); Es = np.array([E_of(x, kap)[0] for x in xs])
        mins = [i for i in range(1, len(xs) - 1) if Es[i] < Es[i - 1] and Es[i] < Es[i + 1]]
        if not mins: return None
        i = mins[0]; r = minimize_scalar(lambda lx: E_of(np.exp(lx), kap)[0], bounds=(np.log(xs[i - 1]), np.log(xs[i + 1])), method="bounded")
        rc = float(np.exp(r.x)); En, p = E_of(rc, kap)
        if abs(En - kap) < 1e-9 * kap: break
        kap = 0.5 * kap + 0.5 * En
    h = 1e-3 * rc; E2 = (E_of(rc + h, kap)[0] - 2 * En + E_of(rc - h, kap)[0]) / h**2
    Q = 2 * np.pi * p["R"] * p["lam"]; g = 4 * np.pi * p["It"] * p["R"]**2 / abs(Q)
    mach = (1 + eps) * p["W0"] * s["peak"]
    out = dict(kh=kh, eps=eps, nu=nu, Gam=Gam, eta=eta, R_over_rc=Rh, rc_xi=rc, R_xi=p["R"], W0=p["W0"], kappa=kap, Q=Q,
               g=g, alpha=Gam * Q**2 / (4 * np.pi * kap), mach=mach, n_minima=len(mins), breathing_E2=E2,
               g_estimate=2 * kh**2 * s["Phi"] * p["R"] / (eps * p["W0"] * abs(s["lam"])))
    # local stability of the line cross-section at this size
    lr = lanalyze(aligned(kh, eps), nu=nu, Gam=Gam)
    lr.update(Omega=p["W0"] / rc, W0=p["W0"], nu=nu, Gam=Gam, rc_xi=rc, breathing_E2=E2, mach_peak=mach)
    st = stability(lr)
    passed = all([st["S1_rayleigh+"], st["S1_rayleigh-"], st["S2_michael+"], st["S2_michael-"],
                  min(st["S3_HG_Ri_min+"], st["S3_HG_Ri_min-"]) >= 0.25, st["S6_KH"] in ("stable", "no sheet"), E2 > 0, mach < 1])
    out.update(stab={k: v for k, v in st.items()}, all_local_pass=bool(passed))
    return out

if __name__ == "__main__":
    rows = []
    for kh in np.round(np.linspace(0.2, 2.3, 22), 3):
        for eps in (0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0):
            for nu in (1, 10, 100):
                for G in (0.01, 1, 100):
                    o = ring(kh, eps, nu, G)
                    if o: rows.append(o)
    json.dump(rows, open("line/aligned_ring.json", "w"), default=float, indent=0)
    ok = [o for o in rows if o["all_local_pass"]]
    print(len(rows), "rings sized;", len(ok), "pass all local tests and subsonic")
    for o in sorted(ok, key=lambda o: abs(np.log(o["g"] / 2)))[:15]:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in o.items() if k not in ("stab",)})
    # trend: g vs kh for passing rings
    import collections; d = collections.defaultdict(list)
    for o in ok: d[o["kh"]].append(o["g"])
    print("kh: min g among passing:", {k: round(min(v), 2) for k, v in sorted(d.items())})
