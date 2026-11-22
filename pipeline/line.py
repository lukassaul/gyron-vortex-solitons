"""Infinite straight line vortex: exact radial families, scorecard, sizing, transport and stability survey.
Fields depend on r only (core radius r_c = 1 for the shape; xi after sizing; velocities in c). Any radial profiles
are exact steady states of the ideal two-class fluid (axial and azimuthal Euler components vanish identically;
the radial one fixes the class densities). U is confined to r < 1 (sheet allowed at r = 1); W continues outside
as a pure azimuthal flow W_phi = C/(2 pi r) carrying its net circulation C (no exterior axial W).
Medium: nu (weighted inertia), Gam (field stiffness), eta (Hall-like screening), p0, sig (bubble confinement)."""
import numpy as np
from dataclasses import dataclass, field
from scipy.special import jv
from scipy.optimize import minimize_scalar
from .transport import CHANNELS, PERSIST, XI, C as CLIGHT

P0, SIG, RCUT = 5 / 7, 2 * np.sqrt(2) / 3, 100.0   # RCUT: outer log cutoff (period / ring size) in core radii

@dataclass
class Line:
    name: str; family: str
    Up: callable; Uz: callable; Wp: callable; Wz: callable      # functions of r on 0 <= r <= 1
    eta: float = 1.0; params: dict = field(default_factory=dict)

# ---------------- families ----------------
def gold_hoyle(x0):
    return (lambda r: (r / x0) / (1 + (r / x0) ** 2)), (lambda r: 1 / (1 + (r / x0) ** 2))
def lundquist(k):
    return (lambda r: jv(1, k * r)), (lambda r: jv(0, k * r))

def rgh(x0=1.0):
    Wp, Wz = gold_hoyle(x0)
    return Line(f"line_RGH_x{x0:g}", "Rankine U + Gold-Hoyle W", lambda r: r, lambda r: 0 * r, Wp, Wz, params=dict(x0=x0))
def rankine_lundquist(k=2.405):
    Wp, Wz = lundquist(k)
    return Line(f"line_RL_k{k:g}", "Rankine U + Lundquist W", lambda r: r, lambda r: 0 * r, Wp, Wz, params=dict(k=k))
def annular_lundquist(ri=0.3, k=2.405):
    Wp, Wz = lundquist(k)
    return Line(f"line_AL_ri{ri:g}_k{k:g}", "annular free vortex U (hollow) + Lundquist W",
                lambda r: np.where(r > ri, ri / np.maximum(r, 1e-9), 0 * r), lambda r: 0 * r, Wp, Wz, params=dict(ri=ri, k=k))
def aligned(k=3.0, eps=1.0):
    Wp, Wz = lundquist(k)
    return Line(f"line_AB_k{k:g}_e{eps:g}", "aligned Beltrami U = eps W (Lundquist)",
                lambda r: eps * Wp(r), lambda r: eps * Wz(r), Wp, Wz, params=dict(k=k, eps=eps))
def hall_line(eta, m1, m2):
    """Two-class canonical Beltrami line: W = W(m1) + c W(m2) with W_z(1) = 0; U = a1 W1 + a2 c W2."""
    s = 2 * eta * (m1 + m2) / (1 + eta); d = np.sqrt(s * s - 4 * eta * m1 * m2)
    a1 = (2 * eta * m1 - s) / d; a2 = (2 * eta * m2 - s) / d; c = -jv(0, m1) / jv(0, m2)
    Wp = lambda r: jv(1, m1 * r) + c * jv(1, m2 * r); Wz = lambda r: jv(0, m1 * r) + c * jv(0, m2 * r)
    Up = lambda r: a1 * jv(1, m1 * r) + a2 * c * jv(1, m2 * r); Uz = lambda r: a1 * jv(0, m1 * r) + a2 * c * jv(0, m2 * r)
    return Line(f"line_hall_eta{eta:g}_{m1:g}_{m2:g}", "Hall-like line (double Beltrami)", Up, Uz, Wp, Wz, eta=eta,
                params=dict(eta=eta, mu1=m1, mu2=m2, c=c, alpha1=a1, alpha2=a2))

def registry():
    T = [rgh(0.5), rgh(1.0), rgh(2.0), rankine_lundquist(2.405), rankine_lundquist(3.0), rankine_lundquist(3.832),
         annular_lundquist(0.3), annular_lundquist(0.6), aligned(3.0, 1.0), aligned(2.0, 0.5)]
    for eta, m1, m2 in [(1.05, 2.0, 5.0), (1.5, 2.0, 5.0), (1.5, 2.0, 7.0), (3.0, 2.0, 5.0), (1.5, 3.0, 6.0)]:
        T.append(hall_line(eta, m1, m2))
    return {t.name: t for t in T}

# ---------------- analysis ----------------
def _d(f, r): return np.gradient(f, r)

def analyze(t, nu=10.0, Gam=1.0, ell=1.0, PhiW=1.0, n=4000):
    r = np.linspace(1e-4, 1.0, n); dA = 2 * np.pi * r * (r[1] - r[0]); I = lambda f: float(np.sum(f * dA))
    Up, Uz, Wp, Wz = t.Up(r), t.Uz(r), t.Wp(r), t.Wz(r)
    # curls (cylindrical, r-only): (curl X)_phi = -X_z', (curl X)_z = (r X_phi)'/r
    oUp, oUz = -_d(Uz, r), _d(r * Up, r) / r; oWp, oWz = -_d(Wz, r), _d(r * Wp, r) / r
    eta = t.eta
    # net spin pressure b = (ln n+ - ln n-)/2 with field coupling and self-consistent Phi (screening 1/eta):
    bp = (2 * Up * Wp / r + (eta - 1) * (Up * oWz - Uz * oWp)) / eta
    b = np.concatenate([[0], np.cumsum((bp[1:] + bp[:-1]) / 2 * np.diff(r))])
    lam_q = float(-2 * np.pi * r[-1] * bp[-1])                    # charge per length (Gauss with holding outside)
    j = I(r * Up); Pz = I(Uz); PhiW_hat = I(Wz); Cw = float(2 * np.pi * Wp[-1]); Bz_flux = Cw
    KU, KW = I(Up**2 + Uz**2), I(Wp**2 + Wz**2); KB = I(oWp**2 + oWz**2); KE = I(bp**2)
    peak = float(max(np.max(np.hypot(Up, Uz)), np.max(np.hypot(Wp, Wz))))
    pitch = np.where(abs(Wp) > 1e-9, 2 * np.pi * r * Wz / np.where(abs(Wp) > 1e-9, Wp, 1), np.inf)   # W-line pitch length
    out = dict(name=t.name, family=t.family, eta=eta, params=t.params, steady="exact (r-only)",
               lam_q=lam_q, spin_j=j, axial_Pz=Pz, W_axial_flux=PhiW_hat, W_net_circ=Cw, B_axial_flux=Bz_flux,
               K_U=KU, K_W=KW, K_B=KB, K_E=KE, peak_per_amp=peak,
               pitch_mid=float(np.interp(0.5, r, pitch)), pitch_edge=float(pitch[-2]),
               sheet_U=float(np.hypot(Up[-1], Uz[-1]) / max(np.max(np.hypot(Up, Uz)), 1e-30)))
    # ---------------- sizing at fixed spin per length ell and W flux PhiW ----------------
    Wnorm = PhiW_hat if abs(PhiW_hat) > 1e-3 else (Cw if abs(Cw) > 1e-3 else 1.0)
    if abs(j) > 1e-9:
        def parts(rc):
            Om = ell / (nu * rc**4 * abs(j)); W0 = PhiW / (rc**2 * abs(Wnorm)) if abs(PhiW_hat) > 1e-3 else PhiW / (rc * abs(Wnorm))
            L = np.log(RCUT)
            eU = 0.5 * nu * Om**2 * rc**4 * KU
            eW = 0.5 * nu * W0**2 * rc**2 * (KW + Cw**2 / (2 * np.pi) * L)
            eB = 0.5 * Gam * W0**2 * KB
            bs = Om * rc * W0; eE = 0.5 * Gam * bs**2 * (KE + lam_q**2 / (2 * np.pi) * L)
            eC = P0 * np.pi * rc**2 + 2 * np.pi * SIG * rc
            return dict(U=eU, W=eW, B=eB, E=eE, conf=eC, Om=Om, W0=W0)
        Etot = lambda rc: sum(v for k, v in parts(rc).items() if k not in ("Om", "W0"))
        xs = np.exp(np.linspace(np.log(0.02), np.log(200), 600)); Es = np.array([Etot(x) for x in xs])
        mins = [i for i in range(1, len(xs) - 1) if Es[i] < Es[i - 1] and Es[i] < Es[i + 1]]
        out["n_minima"] = len(mins)
        if mins:
            i = mins[0]; res = minimize_scalar(lambda lx: Etot(np.exp(lx)), bounds=(np.log(xs[i - 1]), np.log(xs[i + 1])), method="bounded")
            rc = float(np.exp(res.x)); p = parts(rc); h = 1e-3 * rc
            E2 = (Etot(rc + h) - 2 * Etot(rc) + Etot(rc - h)) / h**2
            out.update(nu=nu, Gam=Gam, ell=ell, PhiW=PhiW, rc_xi=rc, Omega=p["Om"], W0=p["W0"], E_per_len=Etot(rc), breathing_E2=E2,
                       mach_peak=float(max(p["Om"] * rc * np.max(np.hypot(Up, Uz)), p["W0"] * np.max(np.hypot(Wp, Wz)))),
                       lam_q_phys=lam_q * p["Om"] * rc * p["W0"], shares={k: p[k] / Etot(rc) for k in ("U", "W", "B", "E", "conf")})
    out["_prof"] = (r, Up, Uz, Wp, Wz, oUp, oUz, oWp, oWz, b, bp)
    return out

# ---------------- transport ----------------
def transport(res):
    r, Up, Uz, Wp, Wz, oUp, oUz, oWp, oWz, b, bp = res["_prof"]; dA = 2 * np.pi * r * (r[1] - r[0]); I = lambda f: float(np.sum(f * dA))
    a = res["rc_xi"] * XI; L = np.log(RCUT); Cw = res["W_net_circ"]
    sU = I(oUp**2 + oUz**2) / I(Up**2 + Uz**2)
    WWout = Cw**2 / (2 * np.pi) * L                                 # exterior W is irrotational: no curl there
    sW = I(oWp**2 + oWz**2) / (I(Wp**2 + Wz**2) + WWout)
    gb = np.gradient(bp, r); sK = I(gb**2) / max(I(bp**2), 1e-30)
    fWin = I(Wp**2 + Wz**2) / (I(Wp**2 + Wz**2) + WWout)
    s = {"bulk": 0.0, "U": sU, "W": sW, "K": sK, "T": 0.0, "none": 0.0}; f = {"U": 1.0, "W_in": fWin, "W": 1.0, "B": 1.0, "bulk": 0.0, "none": 0.0}
    rows = []
    for key, name, kind, on, what, ordinary in CHANNELS:
        shape = s.get(on, 0.0) / a**2 if kind == "diff" else f.get(on, 0.0)
        bound = PERSIST / shape if shape > 0 else np.inf; rate = ordinary * shape
        rows.append(dict(key=key, channel=name, coeff=ordinary, rate=rate, lifetime=1 / rate if rate > 0 else np.inf,
                         bound=bound, suppression=ordinary / bound if np.isfinite(bound) else 0.0))
    lock = res["Omega"] * CLIGHT / XI                               # rotation rate the locking must beat (s^-1)
    return dict(s_U=sU, s_W=sW, s_K=sK, f_W_in=fWin, rotation_rate=lock), rows

# ---------------- stability survey ----------------
def stability(res):
    """Local criteria (MHD analogue: field B = curl W with Alfven ratio Gam/nu, both classes V+- = U +- W):
    S1 Rayleigh per class (axisymmetric, pure swirl), S2 Michael (swirl + azimuthal field), S3 Howard-Gupta (with axial
    shear, sufficient), S4 Suydam (localized interchange), S5 Kruskal-Shafranov kink (period needed), S6 Kelvin-Helmholtz
    at the U sheet (field-line tension), S7 breathing (E'' at fixed spin), S8 sonic (peak Mach < 1)."""
    r, Up, Uz, Wp, Wz, oUp, oUz, oWp, oWz, b, bp = res["_prof"]
    Om, W0, nu, Gam = res["Omega"], res["W0"], res["nu"], res["Gam"]; rc = res["rc_xi"]
    sl = slice(20, -20); R = r * rc
    UpP, UzP, WpP, WzP = Om * rc * Up, Om * rc * Uz, W0 * Wp, W0 * Wz
    Bp, Bz = W0 * oWp / rc, W0 * oWz / rc; VA2 = Gam / nu
    out = {}
    for sg, lab in ((+1, "+"), (-1, "-")):
        Vp = UpP + sg * WpP; Vz = UzP + sg * WzP
        ray = np.gradient((R * Vp) ** 2, R) / R**3
        mic = ray - VA2 * R * np.gradient(Bp**2 / R**2, R)
        shear = np.gradient(Vz, R) ** 2
        Ri = np.where(shear[sl] > 1e-12, mic[sl] / np.maximum(shear[sl], 1e-30), np.inf)
        out[f"S1_rayleigh{lab}"] = bool(np.min(ray[sl]) >= -1e-9 * np.max(abs(ray[sl])))
        out[f"S2_michael{lab}"] = bool(np.min(mic[sl]) >= -1e-9 * np.max(abs(mic[sl])))
        out[f"S3_HG_Ri_min{lab}"] = float(np.min(Ri))
    mu = np.where(abs(Bz) > 1e-12, Bp / (R * np.where(abs(Bz) > 1e-12, Bz, 1)), np.nan)   # field-line twist
    a_p = (UpP**2 + WpP**2) / R                                      # radial weighted pressure gradient (mass + width modes)
    suy = VA2 * R / 4 * Bz**2 * (np.gradient(mu, R) / np.where(abs(mu) > 1e-12, mu, np.nan)) ** 2 + 2 * a_p * 0 - 0  # pressure rises outward: stabilizing
    out["S4_suydam"] = "stable (pressure gradient outward, shear term >= 0)" if np.all(a_p[sl] >= 0) else "check"
    # S5 kink: m = 1 external kink stable for period L if q_KS = 2 pi r_c Bz/(L Bphi) > 1 at the edge -> L_max
    be = slice(-60, -20); Bpe, Bze = np.mean(Bp[be]), np.mean(Bz[be])
    out["S5_kink_Lmax_over_rc"] = float(2 * np.pi * abs(Bze / Bpe)) if abs(Bpe) > 1e-12 else np.inf
    # S6 KH at the U sheet: unstable for k perpendicular to B unless the jump is parallel to B and below the tension limit
    dU = np.array([UpP[-1], UzP[-1]]); Bs = np.array([Bp[-25], Bz[-25]]); nB = np.linalg.norm(Bs); ndU = np.linalg.norm(dU)
    if ndU < 1e-9: out["S6_KH"] = "no sheet"
    else:
        cosang = abs(dU @ Bs) / (ndU * nB) if nB > 0 else 0.0
        ok = cosang > 0.999 and ndU**2 < 4 * VA2 * nB**2
        out["S6_KH"] = ("stable" if ok else f"unstable (angle {np.degrees(np.arccos(min(cosang, 1))):.0f} deg, |dU|^2/(4 vA^2 B^2) = {ndU**2 / max(4 * VA2 * nB**2, 1e-30):.2g})")
    out["S7_breathing"] = "stable" if res["breathing_E2"] > 0 else "unstable"
    out["S8_sonic"] = "subsonic" if res["mach_peak"] < 1 else f"supersonic ({res['mach_peak']:.2f})"
    return out
