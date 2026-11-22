"""Scorecard analysis of a Target. Lengths in units of the core radius a (then in xi after sizing); velocities in c.
Medium: nu (weighted/bulk inertia), Gam (field stiffness), p0 = 5/7, sig = 2 sqrt2/3 (bubble confinement)."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from scipy.interpolate import RegularGridInterpolator as RGI

P0, SIG = 5 / 7, 2 * np.sqrt(2) / 3

class Grid:
    def __init__(self, n=260, L=1.6):
        self.R = np.linspace(L / (2 * n), L, n); self.Z = np.linspace(-L, L, 2 * n)
        self.RR, self.ZZ = np.meshgrid(self.R, self.Z, indexing="ij")
        self.rho = np.hypot(self.RR, self.ZZ); self.s2 = (self.RR / self.rho) ** 2
        self.dr = self.R[1] - self.R[0]; self.dz = self.Z[1] - self.Z[0]
        self.inside = self.rho < 1; self.deep = self.rho < 0.95
        self.dV = 2 * np.pi * self.RR * self.dr * self.dz

def stream(G, H, grid, exterior):
    rho, s2 = grid.rho, grid.s2
    chi = np.where(grid.inside, G(np.minimum(rho, 1)) * s2, (G(1.0) * s2 / rho) if exterior else 0.0)
    I = np.where(grid.inside, H(np.minimum(rho, 1)) * s2, (H(1.0) * s2 / rho) if exterior else 0.0)
    return chi, I

def vec(chi, I, g):
    cr = np.gradient(chi, g.dr, axis=0); cz = np.gradient(chi, g.dz, axis=1)
    return np.stack([-cz / g.RR, I / g.RR, cr / g.RR])          # (r, phi, z)

def curl(F, g):
    Fr, Fp, Fz = F
    return np.stack([-np.gradient(Fp, g.dz, axis=1),
                     np.gradient(Fr, g.dz, axis=1) - np.gradient(Fz, g.dr, axis=0),
                     np.gradient(g.RR * Fp, g.dr, axis=0) / g.RR])

def cross(A, B):
    return np.stack([A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2], A[0] * B[1] - A[1] * B[0]])

def dot(A, B): return np.sum(A * B, axis=0)

def winding(chi, I, g):
    """q on the middle closed surface of W (nan if no closed band)."""
    iz = np.argmin(abs(g.Z)); prof = np.where(g.inside[:, iz], chi[:, iz], np.nan)
    if np.nanmax(prof) <= 0: chi, I, prof = -chi, -I, -prof
    iO = np.nanargmax(prof); chiO = prof[iO]; chiB = max(np.nanmax(np.where(g.rho[:, iz] > 0.98, prof, -np.inf)), 0)
    if chiO <= chiB * 1.001 or iO == 0: return np.nan, np.nan
    lev = chiB + 0.5 * (chiO - chiB)
    cs = plt.figure().add_subplot().contour(g.RR, g.ZZ, np.where(g.inside, chi, np.nan), levels=[lev]); plt.close("all")
    ps = [p for p in cs.allsegs[0] if len(p) > 20 and np.allclose(p[0], p[-1], atol=5e-3)]
    if not ps: return np.nan, g.R[iO]
    p = max(ps, key=len); rr, zz = p[:, 0], p[:, 1]; dl = np.hypot(np.diff(rr), np.diff(zz))
    rm, zm = (rr[1:] + rr[:-1]) / 2, (zz[1:] + zz[:-1]) / 2
    gr = np.gradient(chi, g.dr, axis=0); gz = np.gradient(chi, g.dz, axis=1)
    f = lambda F: RGI((g.R, g.Z), np.nan_to_num(F), bounds_error=False, fill_value=0)(np.c_[rm, zm])
    return np.sum(f(I) / (rm * np.hypot(f(gr), f(gz))) * dl) / (2 * np.pi), g.R[iO]

def analyze(t, nu=10.0, Gam=1.0, grid=None):
    g = grid or Grid()
    chiW, IW = stream(t.GW, t.HW, g, True); W = vec(chiW, IW, g)
    zero = lambda x: 0 * x
    chiU, IU = stream(t.GU or zero, t.HU or zero, g, False); U = vec(chiU, IU, g)
    wW, wU = curl(W, g), curl(U, g)
    m = g.deep & (g.RR > 0.03)
    def resid(L, scale):
        c = curl(L, g); return float(np.sqrt(np.mean(dot(c, c)[m])) / max(scale, 1e-30))
    sc = lambda A, B: np.sqrt(np.mean((dot(A, A) * dot(B, B))[m]))
    R0 = resid(cross(U, wU), sc(U, wU)) if t.GU else 0.0
    R1 = resid(cross(U, wW) + cross(W, wU), sc(U, wW) + sc(W, wU)) if t.GU else 0.0
    R2 = resid(cross(W, wW), sc(W, wW))
    # multi-fluid (Hall-like) steadiness: classes V± = U ± W, canonical momentum P± = V± ± (eta-1) W (h = 1, w A = (eta-1) W)
    Rc = []
    for sgn in (+1, -1):
        V = U + sgn * W; Om = wU + sgn * t.eta * wW
        Rc.append(resid(cross(V, Om), sc(V, Om)))
    q, rO = winding(chiW, IW, g)
    ins = g.inside
    impulse = float(np.sum((U[2] * g.dV)[ins])); j = float(np.sum((IU * g.dV)[ins])) if t.GU else 0.0
    D = float(t.HW(1.0)); G1 = float(t.GW(1.0))
    th = np.linspace(0.02, np.pi - 0.02, 300); pts = np.c_[0.985 * np.sin(th), 0.985 * np.cos(th)]
    uw = RGI((g.R, g.Z), -dot(U, W) / t.eta, bounds_error=False, fill_value=0)(pts)
    Q = float(4 * np.pi * np.trapezoid(uw * 0.5 * np.sin(th), th))
    K_in = float(np.sum(((dot(U, U) + dot(W, W)) * g.dV)[ins]))
    K_out = float(np.sum((dot(W, W) * g.dV)[~ins])) + (4 * np.pi / 3) * (G1**2 + 0.5 * D**2) / 1.6**3  # grid + dipole tail
    Bf = 0.5 * float(np.sum((dot(wW, wW) * g.dV)[m])) + (4 * np.pi / 3) * D**2      # interior + exterior dipole (a = 1)
    peak = float(np.nanmax(np.where(m, np.maximum(np.sqrt(dot(U, U)), np.sqrt(dot(W, W))), 0)))
    rig = np.nan
    if t.GU and abs(j) > 0:
        Om = j / float(np.sum((g.RR**2 * g.dV)[ins])); rig = float(Om**2 * np.sum((g.RR**2 * g.dV)[ins]) / max(np.sum((dot(U, U) * g.dV)[ins]), 1e-30))
    out = dict(name=t.name, family=t.family, eta=t.eta, params=t.params, steady_R0=R0, steady_R1=R1, steady_R2=R2, steady_canon=max(Rc),
               q=q, r_O=rO, impulse=impulse, spin_j=j, moment_D=D, charge_Q=Q, K=K_in + K_out, B_energy=Bf,
               peak_per_A=peak, rigid_fraction=rig, orientation=("right" if D * j * Q > 0 else "wrong" if D * j * Q < 0 else "n/a"))
    if abs(j) > 1e-9 and abs(Q) > 1e-12:
        kap = 50.0
        for _ in range(300):
            E = lambda a: kap**2 * out["K"] / (8 * nu * a**5 * j**2) + Gam * kap**2 * Bf / (4 * nu**2 * a**7 * j**2) + P0 * 4 * np.pi / 3 * a**3 + 4 * np.pi * SIG * a**2
            r = minimize_scalar(lambda la: E(np.exp(la)), bounds=(np.log(1e-3), np.log(1e4)), method="bounded")
            a = float(np.exp(r.x)); En = E(a)
            if abs(En - kap) < 1e-10 * kap: break
            kap = 0.5 * kap + 0.5 * En
        A = kap / (2 * nu * a**4 * abs(j)); gfac = 16 * np.pi * abs(D) * a / (abs(Q) * A)
        Qt = abs(Q) * A**2 * a
        out.update(nu=nu, Gam=Gam, a_xi=a, A=A, kappa=kap, g=gfac, mu_over_muB=gfac / 2, mach_peak=A * peak,
                   valid_subsonic=bool(A * peak < 1), alpha=Gam * Qt**2 / (4 * np.pi * kap))
    out["_fields"] = (g, chiW, IW, chiU, IU)
    return out

def figure(res, path):
    g, chiW, IW, chiU, IU = res["_fields"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 5.6))
    for a_, chi, I, ttl in ((ax[0], chiW, IW, "W: correlation flux (field, moment)"), (ax[1], chiU, IU, "U: mass drift (spin)")):
        Fp = I / g.RR; mm = np.nanpercentile(abs(Fp), 99) or 1
        a_.pcolormesh(g.RR, g.ZZ, Fp, cmap="RdBu_r", vmin=-mm, vmax=mm, shading="auto")
        if np.nanmax(abs(chi)) > 0: a_.contour(g.RR, g.ZZ, chi, levels=22, colors="k", linewidths=0.6)
        t = np.linspace(-np.pi / 2, np.pi / 2, 100); a_.plot(np.cos(t), np.sin(t), "g--", lw=1)
        a_.set_aspect("equal"); a_.set_title(ttl, fontsize=10); a_.set_xlabel("r / a"); a_.set_ylabel("z / a")
    fig.suptitle(res["name"], fontsize=11); plt.tight_layout(); plt.savefig(path, dpi=90); plt.close(fig)
