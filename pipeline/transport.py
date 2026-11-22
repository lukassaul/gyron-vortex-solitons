"""Transport coefficients: decay rates of each candidate, channel by channel (first order: dissipative terms
evaluated on the ideal solution). Units: lengths in xi, time in xi/c. Diffusivities in c*xi, rates in c/xi.
Persistence target: every rate < 2e-18 s^-1 (atoms ~4e17 s)."""
import numpy as np
from .analysis import vec, curl, dot

XI, C = 3.8616e-13, 2.99792458e8          # healing length (m), c (m/s)
T0 = XI / C                                # time unit, s
PERSIST = 2e-18                            # s^-1
ELL = 1e-32                                # collisional mean free path (m)
NU_C = 1e40                                # collision rate (s^-1)

# (key, symbol, kind, acts on, what it does, ordinary collisional value in SI)
CHANNELS = [
    ("nu_b",   "bulk shear viscosity",            "diff", "bulk",  "shears bulk vorticity; ring flow is irrotational (n = 0)", C * ELL),
    ("zeta_b", "bulk volume viscosity",           "diff", "bulk",  "damps compression; zero for a static ring",              C * ELL),
    ("nu_U",   "U shear viscosity",               "diff", "U",     "unwinds the spin-carrying drift",                         C * ELL),
    ("nu_W",   "W diffusivity (resistive)",       "diff", "W",     "diffuses the correlation flux, field and moment",         C * ELL),
    ("D_K",    "spin-pressure diffusion (Fick)",  "diff", "K",     "spreads the charge density inside the core",              C * ELL),
    ("chi_T",  "weighted heat conduction",        "diff", "T",     "isothermal closure: no steady gradient to relax",         C * ELL),
    ("r_Ub",   "U-bulk drag (mutual friction)",   "rate", "U",     "locks U to the bulk; converts it to bulk vorticity",      NU_C),
    ("r_cl",   "inter-class drag (+w vs -w)",     "rate", "W_in",  "damps the class counterflow 2W inside the core",           NU_C),
    ("r_F",    "K-current relaxation (C_F ~ -F)", "rate", "W",     "relaxes the static current; moment and field decay",      NU_C),
    ("r_s",    "transverse spin relaxation",      "rate", "B",     "relaxes the locked spin S = I2 B",                        NU_C),
    ("r_H",    "inclination relaxation",          "rate", "none",  "ring carries no inclination structure at this order",     NU_C),
    ("r_par",  "axial-spin relaxation",           "rate", "none",  "gravity sector; eps < 1e-90",                             1e-90 * NU_C),
]

def shape_factors(res):
    """Dimensionless shape numbers (a = 1): diffusive s = <|curl X|^2>/<|X|^2>, rate fractions f."""
    g, chiW, IW, chiU, IU = res["_fields"]
    W, U = vec(chiW, IW, g), vec(chiU, IU, g); wW, wU = curl(W, g), curl(U, g)
    m = g.deep & (g.RR > 0.03); ins = g.inside; out = ~ins & (g.rho > 1.05)
    I = lambda f, msk: float(np.sum((f * g.dV)[msk]))
    EU, EWin, EWout = I(dot(U, U), ins), I(dot(W, W), ins), I(dot(W, W), out)
    sU = I(dot(wU, wU), m) / max(I(dot(U, U), m), 1e-30)
    sW = (I(dot(wW, wW), m) + I(dot(wW, wW), out)) / max(I(dot(W, W), m) + EWout, 1e-30)
    b = -dot(U, W) / res["eta"]; gb = np.stack([np.gradient(b, g.dr, axis=0), np.gradient(b, g.dz, axis=1)])
    sK = I(np.sum(gb**2, axis=0), m) / max(I(b * b, m), 1e-30)
    iz = np.argmin(abs(g.Z)); k1 = np.searchsorted(g.R, 0.97)          # U tangential jump at the core edge (sheet)
    jump = float(np.sqrt(U[1][k1, iz]**2 + U[2][k1, iz]**2) / max(np.sqrt(np.nanmax(dot(U, U)[m])), 1e-30))
    fWin = EWin / max(EWin + EWout, 1e-30)
    return dict(s_U=sU, s_W=sW, s_K=sK, f_W_in=fWin, sheet_jump=jump, omega_rot=float(np.sqrt(np.nanmax(dot(U, U)[m]))))

def channel_table(res, coeffs=None):
    """Rates (s^-1), lifetimes and persistence bounds for every channel, at the sized ring (a_xi, A)."""
    sf = shape_factors(res); a = res["a_xi"] * XI
    s = {"bulk": 0.0, "U": sf["s_U"], "W": sf["s_W"], "K": sf["s_K"], "T": 0.0, "none": 0.0}
    f = {"U": 1.0, "W_in": sf["f_W_in"], "W": 1.0, "B": 1.0, "bulk": 0.0, "none": 0.0}
    rows = []
    for key, name, kind, on, what, ordinary in CHANNELS:
        val = (coeffs or {}).get(key, ordinary)
        shape = s.get(on, 0.0) / a**2 if kind == "diff" else f.get(on, 0.0)      # rate per unit coefficient
        rate = val * shape
        bound = PERSIST / shape if shape > 0 else np.inf
        rows.append(dict(key=key, channel=name, kind=kind, acts_on=on, effect=what, coeff=val, rate=rate,
                         lifetime=(1 / rate if rate > 0 else np.inf), bound=bound,
                         suppression=(val / bound if np.isfinite(bound) else 0.0)))
    lock = sf["omega_rot"] * res["A"] * C / a            # locking must beat the rotation rate
    return sf, rows, lock
