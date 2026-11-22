"""Where in parameter space does a stable, subsonic, thin Carnot ring exist? Reduced thin-ring model (units xi, c).
Conditions: (1) balance W = D exists; (2) Carnot-stable: a* < a_E; (3) subsonic core: V = n/a* < 1 (choke, a* > n xi);
(4) subsonic toroidal spin flow U_t = J/(2 pi^2 nu_i a*^2 R*^2) < 1; (5) thin ring R*/a* >= 2. (R equilibrium is convex: always stable.)"""
import itertools, json, numpy as np
from scipy.optimize import brentq
P0, SIG = 5 / 7, 2 * np.sqrt(2) / 3
def solve(n, nu_i, nu_v, C, J, p0=P0, sig=SIG):
    Gam = 2 * np.pi * n
    A = lambda a: nu_i * Gam**2 / (16 * np.pi) + p0 * np.pi * a**2 + 2 * np.pi * sig * a
    Req = lambda a: (3 * J**2 / (8 * np.pi**3 * nu_i * a**2 * A(a)))**0.25
    E = lambda a: (8 * np.pi / 3) * Req(a) * A(a); D = lambda a: 2 * nu_i * nu_v * Req(a) * Gam**2 / a**2
    a = np.geomspace(1e-3, 1e3, 3000); Ea = E(a); aE = a[np.argmin(Ea)]; f = lambda x: C * E(x) - D(x)
    sols = []
    for i in np.nonzero(np.diff(np.sign(f(a))))[0]:
        x = brentq(f, a[i], a[i + 1]); h = 1e-5 * x; lam = (f(x + h) - f(x - h)) / (E(x + h) - E(x - h)); R = Req(x)
        sols.append(dict(a=x, R=R, E=E(x), stable=bool(lam < 0), V=n / x, Ut=J / (2 * np.pi**2 * nu_i * x**2 * R**2), aE=aE))
    return sols
grid = dict(n=[1/3, 1, 2], nu_i=[1, 3, 10, 30, 100], nu_v=[1e-4, 1e-3, 1e-2], C=[1e-3, 1e-2, 1e-1, 1], J=[10, 30, 100, 300, 1000, 3000])
rows = []; fails = dict(no_balance=0, unstable=0, supersonic_core=0, supersonic_spin=0, not_thin=0)
for n, nu_i, nu_v, C, J in itertools.product(*grid.values()):
    s = solve(n, nu_i, nu_v, C, J)
    if not s: fails["no_balance"] += 1; continue
    best = None
    for x in s:
        ok = [x["stable"], x["V"] < 1, x["Ut"] < 1, x["R"] / x["a"] >= 2]
        if all(ok): best = x; break
    if best: rows.append(dict(n=n, nu_i=nu_i, nu_v=nu_v, C=C, J=J, **best))
    else:
        x = s[0]
        for k, ok in zip(["unstable", "supersonic_core", "supersonic_spin", "not_thin"], [x["stable"], x["V"] < 1, x["Ut"] < 1, x["R"] / x["a"] >= 2]):
            if not ok: fails[k] += 1; break
total = np.prod([len(v) for v in grid.values()])
print(f"{len(rows)} of {total} parameter points give a stable, subsonic, thin Carnot ring. First failure otherwise: {fails}")
for n in grid["n"]:
    r = [x for x in rows if x["n"] == n]; print(f"  n = {n:.3f}: {len(r)} stable points")
    for x in sorted(r, key=lambda x: x["E"])[:4]:
        print(f"     nu_i={x['nu_i']:g} nu_v={x['nu_v']:g} C={x['C']:g} J={x['J']:g}: a*={x['a']:.2f}xi (choke {n:.2f}, a_E={x['aE']:.2f}), R*={x['R']:.1f}xi, R/a={x['R']/x['a']:.1f}, V={x['V']:.2f}, Ut={x['Ut']:.2f}, E={x['E']:.0f}")
# which condition binds, by n
json.dump(rows, open("carnot/stable_map.json", "w"), default=float)
