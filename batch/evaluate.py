"""Scorecard for every batch solution over a (nu, Gamma) grid. Writes batch/evaluation.json."""
import glob, json, sys, numpy as np
sys.path.insert(0, ".")
from pipeline.targets import double_beltrami
from pipeline.analysis import analyze, Grid, P0, SIG
from scipy.optimize import minimize_scalar
def size(sh, nu, Gam):
    j, K, Bf, D, Q, peak = sh["spin_j"], sh["K"], sh["B_energy"], sh["moment_D"], sh["charge_Q"], sh["peak_per_A"]
    if abs(j) < 1e-9 or abs(Q) < 1e-12: return None
    kap = 50.0
    for _ in range(300):
        E = lambda a: kap**2*K/(8*nu*a**5*j**2) + Gam*kap**2*Bf/(4*nu**2*a**7*j**2) + P0*4*np.pi/3*a**3 + 4*np.pi*SIG*a**2
        r = minimize_scalar(lambda la: E(np.exp(la)), bounds=(np.log(1e-3), np.log(1e4)), method="bounded")
        a = float(np.exp(r.x)); En = E(a)
        if abs(En-kap) < 1e-10*kap: break
        kap = 0.5*kap + 0.5*En
    A = kap/(2*nu*a**4*abs(j)); g = 16*np.pi*abs(D)*a/(abs(Q)*A); Qt = abs(Q)*A**2*a
    return dict(nu=nu, Gam=Gam, a_xi=a, A=A, kappa=kap, g=g, mach_peak=A*peak, alpha=Gam*Qt**2/(4*np.pi*kap))
grid = Grid(); out = []
for f in sorted(glob.glob("batch/solutions_eta*.json")):
    d = json.load(open(f))
    for m1, m2 in d["solutions"]:
        t = double_beltrami(d["eta"], m1, m2); sh = analyze(t, grid=grid); sh.pop("_fields")
        runs = [r for nu in (1, 3, 10, 30, 100) for Gam in (0.01, 0.1, 1, 10, 100) if (r := size(sh, nu, Gam))]
        out.append(dict(shape={k: v for k, v in sh.items() if k not in ("nu","Gam","a_xi","A","kappa","g","mu_over_muB","mach_peak","valid_subsonic","alpha")}, runs=runs))
        print(t.name, f"q={sh['q']:.3f} canon={sh['steady_canon']:.1e} orient={sh['orientation']} g(nu=10,G=1)={sh.get('g',float('nan')):.3g}", flush=True)
json.dump(out, open("batch/evaluation.json", "w"), indent=1, default=float)
