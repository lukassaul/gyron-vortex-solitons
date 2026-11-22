import sys, json, numpy as np
sys.path.insert(0, "."); sys.path.insert(0, "line")
from aligned_ring import ring
rows = []
for kh in (0.3, 0.4, 0.5, 0.6, 0.8):
    for eps in (0.6, 0.8, 1.0, 1.25, 1.5):
        for nu in np.round(np.geomspace(2, 60, 8), 2):
            for G in np.round(np.geomspace(0.03, 30, 7), 3):
                o = ring(kh, eps, nu, G)
                if o: o.pop("stab"); rows.append(o)
json.dump(rows, open("line/aligned_refine.json", "w"), default=float)
ok = [o for o in rows if o["all_local_pass"]]
print(len(rows), len(ok))
for o in sorted(ok, key=lambda o: abs(np.log(o["g"] / 2)))[:12]:
    print(" ".join(f"{k}={o[k]:.3g}" for k in ("kh", "eps", "nu", "Gam", "R_over_rc", "rc_xi", "R_xi", "W0", "mach", "kappa", "g", "alpha")))
print("closest alpha to 1/137 among passing:")
for o in sorted(ok, key=lambda o: abs(np.log(o["alpha"] * 137.036)))[:6]:
    print(" ".join(f"{k}={o[k]:.3g}" for k in ("kh", "eps", "nu", "Gam", "R_xi", "mach", "kappa", "g", "alpha")))
