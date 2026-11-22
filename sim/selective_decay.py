"""T3: selective decay from a random state. Expect E to decay faster than K± and the state to approach canonical
Beltrami (V± x Om± -> 0), i.e. a double-Beltrami relaxed state. Also: what the drag channels do instead."""
import sys, time, json, numpy as np
sys.path.insert(0, ".")
from sim.hall2f import Solver, random_field
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
runs = {"viscous (nu_U = nu_W = 0.01)": dict(nu_U=0.01, nu_W=0.01),
        "viscous + inter-class drag (r_cl = 0.05)": dict(nu_U=0.01, nu_W=0.01, r_cl=0.05),
        "viscous + U-bulk drag (r_Ub = 0.05)": dict(nu_U=0.01, nu_W=0.01, r_Ub=0.05)}
res = {}
for name, co in runs.items():
    S = Solver(N=32, eta=1.5, **co); Uh, Wh = random_field(S, 1, 5, amp=1.0, seed=3); d0 = S.diag(Uh, Wh)
    t, hist, nxt, t0 = 0.0, [], 0.0, time.time()
    while t < 120:
        if t >= nxt: d = S.diag(Uh, Wh); d["t"] = t; hist.append(d); nxt += 2.0
        dt = S.cfl_dt(Uh, Wh, 0.4); Uh, Wh = S.step(Uh, Wh, min(dt, 0.5)); t += min(dt, 0.5)
    d = S.diag(Uh, Wh); d["t"] = t; hist.append(d)
    # shell spectrum of final state
    k = np.sqrt(S.K2); sp = [float(np.sum((abs(Uh)**2 + S.eta * abs(Wh)**2)[:, (k > n - .5) & (k <= n + .5)])) for n in range(1, 6)]
    res[name] = dict(hist=hist, spectrum=sp)
    print(f"{name} ({time.time()-t0:.0f}s): E/E0 = {d['E']/d0['E']:.3e}, K+/K+0 = {d['Kp']/d0['Kp']:.3e}, K-/K-0 = {d['Km']/d0['Km']:.3e}, "
          f"Beltrami residual +/- = {d['beltrami_p']:.1e}/{d['beltrami_m']:.1e}, EW/EU = {d['EW']/max(d['EU'],1e-300):.3g}, shell fractions = {np.round(np.array(sp)/sum(sp),3)}", flush=True)
json.dump(res, open("sim/selective_decay.json", "w"), default=float)
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for name, r in res.items():
    h = r["hist"]; t = [x["t"] for x in h]; E0 = h[0]["E"]
    ax[0].semilogy(t, [x["E"] / E0 for x in h], label=name + ": E")
    ax[0].semilogy(t, [abs(x["Kp"] / h[0]["Kp"]) for x in h], "--", label=name + ": |K+|")
    ax[1].semilogy(t, [x["beltrami_p"] for x in h], label=name)
ax[0].set_xlabel("t"); ax[0].set_title("energy vs canonical helicity"); ax[0].legend(fontsize=6)
ax[1].set_xlabel("t"); ax[1].set_title("canonical Beltrami residual |V+ x Om+|^2"); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig("sim/selective_decay.png", dpi=110)
