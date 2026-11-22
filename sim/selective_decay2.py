import sys, time, json, numpy as np
sys.path.insert(0, ".")
from sim.hall2f import Solver, random_field
def helical(S, Fh):
    k = np.sqrt(S.K2s); return 0.5 * (Fh + S.curl_h(Fh) / k)     # keep the curl = +k part
res = {}
for nu in (0.003, 0.0015):
    S = Solver(N=48, eta=1.5, nu_U=nu, nu_W=nu); Uh, Wh = random_field(S, 1, 6, amp=1.0, seed=5)
    Uh, Wh = helical(S, Uh), helical(S, Wh); f = np.sqrt(S.diag(Uh, Wh)["E"]); Uh, Wh = Uh / f, Wh / f
    d0 = S.diag(Uh, Wh); t, nxt, hist, t0 = 0.0, 0.0, [], time.time()
    while t < 40:
        if t >= nxt: d = S.diag(Uh, Wh); d["t"] = t; hist.append(d); nxt += 1.0
        dt = min(S.cfl_dt(Uh, Wh, 0.4), 0.3); Uh, Wh = S.step(Uh, Wh, dt); t += dt
    d = S.diag(Uh, Wh); d["t"] = t; hist.append(d); res[nu] = hist
    k = np.sqrt(S.K2); sp = [float(np.sum((abs(Uh)**2 + S.eta * abs(Wh)**2)[:, (k > n - .5) & (k <= n + .5)])) for n in range(1, 7)]
    print(f"nu={nu} ({time.time()-t0:.0f}s): E/E0={d['E']/d0['E']:.3e}  K+/K+0={d['Kp']/d0['Kp']:.3e}  K-/K-0={d['Km']/d0['Km']:.3e}  "
          f"Beltrami res +/-: {d0['beltrami_p']:.2f}->{d['beltrami_p']:.1e} / {d0['beltrami_m']:.2f}->{d['beltrami_m']:.1e}  shells={np.round(np.array(sp)/sum(sp),3)}", flush=True)
    for x in hist[::8]: print(f"   t={x['t']:5.1f}  E={x['E']:.3e}  K+={x['Kp']:.3e}  K-={x['Km']:.3e}  bel+={x['beltrami_p']:.2e}")
json.dump({str(k): v for k, v in res.items()}, open("sim/selective_decay2.json", "w"), default=float)
