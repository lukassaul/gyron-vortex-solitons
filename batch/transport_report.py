"""Transport scorecard for the batch candidates (each at its best medium) and the registry targets."""
import json, sys, numpy as np
sys.path.insert(0, ".")
from pipeline.targets import double_beltrami, registry
from pipeline.analysis import analyze, Grid
from pipeline.transport import channel_table, CHANNELS, PERSIST
E = json.load(open("batch/evaluation.json")); grid = Grid()
cand = []
for e in E:
    s = e["shape"]
    if s["steady_canon"] < 1e-2 and np.isfinite(s["q"]) and abs(s["q"] - 1) < 0.03 and s["orientation"] == "right":
        ok = [r for r in e["runs"] if r["mach_peak"] < 1]
        if ok: cand.append((s, min(ok, key=lambda r: abs(np.log(r["g"] / 2)))))
L = ["# Transport scorecard", "",
     "Every transport coefficient, evaluated to first order on each candidate (dissipative terms on the ideal solution) at its energy-selected size. "
     "Persistence requires every rate below 2e-18 s^-1. 'Ordinary' = the collisional value (diffusivities c*l with l = 1e-32 m; rates nu_c = 1e40 s^-1). "
     "'Suppression' = ordinary value / allowed value: the factor by which the kernel must suppress that channel.", ""]
L += ["## Channels", "", "| key | channel | kind | acts on | effect |", "|---|---|---|---|---|"]
L += [f"| {k} | {n} | {'diffusivity' if kd == 'diff' else 'rate'} | {on} | {w} |" for k, n, kd, on, w, _ in CHANNELS]
summary = {}
for s, r in cand:
    t = double_beltrami(s["params"]["eta"], s["params"]["mu1"], s["params"]["mu2"])
    res = analyze(t, nu=r["nu"], Gam=r["Gam"], grid=grid); sf, rows, lock = channel_table(res)
    L += ["", f"## {t.name}  (nu = {r['nu']:g}, Gamma = {r['Gam']:g}, a = {res['a_xi']:.3g} xi, g = {res['g']:.3g})", "",
          f"Shape numbers: s_U = {sf['s_U']:.3g}, s_W = {sf['s_W']:.3g}, s_K = {sf['s_K']:.3g} (per a^2); W energy inside {sf['f_W_in']:.2f}; "
          f"U sheet jump at the edge {sf['sheet_jump']:.2f} of peak; rotation rate {lock:.2e} s^-1 (locking must be faster).", "",
          "| channel | ordinary value | ordinary lifetime (s) | allowed value | suppression needed |", "|---|---|---|---|---|"]
    for x in rows:
        L.append(f"| {x['channel']} | {x['coeff']:.2e} | {x['lifetime']:.2e} | {x['bound']:.2e} | {x['suppression']:.1e} |")
        summary.setdefault(x["channel"], []).append(x["suppression"])
L[L.index("## Channels"):L.index("## Channels")] = ["## Summary over candidates (suppression factor needed, min to max)", "",
     "| channel | min | max |", "|---|---|---|"] + [f"| {k} | {min(v):.1e} | {max(v):.1e} |" for k, v in summary.items()] + [""]
open("batch/TRANSPORT.md", "w").write("\n".join(L)); print("\n".join(L[:40]))
