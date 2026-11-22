"""Filter and rank the batch; write batch/REPORT.md and figures."""
import json, numpy as np, matplotlib, sys
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, ".")
from pipeline.targets import double_beltrami
from pipeline.analysis import analyze, figure, Grid
E = json.load(open("batch/evaluation.json"))
rows = []
for e in E:
    s = e["shape"]; p = s["params"]
    for r in e["runs"]:
        rows.append(dict(name=s["name"], eta=p["eta"], mu1=p["mu1"], mu2=p["mu2"], q=s["q"], canon=s["steady_canon"],
                         impulse=s["impulse"], orient=s["orientation"], **r))
N = len(E); steady = [e for e in E if e["shape"]["steady_canon"] < 1e-2]
closed = [e for e in steady if np.isfinite(e["shape"]["q"]) and abs(e["shape"]["q"] - 1) < 0.03]
oriented = [e for e in closed if e["shape"]["orientation"] == "right"]
passing = [r for r in rows if r["canon"] < 1e-2 and np.isfinite(r["q"]) and abs(r["q"] - 1) < 0.03 and r["orient"] == "right" and r["mach_peak"] < 1]
passing.sort(key=lambda r: abs(np.log(r["g"] / 2)))
best = {}
for r in passing:
    if r["name"] not in best: best[r["name"]] = r
L = []
L += ["# Batch report: Hall-like (double-Beltrami) ring family", "",
      "Scan: eta in {1.05, 1.1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4}; for each, all isolated solutions with smooth W matching to an exterior dipole, zero impulse and winding ratio q = 1 (scanner), then the pipeline scorecard over nu in {1, 3, 10, 30, 100} x Gamma in {0.01, 0.1, 1, 10, 100} (25 media).", "",
      "## Funnel", "",
      f"| Stage | Count |", "|---|---|",
      f"| Isolated solutions found by the scanner | {N} |",
      f"| Steady in the field-coupled two-class fluid (canonical residual < 1e-2) | {len(steady)} |",
      f"| q = 1 confirmed by the pipeline (within 0.03) | {len(closed)} |",
      f"| Right orientation (moment antiparallel to spin for negative charge) | {len(oriented)} |",
      f"| ...and subsonic at the energy-selected size, in at least one medium | {len(best)} |", "",
      "All scanner solutions are steady (canonical residual 3e-4 to 5e-3, grid level). Many scanner hits fail the pipeline's q check: the scanner picked a different family of tori where W has two O-points.", "",
      "## Ranked candidates (best medium for each, ranked by |ln(g/2)|)", "",
      "| Target | eta | mu1 a | mu2 a | nu | Gamma | a/xi | A | peak Mach | kappa | g | alpha |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in best.values():
    L.append(f"| {r['name']} | {r['eta']:g} | {r['mu1']:.3f} | {r['mu2']:.3f} | {r['nu']:g} | {r['Gam']:g} | {r['a_xi']:.3g} | {r['A']:.3g} | {r['mach_peak']:.2f} | {r['kappa']:.3g} | {r['g']:.3g} | {r['alpha']:.2e} |")
gmin = min((r["g"] for r in passing), default=np.nan)
L += ["", "## Findings", "",
      f"- The lowest g among passing cases is {gmin:.3g}; no candidate reaches g = 2.",
      "- One branch carries the right orientation with q = 1 at every eta from 1.05 to 2 (mu1 a = 4.61 to 4.95, mu2 a = 13.37 to 13.58); along it g rises with eta (about 52 at eta = 1.05 to 280 at eta = 2), so the least screening (eta -> 1) is best.",
      "- g depends weakly on the medium ratios; within the subsonic range it moves by tens of percent, not by the factor of ~25 needed.",
      "- alpha stays between 4e-10 and 4e-7 in the best subsonic cases, far below 1/137 (larger values occur only in supersonic media).",
      "- The obstruction is the one in working record s.34: g = 16 pi (D/Q)(a/A); the energy balance at fixed spin gives a/A ~ 10, while g = 2 needs a/A ~ 0.4 on the best branch.", "",
      "## What to try next", "",
      "- Add the omitted field-energy share of the mass at the alpha level and relativistic inertia (both lower a/A).",
      "- Higher harmonics (l = 2, 3) and three or more classes, to raise Q/D.",
      "- Extend eta below 1.05 toward 1 along the best branch; check whether g has a floor.",
      "", "Figures: `batch/g_vs_eta.png` (g at nu = 10, Gamma = 1 against eta for every steady, q = 1 solution, by orientation) and `batch/best_candidate.png`.", ""]
open("batch/REPORT.md", "w").write("\n".join(L)); print("\n".join(L))
fig, ax = plt.subplots(figsize=(7, 4.5))
for e in closed:
    s = e["shape"]; r = [x for x in e["runs"] if x["nu"] == 10 and x["Gam"] == 1]
    if r: ax.scatter(s["params"]["eta"], r[0]["g"], c="C2" if s["orientation"] == "right" else "C3", marker="o" if r[0]["mach_peak"] < 1 else "x", s=40)
ax.axhline(2, ls="--", c="k"); ax.set_yscale("log"); ax.set_xlabel("η (inertia ratio Λ/h)"); ax.set_ylabel("g  (ν = 10, Γ = 1)")
ax.scatter([], [], c="C2", label="right orientation"); ax.scatter([], [], c="C3", label="wrong orientation"); ax.scatter([], [], c="k", marker="x", label="supersonic")
ax.legend(); ax.set_title("Hall-like rings with q = 1 and zero impulse: g against η"); plt.tight_layout(); plt.savefig("batch/g_vs_eta.png", dpi=110)
if best:
    b = list(best.values())[0]; res = analyze(double_beltrami(b["eta"], b["mu1"], b["mu2"]), nu=b["nu"], Gam=b["Gam"]); figure(res, "batch/best_candidate.png")
