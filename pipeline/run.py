"""Run the scorecard on targets.  Usage (from the repository root):
    python -m pipeline.run [--targets NAME ...] [--nu 10] [--gam 1] [--out reports]"""
import argparse, json, os, numpy as np
from .targets import registry
from .analysis import analyze, figure, Grid

COLS = [("name", "{}"), ("q", "{:.3f}"), ("steady_R2", "{:.2g}"), ("steady_R1", "{:.2g}"), ("steady_canon", "{:.2g}"), ("impulse", "{:+.3g}"),
        ("spin_j", "{:+.3g}"), ("charge_Q", "{:+.3g}"), ("moment_D", "{:+.3g}"), ("orientation", "{}"),
        ("rigid_fraction", "{:.2f}"), ("a_xi", "{:.3g}"), ("A", "{:.3g}"), ("mach_peak", "{:.2f}"), ("kappa", "{:.3g}"),
        ("g", "{:.3g}"), ("alpha", "{:.2e}")]

def fmt(v, f):
    if v is None: return "—"
    try:
        if isinstance(v, float) and not np.isfinite(v): return "—"
        return f.format(v)
    except Exception: return str(v)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--targets", nargs="*"); ap.add_argument("--nu", type=float, default=10.0)
    ap.add_argument("--gam", type=float, default=1.0); ap.add_argument("--out", default="reports"); a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True); reg = registry(); names = a.targets or list(reg)
    grid = Grid(); rows = []
    for n in names:
        res = analyze(reg[n], nu=a.nu, Gam=a.gam, grid=grid); figure(res, os.path.join(a.out, f"{n}.png"))
        res.pop("_fields"); rows.append(res)
        with open(os.path.join(a.out, f"{n}.json"), "w") as f: json.dump(res, f, indent=1, default=float)
        print(n, "done")
    hdr = "| " + " | ".join(c for c, _ in COLS) + " |\n|" + "---|" * len(COLS) + "\n"
    body = "".join("| " + " | ".join(fmt(r.get(c), f) for c, f in COLS) + " |\n" for r in rows)
    txt = (f"# Scorecard (nu = {a.nu}, Gamma = {a.gam})\n\n"
           "Steady residuals R0, R1, R2 = relative curl of the Lamb terms of the three Euler orders (0 = exactly steady; "
           "R0, R1 need U); steady_canon = same test for the two-class fluid with field coupling (canonical vorticity), the right test for Hall-like targets). Electron targets: q = 1, impulse 0, orientation right, g = 2, mach_peak < 1, alpha = 1/137.\n\n" + hdr + body)
    open(os.path.join(a.out, "scorecard.md"), "w").write(txt); print(txt)

if __name__ == "__main__": main()
