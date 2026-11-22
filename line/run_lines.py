"""Full line-vortex analysis: scorecard, sizing, transport, stability. Writes line/LINE_REPORT.md."""
import sys, json, numpy as np
sys.path.insert(0, ".")
from pipeline.line import registry, analyze, transport, stability
MEDIA = [(nu, G) for nu in (1, 10, 100) for G in (0.01, 1, 100)]
L = ["# Infinite line vortex: scorecard, transport and stability", "",
     "Fields depend on r only, so every target is an exact steady state of the ideal two-class fluid. Core radius r_c is selected by energy at fixed spin per length (ell = 1) and W flux (Phi_W = 1), with outer log cutoff 100 r_c. Units: xi, c.", ""]
L += ["## Shapes", "", "| target | charge/len | spin j | axial P_z | W flux | W net circ | pitch mid | pitch edge | U sheet |", "|---|---|---|---|---|---|---|---|---|"]
reg = registry(); allres = {}
for n, t in reg.items():
    res = analyze(t)
    L.append(f"| {n} | {res['lam_q']:+.3g} | {res['spin_j']:+.3g} | {res['axial_Pz']:+.3g} | {res['W_axial_flux']:+.3g} | {res['W_net_circ']:+.3g} | {res['pitch_mid']:.3g} | {res['pitch_edge']:.3g} | {res['sheet_U']:.2f} |")
L += ["", "## Sizing over media (nu, Gamma): r_c/xi, minima count, peak Mach", "", "| target | " + " | ".join(f"nu={a},G={b}" for a, b in MEDIA) + " |", "|---" * (len(MEDIA) + 1) + "|"]
best = {}
for n, t in reg.items():
    cells = []
    for nu, G in MEDIA:
        res = analyze(t, nu=nu, Gam=G)
        if "rc_xi" in res:
            cells.append(f"{res['rc_xi']:.2g} ({res['n_minima']}) M{res['mach_peak']:.2f}")
            if res["mach_peak"] < 1 and res["n_minima"] == 1 and (n not in best or res["mach_peak"] > best[n]["mach_peak"]): best[n] = res
        else: cells.append("none")
    L.append(f"| {n} | " + " | ".join(cells) + " |")
L += ["", "Cell: r_c/xi (number of minima) and peak Mach. One minimum with Mach < 1 = single stable subsonic core radius.", ""]
L += ["## Stability survey (each target in its subsonic single-minimum medium closest to Mach 1)", "",
      "| target | medium | S1 Rayleigh +/- | S2 Michael +/- | S3 HG Ri min +/- | S4 Suydam | S5 kink L_max/r_c | S6 KH (sheet) | S7 breathing | S8 sonic |", "|---|---|---|---|---|---|---|---|---|---|"]
tr_rows = {}
for n, res in best.items():
    st = stability(res)
    L.append(f"| {n} | nu={res['nu']:g}, G={res['Gam']:g}, r_c={res['rc_xi']:.2g} | {st['S1_rayleigh+']}/{st['S1_rayleigh-']} | {st['S2_michael+']}/{st['S2_michael-']} | "
             f"{st['S3_HG_Ri_min+']:.2g}/{st['S3_HG_Ri_min-']:.2g} | {st['S4_suydam']} | {st['S5_kink_Lmax_over_rc']:.3g} | {st['S6_KH']} | {st['S7_breathing']} | {st['S8_sonic']} |")
    sf, rows = transport(res); tr_rows[n] = (sf, rows)
L += ["", "Ri < 1/4 does not prove instability (Howard-Gupta is sufficient only). S5: the line is kink-stable only if its period (ring circumference) is shorter than L_max.", ""]
L += ["## Transport: suppression needed (ordinary collisional value / value allowed by persistence)", ""]
chans = [x["channel"] for x in next(iter(tr_rows.values()))[1]]
L += ["| target | " + " | ".join(chans) + " |", "|---" * (len(chans) + 1) + "|"]
for n, (sf, rows) in tr_rows.items():
    L.append(f"| {n} | " + " | ".join(f"{x['suppression']:.1e}" for x in rows) + " |")
L += ["", "| target | s_U | s_W | s_K | W energy inside | rotation rate (s^-1) |", "|---|---|---|---|---|---|"]
for n, (sf, rows) in tr_rows.items():
    L.append(f"| {n} | {sf['s_U']:.3g} | {sf['s_W']:.3g} | {sf['s_K']:.3g} | {sf['f_W_in']:.2f} | {sf['rotation_rate']:.2e} |")
open("line/LINE_REPORT.md", "w").write("\n".join(L)); print("\n".join(L))
