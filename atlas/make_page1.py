import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import carnot_ring as M
A = json.load(open("carnot_ring_atlas.json")); Fs = sorted(set(x["F"] for x in A))
fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
ax[0].bar([str(f) for f in Fs], [sum(1 for x in A if x["F"] == f and x["n_ok"]) / sum(1 for x in A if x["F"] == f) for f in Fs], color="C0")
ax[0].set_xlabel("exterior circulation fraction F"); ax[0].set_ylabel("fraction of (P, S/P, J~) with a stable window"); ax[0].set_title("Driven rings need open circulation")
Ps = sorted(set(x["P"] for x in A)); Js = sorted(set(x["J"] for x in A))
for k, f in enumerate((0.5, 1.0)):
    Z = np.full((len(Ps), len(Js)), np.nan)
    for x in A:
        if x["F"] == f and x["S_over_P"] == 0.3 and x["n_ok"]: Z[Ps.index(x["P"]), Js.index(x["J"])] = np.log10(x["a_max"] / x["a_min"]) if x["a_max"] > x["a_min"] else 0.02
    im = ax[k + 1].pcolormesh(np.log10(Js), np.log10(Ps), Z, cmap="viridis", vmin=0, vmax=1.2, shading="nearest")
    ax[k + 1].set_xlabel("log10 J~ (spin)"); ax[k + 1].set_ylabel("log10 P (confinement)"); ax[k + 1].set_title(f"F = {f}: stable window width log10(a_max/a_min)")
    plt.colorbar(im, ax=ax[k + 1])
plt.tight_layout(); plt.savefig("page1_map.png", dpi=110)
# example cross-sections
fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
for k, f in enumerate((0.5, 1.0)):
    cand = [x for x in A if x["F"] == f and x["n_ok"]]; x = sorted(cand, key=lambda x: x["J"])[len(cand) // 2]; e = x["example"]
    M.F = f; w, Uth, Us, nw, dA, wo, d = M.section(e["a"], e["U0"])
    ax[k].plot(w, Uth, label="poloidal U_th"); ax[k].plot(w, Us, label="toroidal (spin) U_s"); ax[k].plot(w, nw, label="weighted density n_w")
    ax[k].axvline(1, c="k", ls=":"); ax[k].axvline(e["a"], c="gray", ls="--")
    ax[k].set_title(f"F={f}, P={x['P']:.3g}, J~={x['J']:.0e}: a*={e['a']:.2f}, R*={e['R']:.0f}, N={e['N']:.2f}, U_self={e['Uself']:.3f}c", fontsize=9)
    ax[k].set_xlabel("w / r_i (hollow core inside w = 1)"); ax[k].legend(fontsize=8)
plt.tight_layout(); plt.savefig("page1_examples.png", dpi=110)
