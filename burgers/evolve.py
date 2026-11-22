"""Time evolution of an axisymmetric vortex ring with and without a co-moving conversion sink (does a steady core form?).
  d_t om + u.grad om = om u_r/r - om div u + nu (lap om - om/r^2),  u = curl(psi)/r + grad phi
  D*psi = -r om ;  lap phi = -kappa g,  g = exp(-d^2/s^2) centred on the current vorticity centroid (sink follows the core).
Units delta0 = 1 (initial Gaussian core), c = 1. Track circulation, core radius, ring radius, translation."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl, json, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
Lr, Lz, Nr, Nz = 16.0, 32.0, 128, 256
r = np.linspace(Lr / Nr / 2, Lr, Nr); z = np.linspace(-Lz / 2, Lz / 2, Nz); dr, dz = r[1] - r[0], z[1] - z[0]
RR, ZZ = np.meshgrid(r, z, indexing="ij")
def build(stokes):
    N = Nr * Nz; I = lambda i, j: i * Nz + j; rows, cols, vals = [], [], []
    for i in range(Nr):
        for j in range(Nz):
            k = I(i, j)
            if i == Nr - 1 or j in (0, Nz - 1): rows.append(k); cols.append(k); vals.append(1.0); continue
            a = 1 / dr**2; b = (-1 if stokes else 1) / (2 * r[i] * dr)
            for di, wt in ((1, a + b), (-1, a - b)):
                ii = i + di
                if ii < 0:
                    if stokes: continue
                    ii = 0
                rows.append(k); cols.append(I(ii, j)); vals.append(wt)
            rows += [k, k, k]; cols += [I(i, j + 1), I(i, j - 1), k]; vals += [1 / dz**2, 1 / dz**2, -2 / dr**2 - 2 / dz**2]
    return spl.splu(sp.csc_matrix((vals, (rows, cols)), shape=(N, N)))
LU_s, LU_p = build(True), build(False)
def solve(LU, rhs):
    rhs = rhs.copy(); rhs[-1, :] = 0; rhs[:, 0] = 0; rhs[:, -1] = 0; return LU.solve(rhs.ravel()).reshape(Nr, Nz)
gr = lambda f: np.gradient(f, dr, axis=0); gz = lambda f: np.gradient(f, dz, axis=1)
def run(kappa, nu, s=2.0, R0=6.0, z0=-8.0, Gam=3.0, T=120.0, label=""):
    om = Gam / np.pi * np.exp(-((RR - R0)**2 + (ZZ - z0)**2)); t = 0.0; hist = []
    dt = min(0.2 * min(dr, dz)**2 / max(nu, 1e-9), 0.1)
    while t <= T:
        w = np.clip(om, 0, None); W = w.sum(); rc, zc = (w * RR).sum() / W, (w * ZZ).sum() / W
        psi = solve(LU_s, -RR * om)
        if kappa > 0: phi = solve(LU_p, -kappa * np.exp(-((RR - rc)**2 + (ZZ - zc)**2) / s**2))
        else: phi = np.zeros_like(om)
        ur = -gz(psi) / RR + gr(phi); uz = gr(psi) / RR + gz(phi)
        div = gr(RR * ur) / RR + gz(uz)
        if len(hist) == 0 or t >= hist[-1]["t"] + 2.0:
            G = om.sum() * dr * dz; a2 = (w * ((RR - rc)**2 + (ZZ - zc)**2)).sum() / W
            hist.append(dict(t=t, Gam=float(G), core=float(np.sqrt(a2)), R=float(rc), z=float(zc), ommax=float(om.max())))
        # upwind advection + source terms + diffusion (explicit)
        dom_r = np.where(ur > 0, (om - np.roll(om, 1, 0)) / dr, (np.roll(om, -1, 0) - om) / dr)
        dom_z = np.where(uz > 0, (om - np.roll(om, 1, 1)) / dz, (np.roll(om, -1, 1) - om) / dz)
        lap = (np.roll(om, -1, 0) - 2 * om + np.roll(om, 1, 0)) / dr**2 + gr(om) / RR + (np.roll(om, -1, 1) - 2 * om + np.roll(om, 1, 1)) / dz**2
        rhs = -(ur * dom_r + uz * dom_z) + om * ur / RR - om * div + nu * (lap - om / RR**2)
        om = om + dt * rhs; om[0, :] = 0; om[-1, :] = 0; om[:, 0] = 0; om[:, -1] = 0
        t += dt
    return hist, om
out = {}
for kappa, label in ((0.0, "no sink"), (0.1, "sink kappa=0.1, s=2"), (0.4, "sink kappa=0.4, s=2")):
    nu = 0.025
    h, om = run(kappa, nu, label=label); out[label] = h
    print(f"{label:24s}: Gam {h[0]['Gam']:.2f} -> {h[-1]['Gam']:.2f}, core {h[0]['core']:.2f} -> {h[-1]['core']:.2f}, R {h[0]['R']:.2f} -> {h[-1]['R']:.2f}, "
          f"z {h[0]['z']:.1f} -> {h[-1]['z']:.1f}, peak vorticity {h[0]['ommax']:.3f} -> {h[-1]['ommax']:.3f}", flush=True)
json.dump(out, open("burgers/evolve.json", "w"))
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
for k, h in out.items():
    t = [x["t"] for x in h]; ax[0].plot(t, [x["Gam"] for x in h], label=k); ax[1].plot(t, [x["core"] for x in h], label=k); ax[2].plot(t, [x["R"] for x in h], label=k)
for a_, ttl in zip(ax, ("circulation", "core radius (rms)", "ring radius")): a_.set_title(ttl); a_.set_xlabel("t (delta/c)"); a_.legend(fontsize=7)
plt.tight_layout(); plt.savefig("burgers/evolve.png", dpi=110)
