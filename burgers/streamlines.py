"""Streamlines of the sink-driven Burgers ring (meridional plane), and a steadiness test.
Units delta = 1, c = 1. Ring radius R, circulation Gam, sink rate kappa (nu_v = kappa/4 so delta = 1).
Fields: psi from D*psi = -r om_phi (Gaussian toroidal vorticity); phi from lap phi = -kappa g (toroidal sink, g Gaussian
peak 1 at the core). Velocity u = (-psi_z/r, psi_r/r) + grad phi; comoving frame: subtract the measured translation.
Steadiness: residual of the toroidal vorticity equation  (u.grad) om - om u_r/r + om div u - nu (lap om - om/r^2) = 0."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R, Gam, kappa = 8.0, 3.0, 0.1; nu = kappa / 4
Lr, Lz, Nr, Nz = 24.0, 24.0, 241, 481
r = np.linspace(Lr / Nr / 2, Lr, Nr); z = np.linspace(-Lz, Lz, Nz); dr, dz = r[1] - r[0], z[1] - z[0]
RR, ZZ = np.meshgrid(r, z, indexing="ij"); d2 = (RR - R)**2 + ZZ**2
om = Gam / np.pi * np.exp(-d2); g = np.exp(-d2)
idx = lambda i, j: i * Nz + j
def op(stokes):
    rows, cols, vals = [], [], []
    for i in range(Nr):
        for j in range(Nz):
            k = idx(i, j)
            if i == Nr - 1 or j in (0, Nz - 1): rows.append(k); cols.append(k); vals.append(1.0); continue
            ri = r[i]; a = 1 / dr**2; b = (-1 if stokes else 1) / (2 * ri * dr)
            for di, w in ((1, a + b), (-1, a - b)):
                ii = i + di
                if ii < 0:                                   # axis
                    if stokes: continue                      # psi = 0 on the axis
                    ii = 0
                rows.append(k); cols.append(idx(ii, j)); vals.append(w)
            rows += [k, k, k]; cols += [idx(i, j + 1), idx(i, j - 1), k]; vals += [1 / dz**2, 1 / dz**2, -2 / dr**2 - 2 / dz**2]
    return sp.csr_matrix((vals, (rows, cols)), shape=(Nr * Nz, Nr * Nz))
def solve(A, rhs):
    rhs = rhs.copy(); rhs[-1, :] = 0; rhs[:, 0] = 0; rhs[:, -1] = 0
    return spl.spsolve(A.tocsc(), rhs.ravel()).reshape(Nr, Nz)
psi = solve(op(True), -RR * om); phi = solve(op(False), -kappa * g)
ur = -np.gradient(psi, dz, axis=1) / RR + np.gradient(phi, dr, axis=0)
uz = np.gradient(psi, dr, axis=0) / RR + np.gradient(phi, dz, axis=1)
w = om / om.sum(); Uself = float(np.sum(w * (np.gradient(psi, dr, axis=0) / RR)))      # vorticity-weighted swirl uz
Usaff = Gam / (4 * np.pi * R) * (np.log(8 * R) - 0.558)
uzc = uz - Uself
div = np.gradient(RR * ur, dr, axis=0) / RR + np.gradient(uz, dz, axis=1)
lap = lambda f: np.gradient(np.gradient(f, dr, axis=0), dr, axis=0) + np.gradient(f, dr, axis=0) / RR + np.gradient(np.gradient(f, dz, axis=1), dz, axis=1)
adv = ur * np.gradient(om, dr, axis=0) + uzc * np.gradient(om, dz, axis=1)
terms = dict(advection=adv, stretching=-om * ur / RR, dilatation=om * div, viscous=-nu * (lap(om) - om / RR**2))
res = sum(terms.values()); core = d2 < 9
scale = max(np.sqrt(np.mean(t[core]**2)) for t in terms.values())
print(f"translation: measured {Uself:.4f} c, Saffman Gaussian ring {Usaff:.4f} c")
print("rms over the core (d < 3 delta):", {k: f"{np.sqrt(np.mean(t[core]**2)):.2e}" for k, t in terms.items()})
print(f"residual rms / largest term = {np.sqrt(np.mean(res[core]**2))/scale:.2f}")
# Burgers balance check on its own: dilatation + viscous + (radial advection by the sink only)
ur_s, uz_s = np.gradient(phi, dr, axis=0), np.gradient(phi, dz, axis=1); divs = np.gradient(RR * ur_s, dr, axis=0) / RR + np.gradient(uz_s, dz, axis=1)
bal = ur_s * np.gradient(om, dr, axis=0) + uz_s * np.gradient(om, dz, axis=1) + om * divs - nu * (lap(om) - om / RR**2)
print(f"Burgers part alone (sink advection + dilatation + viscous): residual/scale = {np.sqrt(np.mean(bal[core]**2))/max(np.sqrt(np.mean((om*divs)[core]**2)),1e-30):.2f}")
fig, ax = plt.subplots(1, 3, figsize=(17, 6))
for a_, (UR, UZ, ttl) in zip(ax[:2], ((ur, uzc, "comoving frame: swirl + intake"), (ur_s, uz_s, "intake only (sink)"))):
    sp_ = np.hypot(UR, UZ)
    a_.streamplot(z, r, UZ, UR, density=2.2, color=np.log10(sp_ + 1e-6), cmap="viridis", linewidth=0.7)
    a_.contour(z, r, om, levels=[om.max() * 0.1, om.max() * 0.5], colors="r", linewidths=1)
    a_.set_xlim(-16, 16); a_.set_ylim(0, 16); a_.set_aspect("equal"); a_.set_xlabel("z / delta"); a_.set_ylabel("r / delta"); a_.set_title(ttl, fontsize=10)
im = ax[2].pcolormesh(z, r, res / scale, cmap="RdBu_r", vmin=-1, vmax=1, shading="auto"); ax[2].set_xlim(R - 6 - 8, 8 + 6 - 8 + 8) if False else ax[2].set_xlim(-6, 6)
ax[2].set_ylim(R - 6, R + 6); ax[2].set_aspect("equal"); ax[2].set_title("vorticity-equation residual / largest term", fontsize=10); plt.colorbar(im, ax=ax[2])
ax[2].set_xlabel("z / delta")
plt.tight_layout(); plt.savefig("burgers/streamlines.png", dpi=110)
