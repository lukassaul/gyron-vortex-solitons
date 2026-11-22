"""Can hourglass (focusing) contacts keep a circulating ring of momentum collimated?
(A) Fixed-wall limit: closed polygonal orbit of N contacts on hourglass waists (saddles: meridional radius rho, concave;
    azimuthal radius r0, convex). Paraxial 2x2 maps per plane; waist orientation alternating (strong focusing).
(B) Equal-mass impulse chain (Silsbee): angle map theta' = lam*theta + noise, lam = g/(rho1 + rho2) (radii >0 convex, <0 concave).
(C) Persistence requirement vs thermal noise."""
import numpy as np
def drift(l): return np.array([[1, l], [0, 1.]])
def mirror(k): return np.array([[1, 0], [-2 * k, 1.]])   # thin mirror, focusing power 2k (k>0 focusing)
print("(A) fixed walls, N = 12 contacts per turn, incidence phi = pi/2 - pi/N")
N = 12; phi = np.pi / 2 - np.pi / N; l = 1.0
res = np.zeros((40, 40)); A = np.geomspace(0.05, 20, 40); B = np.geomspace(0.05, 20, 40)
for i, rho in enumerate(A):          # meridional (concave) radius / l
    for j, r0 in enumerate(B):       # azimuthal (convex) radius / l
        ok = True
        for plane in (0, 1):
            M = np.eye(2)
            for n in range(N):
                alt = n % 2               # alternate waist orientation
                kt = (1 / rho if (plane ^ alt) == 0 else -1 / r0)
                k = kt / np.cos(phi) if plane == 0 else kt * np.cos(phi)   # tangential / sagittal astigmatism
                M = mirror(k) @ drift(l) @ M
            ok &= abs(np.trace(M)) < 2
        res[i, j] = ok
print(f"  stable fraction of (rho/l, r0/l) in [0.05, 20]^2: {res.mean():.2f}")
ii, jj = np.nonzero(res)
if len(ii): print(f"  stable rho/l in [{A[ii].min():.2f}, {A[ii].max():.2f}], r0/l in [{B[jj].min():.2f}, {B[jj].max():.2f}]")
print("\n(B) equal-mass impulse chain, lam = g/(rho1+rho2); convex nose rho1 = 1, gap g = 1")
for rho2 in (1.0, 3.0, -1.5, -2.0, -3.0, -10.0, -1e6):
    lam = 1.0 / (1.0 + rho2)
    print(f"  target face radius {rho2:+g}: lam = {lam:+.3f} -> {'focusing' if abs(lam) < 1 else 'defocusing'}")
print("\n(C) thermal noise: theta_{n+1} = lam theta_n + s xi;  rms theta = s/sqrt(1-lam^2); loss per link ~ theta^2")
for s in (1.0, 0.1, 1e-3):
    for lam in (-0.5, -0.1):
        th = s / np.sqrt(1 - lam**2); print(f"  noise s={s:g}, lam={lam}: rms angle {th:.3g}, links before momentum leaks out ~ {1/th**2:.3g}")
turns_per_s = 2.998e8 / (2 * np.pi * 3.86e-13); need = 2e-18 / turns_per_s
print(f"\n  persistence: ring turn rate ~ {turns_per_s:.1e} s^-1 -> allowed loss per turn {need:.1e}; with N links/turn, rms angle <~ {np.sqrt(need/12):.1e} rad")
print("  dilute gas: lateral thermal displacement per gap ~ gap (free flight ~ l) -> s ~ 1")
