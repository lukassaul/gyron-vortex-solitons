"""Stability of the Carnot-balanced core size (thin ring, fixed circulation Gam and spin J; units xi, c, rho_m).
Fast: ring radius R at mechanical equilibrium for each core radius a (dE/dR = 0).  Slow: a evolves by the budget
  E_eff'(a) da/dt = W - D,   W = C E_eff (engine, C = (T0/Tpar - 1) eps kappa1),   D = 2 nu_i nu_v R Gam^2 / a^2.
Fixed point a*: W = D. Stable iff (W - D)'/E_eff' < 0 at a*; for D falling with a this means E_eff'(a*) < 0,
i.e. a* below the passive energy-minimum core radius a_E."""
import numpy as np
from scipy.optimize import brentq
P0, SIG = 5 / 7, 2 * np.sqrt(2) / 3
nu_i, nu_v, Gam, J = 10.0, 1e-3, 2 * np.pi, 300.0         # inertia ratio, viscosity (c xi units), circulation quantum, spin
def A(a): return nu_i * Gam**2 / (16 * np.pi) + P0 * np.pi * a**2 + 2 * np.pi * SIG * a      # tension per length
def Req(a): return (3 * J**2 / (8 * np.pi**3 * nu_i * a**2 * A(a)))**0.25
def Eeff(a): return (8 * np.pi / 3) * Req(a) * A(a)
def D(a): return 2 * nu_i * nu_v * Req(a) * Gam**2 / a**2
a = np.geomspace(0.02, 5, 4000); E = Eeff(a); dE = np.gradient(E, a)
aE = a[np.argmin(E)]
print(f"passive energy minimum: a_E = {aE:.3f} xi, R = {Req(aE):.3f} xi (R/a = {Req(aE)/aE:.2f})")
print(f"{'C (1/time)':>12} {'a*':>8} {'R*':>8} {'R*/a*':>7} {'E_eff`(a*)':>11} {'growth rate':>12}  verdict")
for C in [1e-3, 3e-3, 1e-2, 3e-2, 1e-1, 3e-1, 1.0]:
    f = lambda x: C * Eeff(x) - D(x)
    s = np.sign(f(a)); idx = np.nonzero(np.diff(s))[0]
    if not len(idx): print(f"{C:12.1e}   no balance in range"); continue
    for i in idx:
        x = brentq(f, a[i], a[i + 1]); h = 1e-5 * x
        dEx = (Eeff(x + h) - Eeff(x - h)) / (2 * h); dfx = (f(x + h) - f(x - h)) / (2 * h)
        lam = dfx / dEx                     # d(da/dt)/da at the fixed point
        print(f"{C:12.1e} {x:8.3f} {Req(x):8.3f} {Req(x)/x:7.2f} {dEx:+11.3e} {lam:+12.3e}  {'STABLE' if lam < 0 else 'unstable'}")
print(f"\nthreshold for a stable balance: C > C_min = D(a_E)/E_eff(a_E) = {D(aE)/Eeff(aE):.3e}")
