"""Existence of a stable window: some core radius a with n < a < a_E (subsonic and Carnot-stable), R(a)/a >= 2, U_t < 1.
The engine ratio C/nu_v can always be tuned to place a* there (balance W = D gives C/nu_v = D/(nu_v E))."""
import numpy as np
P0, SIG = 5 / 7, 2 * np.sqrt(2) / 3
def window(n, nu_i, J, p0=P0, sig=SIG):
    Gam = 2 * np.pi * n
    A = lambda a: nu_i * Gam**2 / (16 * np.pi) + p0 * np.pi * a**2 + 2 * np.pi * sig * a
    Req = lambda a: (3 * J**2 / (8 * np.pi**3 * nu_i * a**2 * A(a)))**0.25
    E = lambda a: (8 * np.pi / 3) * Req(a) * A(a)
    a = np.geomspace(1e-3, 1e3, 6000); aE = a[np.argmin(E(a))]; R = Req(a); Ut = J / (2 * np.pi**2 * nu_i * a**2 * R**2)
    ok = (a > n) & (a < aE) & (R / a >= 2) & (Ut < 1)
    if not ok.any(): return aE, None
    D = 2 * nu_i * Req(a) * Gam**2 / a**2; r = D / E(a)       # required C/nu_v
    return aE, (a[ok].min(), a[ok].max(), r[ok].min(), r[ok].max())
print("n, nu_i, J -> a_E, stable window in a (xi), required C/nu_v")
for n in (1/3, 1, 2, 3):
    for nu_i in (1, 10, 100, 1000):
        for J in (100, 1000, 1e4, 1e5):
            aE, w = window(n, nu_i, J)
            if w: print(f"  n={n:.2f} nu_i={nu_i:5g} J={J:7g}: a_E={aE:6.2f}, a in [{w[0]:.2f}, {w[1]:.2f}], C/nu_v in [{w[2]:.3g}, {w[3]:.3g}]")
