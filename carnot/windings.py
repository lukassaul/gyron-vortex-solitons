"""Same Carnot balance for rings of different circulation (winding) n: Gam = 2 pi n."""
import numpy as np
from scipy.optimize import brentq
P0, SIG = 5 / 7, 2 * np.sqrt(2) / 3; nu_i, nu_v, J, C = 10.0, 1e-3, 300.0, 0.03
for n in (1/3, 1, 2, 3, 5):
    Gam = 2 * np.pi * n
    A = lambda a: nu_i * Gam**2 / (16 * np.pi) + P0 * np.pi * a**2 + 2 * np.pi * SIG * a
    Req = lambda a: (3 * J**2 / (8 * np.pi**3 * nu_i * a**2 * A(a)))**0.25
    E = lambda a: (8 * np.pi / 3) * Req(a) * A(a); D = lambda a: 2 * nu_i * nu_v * Req(a) * Gam**2 / a**2
    a = np.geomspace(0.01, 20, 4000); aE = a[np.argmin(E(a))]; f = lambda x: C * E(x) - D(x)
    idx = np.nonzero(np.diff(np.sign(f(a))))[0]
    out = []
    for i in idx:
        x = brentq(f, a[i], a[i + 1]); h = 1e-5 * x; lam = ((f(x+h)-f(x-h))/(E(x+h)-E(x-h)))
        out.append(f"a*={x:.3f} R*={Req(x):.2f} E={E(x):.1f} {'stable' if lam < 0 else 'unstable'}")
    print(f"n={n:5.3f}: a_E={aE:.3f}  " + ("; ".join(out) if out else "no balance"))
