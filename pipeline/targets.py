"""Target definitions. Each target gives, inside the unit sphere, radial functions G(rho) and H(rho) for
chi = G(rho) sin^2(theta) and I = H(rho) sin^2(theta), separately for W and U (U may be zero).
Exterior: W continues as a dipole (chi ~ G(1) sin^2/rho, I ~ H(1) sin^2/rho); U is zero outside."""
import numpy as np
from dataclasses import dataclass, field
from scipy.special import spherical_jn as jn
from typing import Callable, Optional

@dataclass
class Target:
    name: str
    family: str
    GW: Callable; HW: Callable
    GU: Optional[Callable] = None; HU: Optional[Callable] = None
    eta: float = 1.0                       # charge screening 1/eta (Hall-like coupling); 1 = none
    params: dict = field(default_factory=dict)
    note: str = ""

g1 = lambda m, x: x * jn(1, m * x)

def spheromak(lam=np.pi, eps=0.0):
    """Beltrami W (I = lam chi); optional aligned U = eps W."""
    GW = lambda x: g1(lam, x); HW = lambda x: lam * g1(lam, x)
    GU = (lambda x: eps * g1(lam, x)) if eps else None; HU = (lambda x: eps * lam * g1(lam, x)) if eps else None
    return Target(f"spheromak_l{lam:.3f}_e{eps:g}", "Beltrami (single)", GW, HW, GU, HU, params=dict(lam=lam, eps=eps))

def hill(lam=3.6331, b=None, eps=1.0):
    """Swirling Hill W: chi = [rho j1(lam rho) + (b/lam^2) rho^2] sin^2, I = lam chi; matched b if None. U = eps W."""
    if b is None: b = -lam**3 * jn(0, lam) / 3
    G = lambda x: g1(lam, x) + b / lam**2 * x**2
    return Target(f"hill_l{lam:.3f}_b{b:.3f}", "Swirling Hill", G, lambda x: lam * G(x),
                  (lambda x: eps * G(x)), (lambda x: eps * lam * G(x)), params=dict(lam=lam, b=b, eps=eps))

def double_beltrami(eta, m1, m2):
    """Hall-like ring: W = W(m1) + r W(m2), U = a1 W(m1) + a2 r W(m2) (see working record s.33)."""
    s = 2 * eta * (m1 + m2) / (1 + eta); d = np.sqrt(s * s - 4 * eta * m1 * m2)
    a1 = (2 * eta * m1 - s) / d; a2 = (2 * eta * m2 - s) / d
    r = -(m1**2 * jn(0, m1)) / (m2**2 * jn(0, m2))
    GW = lambda x: g1(m1, x) + r * g1(m2, x); HW = lambda x: m1 * g1(m1, x) + r * m2 * g1(m2, x)
    GU = lambda x: a1 * g1(m1, x) + a2 * r * g1(m2, x); HU = lambda x: a1 * m1 * g1(m1, x) + a2 * r * m2 * g1(m2, x)
    return Target(f"hall_eta{eta:g}_{m1:.3f}_{m2:.3f}", "Hall-like (double Beltrami)", GW, HW, GU, HU, eta=eta,
                  params=dict(eta=eta, mu1=m1, mu2=m2, ratio=r, alpha1=a1, alpha2=a2))

def rigid_core(base: Target, Omega=1.0):
    """Viscously relaxed variant: keep W of `base`, replace U by rigid rotation U_phi = Omega r (no poloidal U)."""
    return Target(base.name + "_rigidU", base.family + " + rigid U", base.GW, base.HW,
                  lambda x: 0.0 * x, lambda x: Omega * x**2,  # I = r U_phi = Omega r^2 -> H(rho) sin^2 with r = rho sin
                  eta=base.eta, params={**base.params, "Omega": Omega}, note="I_U = Omega r^2 (rigid rotation)")

def registry():
    T = [spheromak(np.pi, 0.5), hill(), hill(eps=0.5)]
    for eta, m1, m2 in [(1.5, 3.8503, 8.8136), (1.5, 4.0723, 8.2376), (1.5, 4.7219, 13.4535), (3.0, 3.4055, 12.1875)]:
        T.append(double_beltrami(eta, m1, m2))
    T.append(rigid_core(double_beltrami(1.5, 4.7219, 13.4535)))
    T.append(rigid_core(spheromak(np.pi)))
    return {t.name: t for t in T}
