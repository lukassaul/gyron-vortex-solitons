import numpy as np
from scipy.special import spherical_jn as jn
from scipy.integrate import quad
from hill_q import setup
# net k inside at O(M^2): k = (delta/c^2)[beta*chi - |F|^2/2], populations F± = ±(1±delta)F/2
# exterior Laplace with k continuous at a=1 -> Q = 4*pi*a*<k(a,theta)>
def Q0(lam):
    beta,G=setup(lam); e=1e-6
    dG=lambda r:(G(r+e)-G(r-e))/(2*e)
    def integrand(th, rho=1-1e-9):
        s=np.sin(th); c=np.cos(th); r=rho*s
        chi=G(rho)*s*s
        chi_r=dG(rho)*s*s; chi_t=G(rho)*2*s*c
        Fpol2=(chi_r**2+(chi_t/rho)**2)/r**2
        Fphi2=(lam*chi/r)**2
        return (beta*chi-0.5*(Fpol2+Fphi2))*0.5*s   # <.> over sphere
    avg=quad(integrand,1e-6,np.pi-1e-6,limit=200)[0]
    return 4*np.pi*avg, beta
for lam in (np.pi,3.6331):
    for sgn in (+1,-1):
        q,b=Q0(abs(lam))   # handedness enters only via lam^2 and |F|^2
        print(f"lambda a={sgn*lam:+.4f} beta={b:+.4f}:  Q = delta * A^2/c^2 * ({q:+.4f})")
