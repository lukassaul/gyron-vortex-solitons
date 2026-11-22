import numpy as np
from scipy.optimize import brentq
from hill_q import setup
from solve_pics2 import Er, m0, q0, cE, rO
lam=3.6331; beta,G=setup(lam); G1=G(1.0); alpha=1/137.036
print(f"exterior dipole coefficient G(1) = {G1:.4f}  -> weighted Lamb impulse P_w = 4*pi*G1*A*a^3 (per unit weighted inertia)")
C=np.pi*rO**2/(2*G1)       # impulse balance: nu*A*a = C
print(f"impulse balance:  nu * A * (a/xi) = pi r_O^2/(2 G1) = {C:.4f}")
# unknowns a, kappa, Gamma for given nu (picture B): E'(a)=0, kappa=E, nu*A*a=C with A^2 a = Q/q0, Q=sqrt(4 pi alpha kappa/Gamma)
E,dE=Er["B"]
def solve(nu):
    out=None
    def F(a):
        # from balance: cE*S/a^2 = dE(a), S=4 pi alpha kappa
        S=dE(a)*a**2/cE; kappa=S/(4*np.pi*alpha)
        A=C/(nu*a); Q=q0*A**2*a; Gam=S/Q**2
        EM=(m0/q0)*Gam*Q; Etot=EM+cE*S/a+E(a)
        return Etot-kappa, (a,kappa,Gam,A)
    as_=np.geomspace(0.8,1e6,4000); vals=[F(a)[0] for a in as_]
    sols=[]
    for i in range(len(as_)-1):
        if np.sign(vals[i])!=np.sign(vals[i+1]):
            a=brentq(lambda x:F(x)[0],as_[i],as_[i+1]); sols.append(F(a)[1])
    return sols
print("\n   nu        a/xi     r_O/xi     kappa      Gamma      A    peak|V|/c")
for nu in (1e-3,1e-2,0.1,0.3,1,3):
    for a,k,Gm,A in solve(nu):
        print(f"{nu:7.0e} {a:10.3f} {rO*a:9.3f} {k:10.4g} {Gm:10.4g} {A:7.3f} {3*A:7.2f}")
