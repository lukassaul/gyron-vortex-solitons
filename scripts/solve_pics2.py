import numpy as np
from scipy.optimize import brentq
m0,q0,cE,rO,bw,hel=np.load("coef.npy"); alpha=1/137.036
p0=5/7; sig=2*np.sqrt(2)/3; ac=1.615
Er={ "B": (lambda a: 2*np.pi**2*rO*a*(np.log(8*rO*a)-ac), lambda a: 2*np.pi**2*rO*(np.log(8*rO*a)-ac+1)),
     "A": (lambda a: p0*4*np.pi/3*a**3+sig*4*np.pi*a**2, lambda a: p0*4*np.pi*a**2+sig*8*np.pi*a)}
def a_of_k(pic,k):
    E,dE=Er[pic]; S=4*np.pi*alpha*k
    f=lambda a: dE(a)-cE*S/a**2
    lo=0.5 if pic=="B" else 1e-6
    while f(lo)>0: lo/=2
    return brentq(f,lo,1e12)
def total(pic,k,Gam):
    E,dE=Er[pic]; a=a_of_k(pic,k); S=4*np.pi*alpha*k
    return (m0/q0)*np.sqrt(Gam*S)+cE*S/a+E(a), a
def solve(pic,Gam):
    ks=np.logspace(-4,12,800); f=[total(pic,k,Gam)[0]-k for k in ks]
    for i in range(len(ks)-1):
        if f[i]>0 and f[i+1]<0:
            k=brentq(lambda k: total(pic,k,Gam)[0]-k,ks[i],ks[i+1]); break
    else: return None
    Etot,a=total(pic,k,Gam); S=4*np.pi*alpha*k; Q=np.sqrt(S/Gam)
    EM=(m0/q0)*np.sqrt(Gam*S); EE=cE*S/a; ER=Er[pic][0](a)
    return dict(a=a,kappa=k,Mach=np.sqrt(Q/(q0*a)),WB_WE=EM/EE,f=(EM/k,EE/k,ER/k),d2=Er[pic][1](a))
for pic in ("B","A"):
    print(f"\n=== picture {pic} ({'thin core, line tension' if pic=='B' else 'bubble, pressure+surface'}) ===")
    print(" Gamma     a/xi    r_O/xi   b/xi      kappa    Mach   W_B/W_E   M / E / restore")
    for Gam in (1e-4,1e-2,1e-1,1,10,1e2,1e4):
        d=solve(pic,Gam)
        if d is None: print(f"{Gam:7.0e}  no self-consistent solution"); continue
        print(f"{Gam:7.0e} {d['a']:8.3f} {rO*d['a']:8.3f} {bw*d['a']:7.3f} {d['kappa']:10.4g} {d['Mach']:7.3g} {d['WB_WE']:8.3g}   {d['f'][0]:.2f}/{d['f'][1]:.2f}/{d['f'][2]:.2f}")
