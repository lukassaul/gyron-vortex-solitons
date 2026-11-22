"""Thin ring from the Rankine/Gold-Hoyle line, with conserved J (toroidal U), Gamma_U (Kelvin circulation of U around the core),
Phi_W (W-line flux) and closure N (twist length r0 = R/N). Units: c=1, xi; energy rho_m c^2 xi^3. Electric term omitted (small here)."""
import numpy as np
from scipy.optimize import minimize_scalar
p0=5/7; sig=2*np.sqrt(2)/3
def parts(R,rc,J,GU,PhiW,N,nu,Gam,s=1.0):
    L=2*np.pi*R; lnN=np.log(1+N**2)
    Upol=nu*GU**2/(16*np.pi)                         # Rankine core at fixed circulation: per length, independent of rc
    Utor=J**2/(8*np.pi**3*nu*rc**2*R**4)             # toroidal flow needed to hold J: U_z = J/(2 pi^2 nu rc^2 R^2)
    W0=PhiW*N**2/(np.pi*R**2*lnN)                    # Gold-Hoyle amplitude from the flux quantum, r0 = R/N
    Wkin=nu*PhiW**2*N**2/(2*np.pi*R**2*lnN)
    Wmag=np.pi*Gam*s**2*W0**2
    bub=p0*np.pi*rc**2+2*np.pi*sig*rc
    return {k:L*v for k,v in dict(Upol=Upol,Utor=Utor,Wkin=Wkin,Wmag=Wmag,bubble=bub).items()}
def E_R(R,*a):
    f=lambda lrc: sum(parts(R,np.exp(lrc),*a).values())
    r=minimize_scalar(f,bounds=(np.log(1e-3),np.log(R*0.5)),method="bounded"); return r.fun,np.exp(r.x)
def solve(GU,PhiW,N,nu,Gam):
    kappa=10.0
    for _ in range(60):                               # self-consistent spin: J = hbar/2 = kappa/2 (units rho_m c xi^4)
        J=kappa/2; Rs=np.geomspace(0.5,500,600); E=np.array([E_R(R,J,GU,PhiW,N,nu,Gam)[0] for R in Rs])
        mins=[k for k in range(1,len(Rs)-1) if E[k]<E[k-1] and E[k]<E[k+1]]
        if not mins: return None
        k=mins[0]; kn=E[k]
        if abs(kn-kappa)<1e-6*kappa: break
        kappa=0.5*kappa+0.5*kn
    R=Rs[k]; E0,rc=E_R(R,J,GU,PhiW,N,nu,Gam)
    return R,rc,kappa,len(mins),parts(R,rc,J,GU,PhiW,N,nu,Gam)
print(" Γ_U   Φ_W  N   ν   Γ  |  R/ξ     r_c/ξ   r_c/R   κ      #min  | shortening (Upol+bubble)  vs  backreaction (Utor+Wkin+Wmag)")
for GU,PhiW,N,nu,Gam in ((2*np.pi,1,1,10,1),(2*np.pi,1,1,10,100),(2*np.pi,1,2,10,1),(2*np.pi,3,1,10,1),(2*np.pi,1,1,100,1),(np.pi,1,1,10,1)):
    out=solve(GU,PhiW,N,nu,Gam)
    if out is None: print(f" {GU:4.2f} {PhiW:4g} {N:2d} {nu:4g} {Gam:4g} | no minimum"); continue
    R,rc,k,nm,P=out; sh=P["Upol"]+P["bubble"]; br=P["Utor"]+P["Wkin"]+P["Wmag"]
    print(f" {GU:4.2f} {PhiW:4g} {N:2d} {nu:4g} {Gam:4g} | {R:7.2f} {rc:7.3f} {rc/R:7.3f} {k:9.3g} {nm:3d}  |  {sh:9.3g}  vs  {br:9.3g}   (Utor {P['Utor']:.3g}, Wkin {P['Wkin']:.3g}, Wmag {P['Wmag']:.3g})")
