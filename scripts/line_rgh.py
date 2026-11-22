"""Rankine-U / Gold-Hoyle-W line vortex: closed forms and the stable core radius.
Units: c=1, xi; energy per length in rho_m c^2 xi^2. s = sigma_w/K0 (charge-mode weight)."""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
p0=5/7; sig=2*np.sqrt(2)/3
def Wphi(r,W0,r0): x=r/r0; return W0*x/(1+x*x)
def energies(rc,ell,nu,Gam,s,W0,r0,R=1e4):
    Om=2*ell/(nu*np.pi*rc**4)
    EU=nu*np.pi*Om**2*rc**4/4                                   # = ell^2/(nu pi rc^4)
    Ein=quad(lambda r:(2*Om*Wphi(r,W0,r0))**2*2*np.pi*r,0,rc)[0]
    Eout=2*np.pi*(2*Om*rc*Wphi(rc,W0,r0))**2*np.log(R/rc)
    Eel=0.5*Gam*s**2*(Ein+Eout)
    Ebub=p0*np.pi*rc**2+2*np.pi*sig*rc
    return EU+Eel+Ebub, Om
def lam_q(rc,ell,nu,s,W0,r0):
    Om=2*ell/(nu*np.pi*rc**4); return -4*np.pi*s*Om*rc*Wphi(rc,W0,r0)
print("closed forms (inside the core, x = r/r0):")
print("  a(r) = Ω² r²/2 ;  b(r) = Ω W0 r0 ln(1+x²) ;  c(r) = W0² x²/(2(1+x²))  ;  pitch of W-lines = 2π r0 (uniform)")
print("  λ_q = -4π s Ω W0 r_c x_c/(1+x_c²) ;  ℓ = ν π Ω r_c⁴/2 ;  E_U = ℓ²/(ν π r_c⁴)")
print("  magnetic energy per length = π Γ s² W0² (independent of r0 and r_c); W kinetic = (ν/2) π r0² W0² ln(1+R²/r0²)\n")
print("   ell    nu   Gam   W0   r0  |  r_c*/xi   #minima   E''>0   lambda_q   r_c leading order (2ℓ²/(νπ²p0))^(1/6)")
for ell,nu,Gam,W0,r0 in ((1,10,1,0.2,1),(1,10,100,0.2,1),(10,10,1,0.2,1),(1,100,1,0.2,1),(1,10,1,0.5,0.3),(0.1,1,1,0.2,1)):
    s=1.0; rcs=np.geomspace(0.02,200,2000); E=np.array([energies(r,ell,nu,Gam,s,W0,r0)[0] for r in rcs])
    mins=[k for k in range(1,len(rcs)-1) if E[k]<E[k-1] and E[k]<E[k+1]]
    k=mins[0]; h=rcs[k]*1e-3; d2=(energies(rcs[k]+h,ell,nu,Gam,s,W0,r0)[0]-2*E[k]+energies(rcs[k]-h,ell,nu,Gam,s,W0,r0)[0])/h**2
    lo=(2*ell**2/(nu*np.pi**2*p0))**(1/6)
    print(f"  {ell:5g} {nu:4g} {Gam:5g} {W0:4g} {r0:4g} |  {rcs[k]:7.3f}     {len(mins)}       {d2>0}    {lam_q(rcs[k],ell,nu,s,W0,r0):+.4f}    {lo:.3f}")
