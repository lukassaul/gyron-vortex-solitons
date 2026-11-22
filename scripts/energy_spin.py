"""Energy at fixed spin J = hbar/2 along the compressible q=1 ring family (one carrying population, exterior counterflow).
Units: rho_m c^2 xi^3 (energy), xi (length). Shape integrals on the unit sphere (a=1).
E(A,a) = nu a^3 e_w(A) + Gam a e_f(A) + p0 (4pi/3) a^3 + 4 pi sig a^2 ;  spin: nu a^4 j(A) = E/2  (hbar/2 = kappa/2)."""
import numpy as np
from scipy.interpolate import RegularGridInterpolator as RGI
exec(open("transonic.py").read().split("res={}")[0])
p0=5/7; sig=2*np.sqrt(2)/3
dV=2*np.pi*RR*dr*dz
def integrals(chi,k,A):
    m=inside
    cx=np.vstack([np.zeros((1,nz)),chi]); Rx=np.concatenate([[0],R])
    gr=np.gradient(cx,Rx,axis=0)[1:]; gz=np.gradient(chi,dz,axis=1)
    kk=np.where(np.isfinite(k),k,1.0)
    Vp2=(gr**2+gz**2)/(kk*RR)**2; Vt2=(lam*chi/RR)**2          # V_pol = |grad chi|/(k r), V_phi = L/r = lam chi/r
    ew=np.sum(np.where(m,0.5*kk*(Vp2+Vt2)+kk*np.log(kk)-kk+1,0)*dV)
    I=kk*lam*chi; j=np.sum(np.where(m,I,0)*dV)
    # field: magnetic interior from curl f, electric interior from grad(k-1); exterior dipole + monopole
    Ix=np.vstack([np.zeros((1,nz)),I]); Ir=np.gradient(Ix,Rx,axis=0)[1:]; Iz=np.gradient(I,dz,axis=1)
    Dstar=np.gradient(np.gradient(cx,Rx,axis=0),Rx,axis=0)[1:]-gr/RR+np.gradient(gz,dz,axis=1)
    deep=m&(RHO<0.95)
    efm=0.5*np.sum(np.where(deep,(Ir**2+Iz**2+Dstar**2)/RR**2,0)*dV)
    kn=kk-1; kr=np.gradient(np.vstack([kn[:1],kn]),Rx,axis=0)[1:]; kz=np.gradient(kn,dz,axis=1)
    efe=0.5*np.sum(np.where(deep,kr**2+kz**2,0)*dV)
    th=np.linspace(0.05,np.pi-0.05,200); pts=np.c_[0.93*np.sin(th),0.93*np.cos(th)]
    kb=RGI((R,Z),kn,bounds_error=False,fill_value=None)(pts); Q=4*np.pi*np.trapezoid(kb*0.5*np.sin(th),th)
    Ib=RGI((R,Z),I,bounds_error=False,fill_value=None)(pts); D=np.trapezoid(Ib*np.sin(th)**3,th)/np.trapezoid(np.sin(th)**5,th)
    ef=efm+efe+ (4*np.pi/3)*D**2 + Q**2/(8*np.pi)      # exterior dipole energy ~ (1/2)int|curl f|^2 = (4pi/3) D^2 at a=1 ; Coulomb Q^2/(8 pi)
    band=None
    iz=np.argmin(abs(Z)); c=chi[:,iz]; ins=RHO[:,iz]<0.999; band=(np.max(np.where(ins,c,-1))-A*G1)/(A*G1)
    return ew,ef,j,Q,band
rows=[]; chi=None; A=0.02
while A<0.5:
    out=run(A,chi)
    if out is None: break
    chi,k,M=out; ew,ef,j,Q,band=integrals(chi,k,A)
    rows.append((A,ew,ef,j,Q,band,np.nanmax(np.where(inside&(RHO<0.95),M,0)))); A*=1.08
rows=np.array(rows); np.save("energy_spin_rows.npy",rows)
print("  A      e_w      e_f      j        Q(a=1)   band   Mach")
for r in rows[::3]: print("  %.3f  %.4f  %.4f  %.4f  %+.4f  %+.3f  %.2f"%tuple(r))
print("choke beyond A =",A)
