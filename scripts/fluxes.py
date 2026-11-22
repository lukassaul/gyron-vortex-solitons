import numpy as np
from scipy.special import spherical_jn as jn
from scipy.integrate import dblquad, solve_ivp
from scipy.optimize import brentq
lam=np.pi
g=lambda rho: rho*jn(1,lam*rho)
xO=brentq(lambda x: x*jn(0,x)-jn(1,x),2,3.1); rO=xO/lam
chiO=g(rO); chiB=g(1.0)             # O-point and boundary (theta=pi/2) values
Psi_closed=2*np.pi*(chiO-chiB)      # poloidal flux in closed-torus band
Psi_tot=2*np.pi*chiO                # poloidal flux through O-ring (axis to O)
# toroidal flux through meridional half-plane inside core: int I/r dr dz
Phi_t=dblquad(lambda th,rho: lam*g(rho)*np.sin(th)**2/(rho*np.sin(th))*rho, 0,1, 0,np.pi)[0]
print(f"chi_O={chiO:.4f} chi_boundary={chiB:.4f}")
print(f"Psi(axis->O)={Psi_tot:.4f}  Psi(closed band)={Psi_closed:.4f}  Phi_t(core)={Phi_t:.4f}")
print(f"Phi_t/Psi = {Phi_t/Psi_tot:.4f};  h_total=26.319 -> h/(2 Phi_t Psi) = {26.319/(2*Phi_t*Psi_tot):.4f}")
# rotational transform q = toroidal turns per poloidal turn on a few tori
def rhs(t,y):
    r,z,ph=y; rho=np.hypot(r,z); e=1e-6
    c=lambda rr,zz: g(np.hypot(rr,zz))*(rr**2/(rr**2+zz**2))
    Fr=-(c(r,z+e)-c(r,z-e))/(2*e)/r; Fz=(c(r+e,z)-c(r-e,z))/(2*e)/r
    Fp=lam*c(r,z)/r
    return [Fr,Fz,Fp/r]
for r0 in np.linspace(rO+0.01,0.99,6):
    ev=lambda t,y: y[1]; ev.direction=1
    sol=solve_ivp(rhs,[0,200],[r0,-1e-9,0],events=ev,rtol=1e-9,max_step=0.01)
    ts=sol.t_events[0]; 
    if len(ts)>2:
        dphi=np.diff(sol.y_events[0][:,2]).mean()
        print(f"r0={r0:.3f}: q = dphi/2pi per poloidal turn = {dphi/(2*np.pi):.4f}")
