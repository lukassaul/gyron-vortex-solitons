import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from hill_q import q_profile, setup
qmid=lambda lam: dict(q_profile(lam)[3])[0.5]-1
lam1=brentq(qmid,3.5,3.8,xtol=1e-4)
beta,G=setup(lam1); rO,chiO,chiB,qs=q_profile(lam1)
print(f"q=1 at lambda a={lam1:.4f}, beta={beta:.4f}, rO={rO:.4f}, q-profile={qs}")
def fields(r,z,s):   # s=+1 / -1 handedness
    rho=np.hypot(r,z)+1e-12; s2=(r/rho)**2
    gin=G(rho); gB=G(1.0); D=1.0*gB*lam1  # I_ext coefficient: I continuous -> D = lam*G(1)
    chi=np.where(rho<1,gin*s2,gB*s2/rho)
    I=np.where(rho<1,lam1*gin*s2,D*s2/rho)
    return chi,s*I
fig=plt.figure(figsize=(15,9))
R=np.linspace(1e-3,1.8,400); Z=np.linspace(-1.8,1.8,800); RR,ZZ=np.meshgrid(R,Z)
for j,s in enumerate((+1,-1)):
    chi,I=fields(RR,ZZ,s); Fp=I/RR
    ax=fig.add_subplot(2,2,1+2*j)
    m=abs(Fp).max(); pc=ax.pcolormesh(RR,ZZ,Fp,cmap="RdBu_r",vmin=-m,vmax=m,shading="auto")
    ax.contour(RR,ZZ,chi,levels=np.linspace(0,chiO,24)[1:],colors="k",linewidths=0.7)
    ax.contour(RR,ZZ,chi,levels=[chiB],colors="m",linewidths=1.6)
    t=np.linspace(-np.pi/2,np.pi/2,200); ax.plot(np.cos(t),np.sin(t),"g--")
    ax.plot(rO,0,"o",color="orange",mec="k",ms=8); ax.set_aspect("equal")
    ax.set_title(f"{'right' if s>0 else 'left'}-handed: χ contours, colour F_φ\nmagenta = edge of closed tori (q = {'+' if s>0 else '−'}1)")
    ax3=fig.add_subplot(2,2,2+2*j,projection="3d")
    def rhs(t,y):
        x,yy,z=y; r=np.hypot(x,yy); e=1e-6
        c=lambda rr,zz: fields(rr,zz,s)[0]
        Fr=-(c(r,z+e)-c(r,z-e))/(2*e)/r; Fz=(c(r+e,z)-c(r-e,z))/(2*e)/r
        Fphi=fields(r,z,s)[1]/r
        return [Fr*x/r-Fphi*yy/r, Fr*yy/r+Fphi*x/r, Fz]
    for k,r0 in enumerate((rO+0.03,rO+0.07)):
        sol=solve_ivp(rhs,[0,40],[r0,0,0],max_step=0.02,rtol=1e-8)
        ax3.plot(*sol.y,lw=1.1,color=["C0","C1"][k])
    ax3.set_title("closed F-lines: one toroidal turn per poloidal turn")
    for f in (ax3.set_xlim,ax3.set_ylim,ax3.set_zlim): f(-1,1)
    ax3.set_box_aspect((1,1,1))
plt.tight_layout(); plt.savefig("hill_q1_both.png",dpi=130)
