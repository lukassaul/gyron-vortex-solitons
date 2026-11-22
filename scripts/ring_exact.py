import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import spherical_jn as jn
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

a=1.0; lam=np.pi/a; A=1.0
def fields(r,z):
    rho=np.sqrt(r**2+z**2)+1e-12; s2=(r/rho)**2
    inside=rho<a
    chi=np.where(inside, A*rho*jn(1,lam*rho)*s2, (A*a/np.pi)*(a/rho)*s2)
    I  =np.where(inside, lam*A*rho*jn(1,lam*rho)*s2, A*a*s2/rho)
    return chi,I
# O-ring: d/dx (x j1) = 0  ->  x j0 = j1
xO=brentq(lambda x: x*jn(0,x)-jn(1,x),2,3.1); rO=xO/lam
print("lambda a =",lam*a," O-ring radius r_O/a =",rO)

n=600; R=np.linspace(1e-3,2.2,n); Z=np.linspace(-2.2,2.2,2*n)
RR,ZZ=np.meshgrid(R,Z); chi,I=fields(RR,ZZ)
Fphi=I/RR

fig,ax=plt.subplots(1,3,figsize=(16,6.2))
# (a) poloidal weighted flux + toroidal swirl
pc=ax[0].pcolormesh(RR,ZZ,Fphi,cmap="RdBu_r",shading="auto",vmin=-abs(Fphi).max(),vmax=abs(Fphi).max())
lev=np.linspace(chi.min(),chi.max(),26)[1:]
ax[0].contour(RR,ZZ,chi,levels=lev,colors="k",linewidths=0.8)
ax[0].set_title("(a) poloidal lines of F (contours of χ)\ncolour: toroidal flux F_φ = I/r")
fig.colorbar(pc,ax=ax[0],fraction=0.046)
# (b) poloidal B lines = contours of I
levI=np.linspace(0,I.max(),22)[1:]
ax[1].contour(RR,ZZ,I,levels=levI,colors="C3",linewidths=0.8)
ax[1].set_title("(b) poloidal B lines (contours of I)\ninside B = (g_Φ/c²)λF; outside dipole")
for k in (0,1):
    t=np.linspace(-np.pi/2,np.pi/2,200)
    ax[k].plot(a*np.cos(t),a*np.sin(t),"g--",lw=1.5)
    ax[k].plot([rO],[0],"o",color="orange",ms=9,mec="k")
    ax[k].set_aspect("equal"); ax[k].set_xlabel("r / a"); ax[k].set_ylabel("z / a")
    ax[k].set_xlim(0,2.2); ax[k].set_ylim(-2.2,2.2)
ax[0].text(rO+0.05,0.08,"O-ring\n(bulk vortex core)",fontsize=8)

# (c) 3D helical F-lines on nested tori (inside)
ax3=fig.add_subplot(1,3,3,projection="3d"); ax[2].remove()
def rhs(t,y):
    x,yy,z=y; r=np.hypot(x,yy); rho=np.sqrt(r*r+z*z)
    # F = grad chi x grad phi + I grad phi  (inside)
    e=1e-6
    c=lambda rr,zz: fields(np.array(rr),np.array(zz))[0]
    dchr=(c(r+e,z)-c(r-e,z))/(2*e); dchz=(c(r,z+e)-c(r,z-e))/(2*e)
    Fr=-dchz/r; Fz=dchr/r; Fp=fields(np.array(r),np.array(z))[1]/r
    cp,sp=x/r,yy/r
    return [Fr*cp-Fp*sp, Fr*sp+Fp*cp, Fz]
cols=["C0","C1","C2"]
for i,r0 in enumerate([rO+0.05,rO+0.15,rO+0.3]):
    sol=solve_ivp(rhs,[0,60],[r0,0,0],max_step=0.02,rtol=1e-8)
    ax3.plot(*sol.y,color=cols[i],lw=0.9)
u,v=np.mgrid[0:2*np.pi:60j,0:np.pi:30j]
ax3.plot_wireframe(a*np.cos(u)*np.sin(v),a*np.sin(u)*np.sin(v),a*np.cos(v),color="g",alpha=0.08)
ax3.set_box_aspect((1,1,1)); ax3.set_title("(c) helical F-lines on nested tori")
for f in (ax3.set_xlim,ax3.set_ylim,ax3.set_zlim): f(-1,1)
plt.tight_layout(); plt.savefig("ring_flow.png",dpi=150)
