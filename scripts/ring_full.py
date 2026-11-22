import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe
from scipy.integrate import solve_ivp
from hill_q import setup
lam=3.6331; beta,G=setup(lam); G1=G(1.0); atil=18.86; A=0.222; rO=0.8469; xi=1/atil
def weighted(r,z):
    ra=np.abs(r)+1e-9; rho=np.hypot(ra,z); s2=(ra/rho)**2; ins=rho<1
    chi=np.where(ins,G(rho)*s2,G1*s2/rho); I=np.where(ins,lam*G(rho)*s2,lam*G1*s2/rho)
    return chi,I
def bulkpsi(r,z):
    ra=np.abs(r)+1e-9; m=np.clip(4*ra*rO/((ra+rO)**2+z**2),0,1-1e-12); k=np.sqrt(m)
    return xi*np.sqrt(ra*rO)*((2/k-k)*ellipk(m)-2/k*ellipe(m))
n=500; X=np.linspace(-1.7,1.7,n); Z=np.linspace(-1.7,1.7,n); XX,ZZ=np.meshgrid(X,Z); dx=X[1]-X[0]
sg=np.sign(XX); sg[sg==0]=1
chi,I=weighted(XX,ZZ); cz,cx=np.gradient(chi,dx,dx)
Vx=A*(-cz/np.abs(XX))*sg; Vz=A*(cx/np.abs(XX))*sg*sg  # poloidal F in meridional plane (mirrored)
Vz=A*sg*cx/np.abs(XX)
Vphi=A*I/np.abs(XX)
psi=bulkpsi(XX,ZZ); pz,px=np.gradient(psi,dx,dx)
ux=-pz/np.abs(XX)*sg; uz=px/np.abs(XX)*sg
d=np.hypot(np.abs(XX)-rO,ZZ)/xi; nb=d**2/(d**2+2)
mask=d<0.8; ux[mask]=np.nan; uz[mask]=np.nan
fig=plt.figure(figsize=(16,15))
def deco(ax):
    t=np.linspace(0,2*np.pi,300); ax.plot(np.cos(t),np.sin(t),"w--",lw=1)
    ax.contour(XX,ZZ,np.where(np.hypot(XX,ZZ)<1,chi,np.nan),levels=[G1],colors="m",linewidths=1.3)
    ax.plot([rO,-rO],[0,0],"o",color="orange",mec="k",ms=7); ax.set_aspect("equal"); ax.set_xlabel("x / a  (meridional plane through the axis)"); ax.set_ylabel("z / a")
ax=fig.add_subplot(2,2,1)
pc=ax.pcolormesh(XX,ZZ,np.log10(np.hypot(ux,uz)+1e-6),cmap="magma",shading="auto",vmin=-3,vmax=0); fig.colorbar(pc,ax=ax,label="log₁₀ |u|/c")
ax.contourf(XX,ZZ,nb,levels=[0,0.5],colors=["#2b6cff"],alpha=0.6)
ax.streamplot(X,Z,ux,uz,color="w",density=2.2,linewidth=0.7,arrowsize=0.8); deco(ax)
ax.set_title("(a) bulk flow through the whole ring\nvortex ring of circulation 2πζ; blue: depleted core (n < ½)")
ax=fig.add_subplot(2,2,2)
m_=np.nanmax(abs(Vphi)); pc=ax.pcolormesh(XX,ZZ,Vphi*sg,cmap="RdBu_r",vmin=-m_,vmax=m_,shading="auto"); fig.colorbar(pc,ax=ax,label="toroidal drift V_φ/c (out of / into page)")
ax.streamplot(X,Z,Vx,Vz,color="k",density=2.2,linewidth=0.7,arrowsize=0.8); deco(ax)
ax.set_title("(b) weighted drift V of the carrying population\nstreamlines: poloidal; colour: toroidal component")
# 3D views
def rhs_w(t,y):
    x,yy,z=y; r=np.hypot(x,yy); e=1e-6
    c=lambda rr,zz: weighted(rr,zz)[0]
    Fr=-(c(r,z+e)-c(r,z-e))/(2*e)/r; Fz=(c(r+e,z)-c(r-e,z))/(2*e)/r; Fp=weighted(r,z)[1]/r
    return [Fr*x/r-Fp*yy/r, Fr*yy/r+Fp*x/r, Fz]
def torus(ax):
    u,v=np.mgrid[0:2*np.pi:80j,0:2*np.pi:30j]; R0=rO; rr=0.15
    ax.plot_surface((R0+rr*np.cos(v))*np.cos(u),(R0+rr*np.cos(v))*np.sin(u),rr*np.sin(v),color="m",alpha=0.12,linewidth=0)
    t=np.linspace(0,2*np.pi,200); ax.plot(rO*np.cos(t),rO*np.sin(t),0*t,color="orange",lw=2)
ax3=fig.add_subplot(2,2,3,projection="3d"); torus(ax3)
for ph in np.linspace(0,2*np.pi,8,endpoint=False):
    for r0 in (0.25,0.5,1.1,1.4):
        t=np.linspace(-1.6,1.6,1)
        # bulk streamlines lie in meridional planes: trace contour psi=const in (r,z)
    pass
Rg=np.linspace(0.01,1.7,300); Zg=np.linspace(-1.7,1.7,600); RG,ZG=np.meshgrid(Rg,Zg); PS=bulkpsi(RG,ZG)
cs=plt.figure().add_subplot().contour(RG,ZG,PS,levels=np.geomspace(PS[PS>0].min()*50,PS.max()*0.9,9)); plt.close()
for ph in np.linspace(0,2*np.pi,10,endpoint=False):
    for seg in cs.allsegs:
        for p in seg:
            if len(p)<10: continue
            ax3.plot(p[:,0]*np.cos(ph),p[:,0]*np.sin(ph),p[:,1],color="C3",lw=0.6,alpha=0.8)
ax3.set_title("(c) bulk streamlines in 3D: loops threading the hole,\nwrapping the vortex core (orange); magenta: closed-tori band")
ax4=fig.add_subplot(2,2,4,projection="3d"); torus(ax4)
for r0,col,T in ((rO+0.04,"C0",30),(rO+0.10,"C1",30),(0.45,"C2",12),(0.2,"C4",12)):
    sol=solve_ivp(rhs_w,[0,T],[r0,0,0],max_step=0.02,rtol=1e-8); y=sol.y
    keep=np.hypot(y[0],np.hypot(y[1],y[2]))<1.7; ax4.plot(y[0][keep],y[1][keep],y[2][keep],color=col,lw=1)
    sol=solve_ivp(rhs_w,[0,-T],[r0,0,0],max_step=0.02,rtol=1e-8); y=sol.y
    keep=np.hypot(y[0],np.hypot(y[1],y[2]))<1.7; ax4.plot(y[0][keep],y[1][keep],y[2][keep],color=col,lw=1)
ax4.set_title("(d) weighted streamlines in 3D: blue/orange closed (1,1) on the band;\ngreen/purple open, spiralling up the hole and out as a dipole")
for a_ in (ax3,ax4):
    for f in (a_.set_xlim,a_.set_ylim,a_.set_zlim): f(-1.4,1.4)
    a_.set_box_aspect((1,1,1)); a_.view_init(elev=22,azim=35)
fig.suptitle("Quantized ring q = 1, δ = +1, Γ = 10³ (thin core): a = 18.9ξ, O-ring at 0.85a",fontsize=13)
plt.tight_layout(); plt.savefig("ring_full.png",dpi=115)
