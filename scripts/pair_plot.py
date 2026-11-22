import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import spherical_jn as jn
curve=np.load("q1curve.npy")
def member(l): i=np.argmin(abs(curve[:,0]-l)); return curve[i,0],curve[i,1],curve[i,2]
l1,b1,G1=member(5.4); l2,b2,G2=member(5.2); L2=-l2
A1=1/4.03; A2=-A1*G1/G2      # A1 at peak |V| = c
def fields(lam,b,A,G1_,XX,ZZ):
    ra=np.abs(XX)+1e-9; rho=np.hypot(ra,ZZ); s2=(ra/rho)**2; ins=rho<1
    G=lambda r: r*jn(1,abs(lam)*r)+b/lam**2*r**2
    chi=A*np.where(ins,G(rho)*s2,G1_*s2/rho); I=lam*chi
    return chi,I
n=500; X=np.linspace(-1.6,1.6,n); Z=X.copy(); XX,ZZ=np.meshgrid(X,Z); dx=X[1]-X[0]; sg=np.sign(XX)
def flow(chi):
    cz,cx=np.gradient(chi,dx,dx); return -cz/np.abs(XX)*sg, sg*cx/np.abs(XX)
c1,I1=fields(l1,b1,A1,G1,XX,ZZ); c2,I2=fields(L2,b2,A2,G2,XX,ZZ)
fig,ax=plt.subplots(1,3,figsize=(18,6.4))
for a_,(chi,I,tt) in zip(ax,[(c1,I1,f"(a) population + : λa = {l1:.1f} (right-handed), amplitude {A1:.3f}"),
                            (c2,I2,f"(b) population − : λa = −{l2:.1f} (left-handed), counter-streaming"),
                            (c1-c2,I1-I2,"(c) net flux F = F₊ − F₋ (sets E and B);\nmass flux F₊ + F₋ has no exterior dipole: zero impulse")]):
    u,w=flow(chi); Vp=I/np.abs(XX)*sg; m=np.nanmax(abs(Vp))
    pc=a_.pcolormesh(XX,ZZ,Vp,cmap="RdBu_r",vmin=-m,vmax=m,shading="auto"); fig.colorbar(pc,ax=a_,fraction=0.046,label="toroidal component / c")
    a_.streamplot(X,Z,u,w,color="k",density=1.8,linewidth=0.6,arrowsize=0.7)
    t=np.linspace(0,2*np.pi,200); a_.plot(np.cos(t),np.sin(t),"g--",lw=1); a_.set_aspect("equal"); a_.set_title(tt,fontsize=10)
    a_.set_xlabel("x / a"); a_.set_ylabel("z / a")
fig.suptitle("Counter-streaming pair of closed q = 1 rings (no bulk circulation): zero impulse, net charge and small net moment",fontsize=12)
plt.tight_layout(); plt.savefig("pair_ring.png",dpi=115)
