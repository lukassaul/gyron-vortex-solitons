import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import spherical_jn as jn
exec(open("double_beltrami.py").read().split("results={}")[0])
sols=[(1.5,3.8503,8.8136),(3.0,3.4055,12.1875)]
fig,ax=plt.subplots(2,2,figsize=(12,11))
X=np.linspace(-1.5,1.5,500); XX,ZZ2=np.meshgrid(X,X); rr=np.abs(XX)+1e-6; rh=np.hypot(rr,ZZ2); sq=(rr/rh)**2; sg=np.sign(XX); sg[sg==0]=1
for row,(eta,m1,m2) in enumerate(sols):
    a1,a2,r=coeffs(m1,m2,eta)
    g=lambda m,x: x*jn(1,m*x)
    G1=g(m1,1.0)+r*g(m2,1.0); D=m1*g(m1,1.0)+r*m2*g(m2,1.0)
    chiW=np.where(rh<1,(g(m1,rh)+r*g(m2,rh))*sq,G1*sq/rh); IW=np.where(rh<1,(m1*g(m1,rh)+r*m2*g(m2,rh))*sq,D*sq/rh)
    chiU=np.where(rh<1,(a1*g(m1,rh)+a2*r*g(m2,rh))*sq,0.0); IU=np.where(rh<1,(a1*m1*g(m1,rh)+a2*r*m2*g(m2,rh))*sq,0.0)
    for col,(chi,I,name) in enumerate(((chiW,IW,"W (correlation flux: field, moment)"),(chiU,IU,"U (mass drift: spin; confined, zero impulse)"))):
        a_=ax[row,col]; Fp=I/rr*sg; m=np.nanpercentile(abs(Fp),99)
        a_.pcolormesh(XX,ZZ2,Fp,cmap="RdBu_r",vmin=-m,vmax=m,shading="auto")
        a_.contour(XX,ZZ2,chi,levels=24,colors="k",linewidths=0.6)
        t=np.linspace(0,2*np.pi,200); a_.plot(np.cos(t),np.sin(t),"g--",lw=1); a_.set_aspect("equal")
        a_.set_title(f"η = {eta}: μ₁a = {m1:.2f}, μ₂a = {m2:.2f}\n{name}",fontsize=10); a_.set_xlabel("x / a"); a_.set_ylabel("z / a")
fig.suptitle("Isolated Hall-like (double-Beltrami) rings: q = 1, zero impulse, smooth matching to an exterior dipole\n(lines: poloidal; colour: toroidal component)",fontsize=12)
plt.tight_layout(); plt.savefig("double_beltrami.png",dpi=105)
