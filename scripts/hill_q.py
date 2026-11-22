import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import spherical_jn as jn
from scipy.optimize import brentq, minimize_scalar
# interior chi = [A g(rho) + (b/lam^2) rho^2] sin^2, g = rho j1(lam rho); I = lam chi
# matching I, dI continuous at a=1  ->  beta = -A lam^3 j0(lam)/3
A=1.0
def setup(lam):
    beta=-A*lam**3*jn(0,lam)/3
    G=lambda rho: A*rho*jn(1,lam*rho)+beta/lam**2*rho**2
    return beta,G
def q_profile(lam,n=500):
    beta,G=setup(lam)
    R=np.linspace(1e-3,1,n); Z=np.linspace(-1,1,2*n); RR,ZZ=np.meshgrid(R,Z)
    rho=np.hypot(RR,ZZ); s2=(RR/np.maximum(rho,1e-9))**2
    chi=np.where(rho<1,G(rho)*s2,np.nan)
    gz,gr=np.gradient(chi,Z,R)
    # O-point on midplane
    res=minimize_scalar(lambda r:-G(r),bounds=(0.05,0.999),method="bounded")
    rO,chiO=res.x,G(res.x); chiB=G(1.0)  # boundary max at theta=pi/2
    if chiO<=chiB: return rO,chiO,chiB,[]
    out=[]
    for frac in (0.15,0.5,0.85):
        lev=chiB+frac*(chiO-chiB)
        cs=plt.contour(RR,ZZ,chi,levels=[lev]); plt.close()
        paths=[p for p in cs.allsegs[0] if np.allclose(p[0],p[-1],atol=5e-3) and len(p)>20]
        if not paths: continue
        p=max(paths,key=len); r=p[:,0]; z=p[:,1]
        dl=np.hypot(np.diff(r),np.diff(z)); rm=0.5*(r[1:]+r[:-1]); zm=0.5*(z[1:]+z[:-1])
        from scipy.interpolate import RegularGridInterpolator as RGI
        f=lambda F: RGI((Z,R),np.nan_to_num(F))(np.c_[zm,rm])
        gch=np.hypot(f(gr),f(gz)); I=lam*lev
        out.append((frac, np.sum(I/(rm*gch)*dl)/(2*np.pi)))
    return rO,chiO,chiB,out
for lam in (2.6,2.9,np.pi,3.4,3.7,4.0):
    rO,chiO,chiB,qs=q_profile(lam)
    print(f"lam={lam:.3f} beta={setup(lam)[0]:+.3f} rO={rO:.3f} closed band={(chiO-chiB)/chiO:.3f}  q={[(f,round(q,4)) for f,q in qs]}")
