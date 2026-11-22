import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import spherical_jn as jn
from scipy.optimize import minimize_scalar, brentq
from scipy.interpolate import RegularGridInterpolator as RGI
# shape: chi = [rho j1(lam rho) + (b/lam^2) rho^2] sin^2 theta, I = lam chi ; b free (surface current allowed)
def Gf(lam,b): return lambda r: r*jn(1,lam*r)+b/lam**2*r**2
n=360; R=np.linspace(1e-3,1,n); Z=np.linspace(-1,1,2*n); RR,ZZ=np.meshgrid(R,Z)
rho=np.hypot(RR,ZZ); s2=(RR/np.maximum(rho,1e-9))**2
def qmid(lam,b):
    G=Gf(lam,b); chi=np.where(rho<1,G(rho)*s2,np.nan)
    res=minimize_scalar(lambda r:-G(r),bounds=(0.05,0.999),method="bounded"); chiO=G(res.x); chiB=G(1.0)
    if chiO<=chiB*1.002 or chiO<=0: return np.nan
    gz,gr=np.gradient(chi,Z,R); lev=chiB+0.5*(chiO-chiB)
    cs=plt.contour(RR,ZZ,chi,levels=[lev]); plt.close()
    ps=[p for p in cs.allsegs[0] if np.allclose(p[0],p[-1],atol=5e-3) and len(p)>20]
    if not ps: return np.nan
    p=max(ps,key=len); r=p[:,0]; z=p[:,1]; dl=np.hypot(np.diff(r),np.diff(z)); rm=(r[1:]+r[:-1])/2; zm=(z[1:]+z[:-1])/2
    f=lambda F: RGI((Z,R),np.nan_to_num(F))(np.c_[zm,rm])
    return np.sum(lam*lev/(rm*np.hypot(f(gr),f(gz)))*dl)/(2*np.pi)
curve=[]
for lam in np.linspace(2.4,5.6,17):
    bs=np.linspace(-12,12,49); qs=np.array([qmid(lam,b) for b in bs])
    for i in range(len(bs)-1):
        if np.isfinite(qs[i]) and np.isfinite(qs[i+1]) and (qs[i]-1)*(qs[i+1]-1)<0:
            b=brentq(lambda x: qmid(lam,x)-1,bs[i],bs[i+1],xtol=1e-3); G=Gf(lam,b)
            rO=minimize_scalar(lambda r:-G(r),bounds=(0.05,0.999),method="bounded").x
            curve.append((lam,b,G(1.0),rO)); print(f"lam a={lam:.2f}  b={b:+.3f}  G(1)={G(1.0):+.4f}  r_O={rO:.3f}  matched-b={-lam**3*jn(0,lam)/3:+.3f}")
np.save("q1curve.npy",np.array(curve))
