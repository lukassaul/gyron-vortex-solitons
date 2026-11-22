"""Simplest Hall-like (double-Beltrami) spherical ring, a = 1, l = 1, two classes.
W = c1 W(mu1) + c2 W(mu2), W(mu): chi = rho j1(mu rho) sin^2, I = mu chi.  U = a1 c1 W(mu1) + a2 c2 W(mu2).
alpha_i = (2 eta mu_i - s)/d, s = 2 eta (mu1+mu2)/(1+eta), d = sqrt(s^2 - 4 eta mu1 mu2)   (h = 1, Lambda = eta)
Conditions: (1) W matches an exterior dipole smoothly: c1 mu1^2 j0(mu1) + c2 mu2^2 j0(mu2) = 0
            (2) zero impulse: U has no exterior dipole: a1 c1 j1(mu1) + a2 c2 j1(mu2) = 0
            (3) closure q_W = 1 on the closed tori."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import spherical_jn as jn
from scipy.optimize import brentq, minimize_scalar
from scipy.interpolate import RegularGridInterpolator as RGI
def coeffs(m1,m2,eta):
    s=2*eta*(m1+m2)/(1+eta); disc=s*s-4*eta*m1*m2
    if disc<=0: return None
    d=np.sqrt(disc); a1=(2*eta*m1-s)/d; a2=(2*eta*m2-s)/d
    r=-(m1**2*jn(0,m1))/(m2**2*jn(0,m2))          # c2/c1 from (1)
    return a1,a2,r
def impulse(m1,m2,eta):
    c=coeffs(m1,m2,eta)
    if c is None: return np.nan
    a1,a2,r=c; return a1*jn(1,m1)+a2*r*jn(1,m2)
n=300; R=np.linspace(1e-3,1,n); Z=np.linspace(-1,1,2*n); RR,ZZ=np.meshgrid(R,Z); rho=np.hypot(RR,ZZ); s2=(RR/np.maximum(rho,1e-9))**2
def q_of(m1,m2,eta):
    a1,a2,r=coeffs(m1,m2,eta)
    G=lambda x: x*jn(1,m1*x)+r*x*jn(1,m2*x); H=lambda x: m1*x*jn(1,m1*x)+r*m2*x*jn(1,m2*x)
    sgn=np.sign(G(0.5)) or 1.0
    chi=np.where(rho<1,sgn*G(rho)*s2,np.nan); I=np.where(rho<1,sgn*H(rho)*s2,np.nan)
    rs=np.linspace(0.02,0.999,800); gv=sgn*G(rs); iO=np.argmax(gv); chiO=gv[iO]; chiB=sgn*G(1.0)
    if iO in (0,len(rs)-1) or chiO<=max(chiB,0)*1.001: return np.nan, None
    gz,gr=np.gradient(chi,Z,R); lev=max(chiB,0)+0.5*(chiO-max(chiB,0))
    cs=plt.contour(RR,ZZ,chi,levels=[lev]); plt.close()
    ps=[p for p in cs.allsegs[0] if np.allclose(p[0],p[-1],atol=5e-3) and len(p)>20]
    if not ps: return np.nan, None
    p=max(ps,key=len); rr=p[:,0]; zz=p[:,1]; dl=np.hypot(np.diff(rr),np.diff(zz)); rm=(rr[1:]+rr[:-1])/2; zm=(zz[1:]+zz[:-1])/2
    f=lambda F: RGI((Z,R),np.nan_to_num(F))(np.c_[zm,rm])
    q=np.sum(f(I)/(rm*np.hypot(f(gr),f(gz)))*dl)/(2*np.pi)
    return q, rs[iO]
results={}
for eta in (1.5,3.0,10.0):
    curve=[]
    for m1 in np.linspace(0.5,6,56):
        m2s=np.linspace(m1+0.05,14,400); v=np.array([impulse(m1,m2,eta) for m2 in m2s])
        for i in range(len(m2s)-1):
            if np.isfinite(v[i]) and np.isfinite(v[i+1]) and v[i]*v[i+1]<0 and abs(jn(0,m2s[i]))>1e-3:
                m2=brentq(lambda x: impulse(m1,x,eta),m2s[i],m2s[i+1]); q,rO=q_of(m1,m2,eta)
                curve.append((m1,m2,q,rO if rO else np.nan))
    curve=np.array(curve); results[eta]=curve
    ok=np.isfinite(curve[:,2]) if len(curve) else []
    print(f"eta = {eta}: zero-impulse points {len(curve)}, with a closed band {int(np.sum(ok))};  q range on them: {np.nanmin(curve[:,2]) if np.any(ok) else None} .. {np.nanmax(curve[:,2]) if np.any(ok) else None}")
np.save("db_results.npy",results,allow_pickle=True)
