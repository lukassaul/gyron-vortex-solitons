"""Find isolated Hall-like (double-Beltrami) solutions: smooth W matching, zero impulse, q_W = 1, for each eta.
Usage: python batch/scan_hall.py eta1 eta2 ...   -> batch/solutions_eta<eta>.json"""
import sys, json, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import spherical_jn as jn
from scipy.optimize import brentq
from scipy.interpolate import RegularGridInterpolator as RGI
def coeffs(m1,m2,eta):
    s=2*eta*(m1+m2)/(1+eta); disc=s*s-4*eta*m1*m2
    if disc<=0: return None
    d=np.sqrt(disc); return (2*eta*m1-s)/d,(2*eta*m2-s)/d,-(m1**2*jn(0,m1))/(m2**2*jn(0,m2))
def impulse(m1,m2,eta):
    c=coeffs(m1,m2,eta)
    return np.nan if c is None else c[0]*jn(1,m1)+c[1]*c[2]*jn(1,m2)
n=240; R=np.linspace(1e-3,1,n); Z=np.linspace(-1,1,2*n); RR,ZZ=np.meshgrid(R,Z); rho=np.hypot(RR,ZZ); s2=(RR/np.maximum(rho,1e-9))**2
def q_of(m1,m2,eta):
    c=coeffs(m1,m2,eta)
    if c is None: return np.nan
    a1,a2,r=c; G=lambda x: x*jn(1,m1*x)+r*x*jn(1,m2*x); H=lambda x: m1*x*jn(1,m1*x)+r*m2*x*jn(1,m2*x)
    rs=np.linspace(0.02,0.999,600); gv=G(rs); sg=np.sign(gv[np.argmax(abs(gv))]) or 1; gv=sg*gv
    iO=np.argmax(gv); chiO=gv[iO]; chiB=max(sg*G(1.0),0)
    if iO in (0,len(rs)-1) or chiO<=chiB*1.001: return np.nan
    chi=np.where(rho<1,sg*G(rho)*s2,np.nan); I=np.where(rho<1,sg*H(rho)*s2,np.nan)
    lev=chiB+0.5*(chiO-chiB); cs=plt.figure().add_subplot().contour(RR,ZZ,chi,levels=[lev]); plt.close("all")
    ps=[p for p in cs.allsegs[0] if len(p)>20 and np.allclose(p[0],p[-1],atol=5e-3)]
    if not ps: return np.nan
    p=max(ps,key=len); rr,zz=p[:,0],p[:,1]; dl=np.hypot(np.diff(rr),np.diff(zz)); rm,zm=(rr[1:]+rr[:-1])/2,(zz[1:]+zz[:-1])/2
    gz,gr=np.gradient(chi,Z,R); f=lambda F: RGI((Z,R),np.nan_to_num(F))(np.c_[zm,rm])
    return np.sum(f(I)/(rm*np.hypot(f(gr),f(gz)))*dl)/(2*np.pi)
def branch_m2(m1,eta,m2s):
    v=np.array([impulse(m1,x,eta) for x in m2s]); out=[]
    for i in range(len(m2s)-1):
        if np.isfinite(v[i]) and np.isfinite(v[i+1]) and v[i]*v[i+1]<0 and abs(jn(0,m2s[i]))>1e-3:
            out.append(brentq(lambda x: impulse(m1,x,eta),m2s[i],m2s[i+1]))
    return out
for eta in map(float,sys.argv[1:]):
    m1s=np.linspace(0.6,7,65); m2s=np.linspace(0.7,16,500); pts=[]
    for m1 in m1s:
        for m2 in branch_m2(m1,eta,m2s[m2s>m1+0.05]): pts.append((m1,m2,q_of(m1,m2,eta)))
    pts=np.array(pts) if pts else np.zeros((0,3)); sols=[]
    def follow(m1,m2guess):
        cands=branch_m2(m1,eta,np.linspace(max(m1+0.05,m2guess-0.5),m2guess+0.5,60))
        if not cands: return np.nan,None
        m2=min(cands,key=lambda x: abs(x-m2guess)); return q_of(m1,m2,eta)-1,m2
    for i in range(len(pts)):
        for j in range(len(pts)):
            if pts[j,0]-pts[i,0]>0 and pts[j,0]-pts[i,0]<0.11 and abs(pts[i,1]-pts[j,1])<0.5 and np.isfinite(pts[i,2]) and np.isfinite(pts[j,2]) and (pts[i,2]-1)*(pts[j,2]-1)<0:
                try:
                    g0=pts[i,1]; fa,_=follow(pts[i,0],g0); fb,_=follow(pts[j,0],pts[j,1])
                    if not(np.isfinite(fa) and np.isfinite(fb) and fa*fb<0): continue
                    m1=brentq(lambda x: follow(x,g0+(x-pts[i,0])*(pts[j,1]-g0)/(pts[j,0]-pts[i,0]))[0],pts[i,0],pts[j,0],xtol=2e-4)
                    _,m2=follow(m1,g0+(m1-pts[i,0])*(pts[j,1]-g0)/(pts[j,0]-pts[i,0]))
                    if m2 and not any(abs(m1-s[0])<0.03 and abs(m2-s[1])<0.03 for s in sols): sols.append((m1,m2))
                except Exception: pass
    json.dump(dict(eta=eta,branch_points=len(pts),solutions=sols),open(f"batch/solutions_eta{eta:g}.json","w"))
    print(f"eta={eta}: {len(pts)} zero-impulse points, {len(sols)} isolated q=1 solutions: {[(round(a,3),round(b,3)) for a,b in sols]}",flush=True)
