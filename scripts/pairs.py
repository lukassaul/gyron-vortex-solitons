import numpy as np
from scipy.special import spherical_jn as jn
from scipy.integrate import quad
curve=np.load("q1curve.npy")
def shape(lam,b):
    G=lambda r: r*jn(1,abs(lam)*r)+b/lam**2*r**2; e=1e-6; dG=lambda r:(G(r+e)-G(r-e))/(2*e)
    return G,dG
def bnd(lam,b,th):     # on rho=1-: chi, |f|^2 (A=1)
    G,dG=shape(lam,b); rho=1-1e-9; s=np.sin(th); c=np.cos(th); r=rho*s
    chi=G(rho)*s*s; cr=dG(rho)*s*s; ct=G(rho)*2*s*c
    f2=(cr**2+(ct/rho)**2+(lam*chi)**2)/r**2
    return chi,f2
RHO,TH=np.meshgrid(np.linspace(0.02,0.999,160),np.linspace(0.02,np.pi/2,80),indexing="ij")
_pk={}
def peak(lam,b):       # peak |f| inside (A=1), vectorized, cached
    key=(round(abs(lam),4),round(b,4))
    if key in _pk: return _pk[key]
    G,dG=shape(lam,b); s=np.sin(TH); c=np.cos(TH); r=RHO*s
    f=np.sqrt((dG(RHO)*s*s)**2+(G(RHO)*2*s*c/RHO)**2+(lam*G(RHO)*s*s)**2)/r
    _pk[key]=f.max(); return _pk[key]
rows=[]
for (l1,b1,G1,r1) in curve:
    for (l2,b2,G2,r2) in curve:
        for hand in (+1,-1):                # population - same (+1) or opposite (-1) handedness
            L2=hand*l2
            A1=1.0; A2=-A1*G1/G2           # zero weighted impulse: A1 G1 + A2 G2 = 0
            def k(th):
                c1,f1=bnd(l1,b1,th); c2,f2=bnd(L2,b2,th)
                return (A1**2*(b1*c1-0.5*f1)-A2**2*(b2*c2-0.5*f2))*0.5*np.sin(th)
            th=np.linspace(1e-4,np.pi-1e-4,801); Q2=4*np.pi*np.trapezoid(k(th),th)
            D1=l1*A1*G1-L2*A2*G2
            if abs(D1)<1e-6: continue
            pk=max(peak(l1,b1)*abs(A1),peak(L2,b2)*abs(A2))
            Amax=1/pk; a_max=Amax*abs(Q2)/(8*np.pi*abs(D1))
            rows.append((l1,L2,Q2,D1,pk,a_max))
rows.sort(key=lambda x:-x[5])
print(" lam+    lam-     Q/A^2     D/A    peak|f|/A   a_max/xi (spin mu=mu_B, peak V=c)")
for r in rows[:12]: print(" %.2f  %+.2f  %+8.3f  %+7.3f  %7.2f  %8.4f"%r)
prop=[r for r in rows if abs(r[0]-3.6)<0.01 and abs(r[1]+3.6)<0.01]
print("\nreference (same shape, opposite handedness):",prop[:1])
