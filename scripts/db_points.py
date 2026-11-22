import numpy as np
from scipy.optimize import brentq
from scipy.special import spherical_jn as jn
exec(open("double_beltrami.py").read().split("results={}")[0])
res=np.load("db_results.npy",allow_pickle=True).item()
for eta in (1.5,3.0):
    c=res[eta]; c=c[np.isfinite(c[:,2])]; found=[]
    for i in range(len(c)):
        for j in range(i+1,len(c)):
            if abs(c[i,0]-c[j,0])<0.11 and abs(c[i,1]-c[j,1])<0.6 and (c[i,2]-1)*(c[j,2]-1)<0:
                # refine along the zero-impulse branch: parametrize by mu1, re-solve mu2 near the pair
                def qres(m1):
                    lo,hi=min(c[i,1],c[j,1])-0.4,max(c[i,1],c[j,1])+0.4
                    xs=np.linspace(lo,hi,80); v=[impulse(m1,x,eta) for x in xs]
                    for k in range(79):
                        if np.isfinite(v[k]) and np.isfinite(v[k+1]) and v[k]*v[k+1]<0:
                            m2=brentq(lambda x: impulse(m1,x,eta),xs[k],xs[k+1]); return q_of(m1,m2,eta)[0]-1, m2
                    return np.nan,None
                try:
                    a,b=sorted((c[i,0],c[j,0]))
                    if a==b: continue
                    fa,_=qres(a); fb,_=qres(b)
                    if np.isfinite(fa) and np.isfinite(fb) and fa*fb<0:
                        m1=brentq(lambda x: qres(x)[0],a,b,xtol=1e-3); _,m2=qres(m1)
                        if not any(abs(m1-f[0])<0.05 and abs(m2-f[1])<0.05 for f in found):
                            a1,a2,r=coeffs(m1,m2,eta); q,rO=q_of(m1,m2,eta)
                            G1=jn(1,m1)+r*jn(1,m2)                      # W exterior dipole (chi coefficient)
                            found.append((m1,m2,r,a1,a2,q,rO,G1))
                except Exception: pass
    print(f"\neta = {eta}: isolated solutions (zero impulse, smooth W matching, q = 1):")
    print("   mu1 a    mu2 a    c2/c1     alpha1   alpha2    q       r_O/a   W exterior dipole")
    for f in sorted(found): print("  %7.4f  %7.4f  %+8.4f  %+7.3f  %+7.3f  %.4f  %.3f   %+.4f"%f)
