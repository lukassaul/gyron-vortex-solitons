import numpy as np
from scipy.special import spherical_jn as jn
curve=np.load("q1curve.npy")
nr,nt=300,300
rho=np.linspace(1e-3,0.9995,nr); th=np.linspace(1e-3,np.pi-1e-3,nt); RH,TH=np.meshgrid(rho,th,indexing="ij")
s=np.sin(TH); c=np.cos(TH); r=RH*s; w=2*np.pi*RH**2*s; dV=(rho[1]-rho[0])*(th[1]-th[0])
def shp(lam,b):
    G=lambda x: x*jn(1,abs(lam)*x)+b/lam**2*x**2; e=1e-6; dG=lambda x:(G(x+e)-G(x-e))/(2*e)
    chi=G(RH)*s*s; cr=dG(RH)*s*s; ct=G(RH)*2*s*c; I=lam*chi
    f2=(cr**2+(ct/RH)**2+I**2)/r**2
    return dict(chi=chi,I=I,f2=f2,G1=G(1.0))
res=[]
for (l1,b1,_,_) in curve:
    for (l2,b2,_,_) in curve:
        for hand in (+1,-1):
            L2=hand*l2; P=shp(l1,b1); M=shp(L2,b2)
            A1=1.0; A2=-A1*P["G1"]/M["G1"]
            # charge (boundary average of net dk), net moment D, interior J and kinetic K
            kb=(A1**2*(b1*P["chi"]-0.5*P["f2"])-A2**2*(b2*M["chi"]-0.5*M["f2"]))[-1,:]
            Q=4*np.pi*np.trapezoid(kb*0.5*np.sin(th),th)
            D=l1*A1*P["G1"]-L2*A2*M["G1"]
            j=np.sum((A1*P["I"]+A2*M["I"])*w)*dV
            K_in=np.sum((A1**2*P["f2"]+A2**2*M["f2"])*w)*dV
            # exterior: net dipole (chi coef 2*A1*G1, I coef D); each population carries half the net flux
            Cc=A1*P["G1"]-A2*M["G1"]
            K_out=0.5*(4*np.pi/3)*(4*Cc**2*0+0)  # placeholder replaced below
            # exterior kinetic: populations each |F_net|/2 -> K_out = 1/2 * int |f_net,ext|^2
            ro=np.linspace(1,300,6000); RO,TO=np.meshgrid(ro,th,indexing="ij"); so=np.sin(TO); co=np.cos(TO)
            fr=Cc*2*co/RO**3; ft=Cc*so/RO**3; fp=D*so/RO**2
            K_out=0.5*np.sum((fr**2+ft**2+fp**2)*2*np.pi*RO**2*so)*(ro[1]-ro[0])*(th[1]-th[0])
            K=K_in+K_out
            if abs(Q)<1e-6 or abs(j)<1e-9: continue
            g=4*np.pi*D*K/(Q*j)
            res.append((l1,L2,Q,D,j,K,g))
res=np.array(res)
print(" lam+    lam-       Q        D        J~       K~       g")
for row in sorted(res,key=lambda x:abs(abs(x[6])-2))[:12]: print(" %.2f  %+.2f  %+8.3f %+8.4f %+8.4f %8.3f  %+8.3f"%tuple(row))
print("\nbest-size pairs (from moment condition):")
for row in res:
    if (abs(row[0]-5.4)<0.01 and abs(row[1]+5.2)<0.01) or (abs(row[0]-5.2)<0.01 and abs(row[1]+5.0)<0.01):
        print(" %.2f  %+.2f  %+8.3f %+8.4f %+8.4f %8.3f  %+8.3f"%tuple(row))
print("\nrange of g over all pairs: same-hand", np.nanmin(abs(res[res[:,1]>0,6])),np.nanmax(abs(res[res[:,1]>0,6])),
      " opposite-hand",np.nanmin(abs(res[res[:,1]<0,6])),np.nanmax(abs(res[res[:,1]<0,6])))
