import numpy as np
exec(open("pairs.py").read().split("rows=[]")[0])      # shape(), bnd(), peak()
curve=np.load("q1curve.npy")
def member(l): i=np.argmin(abs(curve[:,0]-l)); return curve[i]
def props(lp,lm):
    l1,b1,G1,_=member(lp); l2,b2,G2,_=member(abs(lm)); L2=np.sign(lm)*l2
    A1=1.0; A2=-A1*G1/G2
    th=np.linspace(1e-4,np.pi-1e-4,801)
    c1,f1=bnd(l1,b1,th); c2,f2=bnd(L2,b2,th)
    Q2=4*np.pi*np.trapezoid((A1**2*(b1*c1-0.5*f1)-A2**2*(b2*c2-0.5*f2))*0.5*np.sin(th),th)
    D1=l1*A1*G1-L2*A2*G2
    pk=max(peak(l1,b1)*abs(A1),peak(L2,b2)*abs(A2))
    A=1/pk                       # carrying flow at the luminal limit
    a=A*abs(Q2)/(8*np.pi*abs(D1))  # size from spin + mu = mu_B
    Qsat=abs(Q2)*A**2*a          # charge in units of e_medium = eps* g_Phi K0 xi
    return a,Qsat,A
print(" pair (λ+, λ-)    g      a/ξ    A_lum   Q~ at luminal limit")
for lp,lm,g in ((5.4,-3.6,2.009),(5.6,-3.6,1.926),(5.2,-3.8,1.931),(5.0,-4.0,1.859),(4.8,-4.0,2.125),(4.6,-4.2,2.065),(4.2,-4.6,2.065),(3.6,-5.4,2.009)):
    a,Q,A=props(lp,lm); print(f" ({lp:.1f}, {lm:+.1f})   {g:.3f}  {a:6.3f}  {A:6.3f}   {Q:.4f}")
