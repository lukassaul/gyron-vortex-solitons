import numpy as np
from scipy.special import spherical_jn as jn, eval_legendre as P
from scipy.optimize import brentq, minimize_scalar
from hill_q import setup
lam=3.6331; beta,G=setup(lam); e=1e-6
dG=lambda r:(G(r+e)-G(r-e))/(2*e); d2G=lambda r:(G(r+e)-2*G(r)+G(r-e))/e**2
# grid (rho,theta) interior
nr,nt=500,400
rho=np.linspace(1e-3,1,nr); th=np.linspace(1e-3,np.pi-1e-3,nt)
RH,TH=np.meshgrid(rho,th,indexing="ij"); s=np.sin(TH); c=np.cos(TH); r=RH*s
chi=G(RH)*s*s; chi_r=dG(RH)*s*s; chi_t=G(RH)*2*s*c
I=lam*chi; I_r=lam*chi_r; I_t=lam*chi_t
Dchi=beta*r**2-lam**2*chi
w=2*np.pi*RH**2*s
dV=(rho[1]-rho[0])*(th[1]-th[0])
F2=(chi_r**2+(chi_t/RH)**2+I**2)/r**2
curl2=(I_r**2+(I_t/RH)**2+Dchi**2)/r**2
Em_in=0.5*np.sum(curl2*w)*dV
# exterior magnetic: chi=G1 s^2/rho, I=lam G1 s^2/rho
G1=G(1.0); ro=np.linspace(1,400,40000); RO,TO=np.meshgrid(ro,th,indexing="ij"); so=np.sin(TO); co=np.cos(TO)
Iext_r=-lam*G1*so**2/RO**2; Iext_t=lam*G1*2*so*co/RO
curl2o=(Iext_r**2+(Iext_t/RO)**2)/(RO*so)**2
Em_out=0.5*np.sum(curl2o*2*np.pi*RO**2*so)*(ro[1]-ro[0])*(th[1]-th[0])
m0=Em_in+Em_out
# electric: k = beta*chi - |F|^2/2 (delta=1, A=1, c=1)
k=beta*chi-0.5*F2
kr=np.gradient(k,rho,axis=0); kt=np.gradient(k,th,axis=1)
Ee_in=0.5*np.sum((kr**2+(kt/RH)**2)*w)*dV
kb=k[-1,:]; Ee_out=0; mono=None
for l in range(0,30):
    kl=(2*l+1)/2*np.trapezoid(kb*P(l,np.cos(th))*np.sin(th),th)
    if l==0: mono=kl
    Ee_out+=0.5*kl**2*(l+1)*4*np.pi/(2*l+1)
q0=abs(4*np.pi*mono); e0=Ee_in+Ee_out; cE=e0/q0**2
hel=-2*np.sum(I*Dchi/r**2*w)*dV  # interior helicity (exterior small)
# closed band half-width at midplane
rO=minimize_scalar(lambda x:-G(x),bounds=(0.3,0.999),method="bounded").x
rin=brentq(lambda x:G(x)-G1,0.3,rO); bw=(1-rin)/2
print(f"m0={m0:.3f} (in {Em_in:.3f}, out {Em_out:.3f})  q0={q0:.3f}  e0={e0:.4f}  c_E={cE:.4f} (point charge 1/8pi={1/(8*np.pi):.4f})")
print(f"h_in={hel:.3f}  rO/a={rO:.4f}  closed band r in [{rin:.4f},1.0], half-width b/a={bw:.4f}")
np.save("coef.npy",[m0,q0,cE,rO,bw,hel])
