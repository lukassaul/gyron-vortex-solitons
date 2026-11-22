"""Spin, charge, moment and g for the isolated Hall-like rings; size from energy at fixed spin J = hbar/2.
Shape units a=1, amplitude A (velocities in c). Charge: K = -2 w n U.W/eta (per-class Bernoulli with Phi = g_Phi K, holding outside).
J = nu A a^4 j_s ; E_kin = (1/2) nu A^2 a^3 K_s ; bubble p0 (4pi/3) a^3 + 4 pi sig a^2 ; hbar = kappa (units rho_m c xi^4).
mu/mu_B = 8 pi D_s a/(Q_s A) ; g = 2 mu/mu_B at J = hbar/2."""
import numpy as np
from scipy.special import spherical_jn as jn
from scipy.optimize import minimize_scalar
exec(open("double_beltrami.py").read().split("def impulse")[0])
p0=5/7; sig=2*np.sqrt(2)/3
rho=np.linspace(1e-3,0.9995,500); th=np.linspace(1e-3,np.pi-1e-3,400); RH,TH=np.meshgrid(rho,th,indexing="ij")
s=np.sin(TH); c=np.cos(TH); r=RH*s; w=2*np.pi*RH**2*s; dV=(rho[1]-rho[0])*(th[1]-th[0])
def props(eta,m1,m2):
    a1,a2,rr=coeffs(m1,m2,eta)
    g=lambda m,x: x*jn(1,m*x); gp=lambda m,x: jn(1,m*x)+m*x*jn(1,m*x,derivative=True)
    GW=g(m1,RH)+rr*g(m2,RH); GWp=gp(m1,RH)+rr*gp(m2,RH); HW=m1*g(m1,RH)+rr*m2*g(m2,RH)
    GU=a1*g(m1,RH)+a2*rr*g(m2,RH); GUp=a1*gp(m1,RH)+a2*rr*gp(m2,RH); HU=a1*m1*g(m1,RH)+a2*rr*m2*g(m2,RH)
    def comps(G,Gp,H):   # poloidal (r, theta) and toroidal components times r
        return Gp*s*s, G*2*s*c/RH, H*s*s
    Wr,Wt,Wp=comps(GW,GWp,HW); Ur,Ut,Up=comps(GU,GUp,HU)
    UW=(Ur*Wr+Ut*Wt+Up*Wp)/r**2; U2=(Ur**2+Ut**2+Up**2)/r**2; W2=(Wr**2+Wt**2+Wp**2)/r**2
    Qs=4*np.pi*np.trapezoid(-UW[-1,:]/eta*0.5*s[-1,:],th)
    js=np.sum(HU*s*s*w)*dV
    G1=g(m1,1.0)+rr*g(m2,1.0); D=m1*g(m1,1.0)+rr*m2*g(m2,1.0)
    ro=np.linspace(1,200,4000); RO,TO=np.meshgrid(ro,th,indexing="ij"); so=np.sin(TO); co=np.cos(TO)
    Wext2=((G1*2*co/RO**3)**2+(G1*so/RO**3)**2+(D*so/RO**2)**2)
    Ks=np.sum((U2+W2)*w)*dV+np.sum(Wext2*2*np.pi*RO**2*so)*(ro[1]-ro[0])*(th[1]-th[0])
    return dict(Qs=Qs,js=js,D=D,Ks=Ks)
def size(P,nu):
    Qs,js,D,Ks=P["Qs"],P["js"],P["D"],P["Ks"]; kap=50.0
    for _ in range(200):
        E=lambda a: kap**2*Ks/(8*nu*a**5*js**2)+p0*4*np.pi/3*a**3+4*np.pi*sig*a**2
        res=minimize_scalar(lambda la: E(np.exp(la)),bounds=(np.log(1e-3),np.log(1e4)),method="bounded"); a=np.exp(res.x); En=E(a)
        if abs(En-kap)<1e-9*kap: break
        kap=0.5*kap+0.5*En
    A=kap/(2*nu*a**4*abs(js)); g=16*np.pi*D*a/(Qs*A)
    return a,A,kap,g
sols=[(1.5,3.8503,8.8136),(1.5,4.0723,8.2376),(1.5,4.7219,13.4535),(3.0,3.4055,12.1875)]
print(" eta   mu1a    mu2a  |  Q_s      j_s       D_s      K_s   | sign(D j Q) | g at nu = 0.1, 1, 10, 100   (a/xi, A at nu=10)")
for eta,m1,m2 in sols:
    P=props(eta,m1,m2); gs=[]; 
    for nu in (0.1,1,10,100):
        a,A,kap,g=size(P,nu); gs.append(g)
        if nu==10: a10,A10=a,A
    print(f" {eta:3.1f} {m1:6.3f} {m2:7.3f} | {P['Qs']:+.4f} {P['js']:+.5f} {P['D']:+.4f} {P['Ks']:7.3f} |    {'+' if P['D']*P['js']*P['Qs']>0 else '-'}       | "+"  ".join(f"{x:+8.3f}" for x in gs)+f"   ({a10:.3f}, {A10:.3f})")
np.save("db_props.npy",np.array([[*s_,*props(*s_).values()] for s_ in sols]))
