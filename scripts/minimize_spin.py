import numpy as np
from scipy.optimize import brentq
from scipy.interpolate import interp1d
rows=np.load("energy_spin_rows.npy"); A_,ew_,ef_,j_,Q_,band_,M_=rows.T
p0=5/7; sig=2*np.sqrt(2)/3
f=lambda y: interp1d(A_,y,kind="cubic")
EW,EF,J,QQ=f(ew_),f(ef_),f(j_),f(Q_)
Aband=A_[np.argmax(band_<0)] if (band_<0).any() else A_[-1]
def at(A,nu,Gam,restore=True):
    E=lambda a: nu*a**3*EW(A)+Gam*a*EF(A)+(p0*4*np.pi/3*a**3+4*np.pi*sig*a**2 if restore else 0)
    g=lambda a: nu*a**4*J(A)-E(a)/2
    a=brentq(g,1e-3,1e6); return a,E(a)
As=np.linspace(A_[0],A_[-1],300)
print(f"closed band lost at A ≈ {Aband:.3f}; choke at A ≈ {A_[-1]:.3f}\n")
print("   nu     Gamma  |  A* (min E)     a/xi     kappa=E    Q~=|Q| a   where")
for nu in (0.1,1,10,100):
    for Gam in (0.01,1,100):
        vals=[at(A,nu,Gam) for A in As]; Es=np.array([v[1] for v in vals]); i=np.argmin(Es)
        a,E=vals[i]; where=("interior minimum" if 0<i<len(As)-1 else ("at the choke" if i==len(As)-1 else "at A -> 0"))
        print(f"  {nu:6g} {Gam:6g}  |  {As[i]:.3f}       {a:8.3f}  {E:10.4g}  {abs(QQ(As[i]))*a:8.4f}   {where}{' (band lost)' if As[i]>Aband else ''}")

print("\npredicted alpha = Gamma Q~^2/(4 pi kappa) at the minimum, and Gamma giving alpha = 1/137:")
from scipy.optimize import brentq as bq
def minimum(nu,Gam):
    vals=[at(A,nu,Gam) for A in As]; Es=np.array([v[1] for v in vals]); i=np.argmin(Es); a,E=vals[i]
    return As[i],a,E,abs(QQ(As[i]))*a, 0<i<len(As)-1
for nu in (3,10,30,100,300):
    out=[]
    for Gam in (1,10,100,1e3,1e4):
        A,a,E,Q,interior=minimum(nu,Gam); out.append(f"Γ={Gam:g}: α={Gam*Q**2/(4*np.pi*E):.2e}{'' if interior else '*'}")
    print(f"  nu={nu:4g}  "+"  ".join(out))
print("  (* = minimum at the choke, not interior)")
