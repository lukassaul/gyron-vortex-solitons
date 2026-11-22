"""q=1 quantized ring (lambda a=3.633, delta=+1), Gamma=1e5 picture B: a=149.8 xi, A=0.078.
Lengths in units of a. Velocities in units of c."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe, eval_legendre as P
from hill_q import setup
lam=3.6331; beta,G=setup(lam); G1=G(1.0); atil=18.86; Mach=0.222; rO=0.8469
xi=1/atil
n=420; R=np.linspace(0.02,1.6,n); Z=np.linspace(-1.6,1.6,2*n); RR,ZZ=np.meshgrid(R,Z)
rho=np.hypot(RR,ZZ); s2=(RR/rho)**2; inside=rho<1
chi=np.where(inside,G(rho)*s2,G1*s2/rho); I=np.where(inside,lam*G(rho)*s2,lam*G1*s2/rho)
dz,dr=Z[1]-Z[0],R[1]-R[0]
chi_z,chi_r=np.gradient(chi,dz,dr); I_z,I_r=np.gradient(I,dz,dr)
# weighted flux F (shape units) -> drift of carrying population V = Mach*F/k+ (k+~1)
Fr,Fz,Fp=-chi_z/RR, chi_r/RR, I/RR
Vr,Vz,Vp=Mach*Fr,Mach*Fz,Mach*Fp
# spin density ~ B ~ curl F : poloidal B = grad I x grad phi, toroidal -Delta*chi/r
Br,Bz=-I_z/RR, I_r/RR
Dchi=np.where(inside,beta*RR**2-lam**2*chi,0.0); Bp=-Dchi/RR
# carrying-population density and net spin pressure K = k+ - k- (inside), Laplace outside
F2=Fr**2+Fz**2+Fp**2
kin=Mach**2*(beta*chi-0.5*F2)
th=np.linspace(1e-3,np.pi-1e-3,400); rb=0.999
def kb_of(t):
    r_=rb*np.sin(t); z_=rb*np.cos(t); i=np.clip(((z_-Z[0])/dz).astype(int),0,2*n-1); j=np.clip(((r_-R[0])/dr).astype(int),0,n-1)
    return kin[i,j]
kb=kb_of(th); kl=[(2*l+1)/2*np.trapezoid(kb*P(l,np.cos(th))*np.sin(th),th) for l in range(20)]
ct=ZZ/rho; kout=sum(kl[l]*rho**(-(l+1))*P(l,ct) for l in range(20))
K=np.where(inside,kin,kout); kplus=1+np.where(inside,kin,0)
# bulk vortex ring at r_O: Stokes stream function, circulation 2*pi*zeta = 2*pi*xi (units a, c)
Gb=2*np.pi*xi
m=np.clip(4*RR*rO/((RR+rO)**2+ZZ**2),0,1-1e-12); k_=np.sqrt(m)
psi=Gb/(2*np.pi)*np.sqrt(RR*rO)*((2/k_-k_)*ellipk(m)-2/k_*ellipe(m))
psi_z,psi_r=np.gradient(psi,dz,dr); ur,uz=-psi_z/RR, psi_r/RR
d=np.hypot(RR-rO,ZZ)/xi; nb=d**2/(d**2+2)        # GP core density (Pade approx)
ur,uz=np.where(d<0.6,np.nan,ur),np.where(d<0.6,np.nan,uz)
circ=lambda ax: (ax.plot(np.cos(np.linspace(-np.pi/2,np.pi/2,200)),np.sin(np.linspace(-np.pi/2,np.pi/2,200)),"w--",lw=1),
                 ax.contour(RR,ZZ,np.where(inside,chi,np.nan),levels=[G1],colors="m",linewidths=1.2))
fig,ax=plt.subplots(2,3,figsize=(17,12))
a0=ax[0,0]; sp=np.hypot(np.hypot(Vr,Vz),Vp)
pc=a0.pcolormesh(RR,ZZ,sp,cmap="viridis",shading="auto"); fig.colorbar(pc,ax=a0,label="|V| / c")
a0.streamplot(R,Z,Vr,Vz,color="w",density=1.6,linewidth=0.7,arrowsize=0.8); circ(a0)
a0.set_title("(a) weighted drift V of the carrying population\nstreamlines: poloidal; colour: |V|/c (incl. toroidal)")
a1=ax[0,1]; pc=a1.pcolormesh(RR,ZZ,nb,cmap="Greys_r",shading="auto",vmin=0,vmax=1); fig.colorbar(pc,ax=a1,label="bulk density n/n₀")
a1.streamplot(R,Z,ur,uz,color=np.hypot(ur,uz),cmap="autumn",density=1.8,linewidth=0.8,arrowsize=0.8); circ(a1)
a1.set_title("(b) bulk flow (vortex ring, circulation 2πζ at the O-ring)\ngrey: bulk density (GP core, width ξ)")
a2=ax[0,2]; m_=np.nanmax(abs(Bp)); pc=a2.pcolormesh(RR,ZZ,Bp,cmap="RdBu_r",vmin=-m_,vmax=m_,shading="auto"); fig.colorbar(pc,ax=a2,label="toroidal spin density (∝ B_φ)")
a2.streamplot(R,Z,Br,Bz,color="k",density=1.6,linewidth=0.7,arrowsize=0.8); circ(a2)
a2.set_title("(c) spin streamlines: S = I_s B ∝ ∇×F\nstreamlines: poloidal; colour: toroidal component")
a3=ax[1,0]; m_=np.nanmax(abs(K)); pc=a3.pcolormesh(RR,ZZ,K,cmap="PuOr_r",vmin=-m_,vmax=m_,shading="auto"); fig.colorbar(pc,ax=a3,label="net spin pressure K/K₀")
a3.contour(RR,ZZ,K,levels=15,colors="k",linewidths=0.4); circ(a3)
a3.set_title("(d) net weighted spin pressure K = K₊ − K₋\n(Laplace outside: Coulomb + multipoles)")
for a_ in (a0,a1,a2,a3): a_.set_aspect("equal"); a_.set_xlabel("r / a"); a_.set_ylabel("z / a"); a_.plot(rO,0,"o",color="orange",mec="k",ms=6)
# 1D profiles
iz0=np.argmin(abs(Z)); ir0=0
a4=ax[1,1]
a4.plot(R,nb[iz0],"k",label="bulk density n/n₀")
a4.plot(R,kplus[iz0],"C0",label="carrying population K₊/K₀")
a4.plot(R,K[iz0],"C1",label="net spin pressure K/K₀ ")
a4.plot(R,sp[iz0],"C2",label="|V|/c (weighted drift)")
a4.plot(R,np.clip(np.hypot(ur,uz)[iz0],0,1),"C3",ls="--",label="|u|/c (bulk)")
a4.axvline(rO,color="orange",ls=":"); a4.axvline(1,color="g",ls=":"); a4.axvspan(0.697,1,color="m",alpha=0.08)
a4.set_title("(e) midplane profiles (z = 0)"); a4.set_xlabel("r / a"); a4.legend(fontsize=8); a4.set_ylim(-0.25,1.25)
a5=ax[1,2]; ZZa=Z; ja=0
a5.plot(Z,nb[:,ja],"k",label="bulk density"); a5.plot(Z,kplus[:,ja],"C0",label="K₊/K₀")
a5.plot(Z,K[:,ja],"C1",label="K/K₀ "); a5.plot(Z,sp[:,ja],"C2",label="|V|/c")
a5.plot(Z,np.hypot(ur,uz)[:,ja],"C3",ls="--",label="|u|/c (bulk)")
a5.set_title("(f) profiles along the axis (r → 0)"); a5.set_xlabel("z / a"); a5.legend(fontsize=8); a5.set_ylim(-0.25,1.25)
fig.suptitle("Quantized ring, q = 1 (λa = 3.633), δ = +1, Γ = 10³ (thin core): a = 18.9ξ, r_O = 0.85a, A = 0.222.  Magenta: closed (1,1) tori",fontsize=12)
plt.tight_layout(); plt.savefig("ring_profiles.png",dpi=120)
sh=np.nanmax(sp); print("peak|V|/c",sh,"per A",sh/Mach,"  min K",np.nanmin(K),"per A^2",np.nanmin(K)/Mach**2,"  min K+",np.nanmin(kplus))
