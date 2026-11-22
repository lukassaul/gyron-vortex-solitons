"""Compressible isothermal Bragg-Hawthorne ring (one population, q=1 shape, lam a=3.633).
div(grad chi/(k r^2)) = k (beta - lam^2 chi / r^2),   ln k = beta chi - (|grad chi|^2/(k^2 r^2) + lam^2 chi^2/r^2)/2   (c=1)
Boundary: chi = A G1 sin^2 on rho=1 (exterior dipole), chi=0 on axis. Raise A until no subsonic solution (choke)."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from hill_q import setup
lam=3.6331; b,G=setup(lam); G1=G(1.0)
nr,nz=70,140; R=np.linspace(0,1.05,nr+1)[1:]; Z=np.linspace(-1.05,1.05,nz); dr=R[1]-R[0]; dz=Z[1]-Z[0]
RR,ZZ=np.meshgrid(R,Z,indexing="ij"); RHO=np.hypot(RR,ZZ); inside=RHO<1
idx=-np.ones(RR.shape,int); idx[inside]=np.arange(inside.sum())
def bc(A): return A*G1*(RR/np.maximum(RHO,1e-9))**2/np.maximum(RHO,1e-9)   # exterior values
def solve_k(chi,A):
    cx=np.vstack([np.zeros((1,nz)),chi]); Rx=np.concatenate([[0],R])        # include the axis (chi=0)
    gr=np.gradient(cx,Rx,axis=0)[1:]; gz=np.gradient(chi,dz,axis=1); g2=(gr**2+gz**2)/RR**2; p2=(lam*chi/RR)**2
    C=A*b*chi-0.5*p2       # beta = A*b (beta scales with amplitude)
    k=np.ones_like(chi); ok=np.ones_like(chi,bool)
    for _ in range(40):
        f=np.log(k)+0.5*g2/k**2-C; fp=1/k-g2/k**3
        k=np.clip(k-f/np.where(abs(fp)<1e-9,1e-9,fp),1e-3,5)
    g=np.sqrt(g2); sub=(C>=np.log(np.maximum(g,1e-12))+0.5)|(g<1e-9)   # subsonic root exists
    deep=RHO<0.95                 # boundary-adjacent nodes straddle the vortex sheet: exclude from the choke test
    k=np.where((k<g)&deep,np.nan,np.maximum(k,np.minimum(g,1.0)))
    return k, (sub|~deep)&np.isfinite(k), np.where(deep,np.sqrt(g2)/k,0)
def assemble(k,A):
    N=inside.sum(); rows=[];cols=[];vals=[]; rhs=np.zeros(N); chib=bc(A)
    kk=np.where(np.isfinite(k),k,1.0)
    for i,j in zip(*np.nonzero(inside)):
        n=idx[i,j]; diag=0.0
        for (di,dj,h) in ((1,0,dr),(-1,0,dr),(0,1,dz),(0,-1,dz)):
            ii,jj=i+di,j+dj
            rf=R[i]+di*dr/2; kf=0.5*(kk[i,j]+(kk[ii,jj] if 0<=ii<nr and 0<=jj<nz else 1.0))
            coef=(1.0/(kf*R[i]*rf*h*h) if di!=0 else 1.0/(kf*R[i]**2*h*h)) if rf>0 else 0.0
            if ii<0:   # axis: chi=0
                diag-=coef; continue
            diag-=coef
            if 0<=ii<nr and 0<=jj<nz and inside[ii,jj]: rows.append(n); cols.append(idx[ii,jj]); vals.append(coef)
            else: rhs[n]-=coef*(chib[ii,jj] if 0<=ii<nr and 0<=jj<nz else 0.0)
        diag+=kk[i,j]*lam**2/R[i]**2      # +k lam^2 chi / r^2 moved to the left
        rows.append(n); cols.append(n); vals.append(diag); rhs[n]+=kk[i,j]*A*b
    return sp.csr_matrix((vals,(rows,cols)),shape=(N,N)),rhs
def run(A,chi0=None):
    chi=bc(A).copy() if chi0 is None else chi0.copy()
    if chi0 is None:
        chi[inside]=A*(G(RHO[inside])*(RR[inside]/RHO[inside])**2)
    for it in range(60):
        k,ok,M=solve_k(chi,A)
        if not ok[inside].all(): return None
        Mx,rhs=assemble(k,A); x=spl.spsolve(Mx,rhs)
        new=chi.copy(); new[inside]=x; d=np.max(abs(new-chi)); chi=0.5*chi+0.5*new
        if d<1e-7*max(A,1e-9): break
    k,ok,M=solve_k(chi,A)
    if not ok[inside].all(): return None
    return chi,k,M
res={}; chi=None; A=0.02; last=None
while A<1.0:
    out=run(A,chi)
    if out is None: break
    chi,k,M=out; Mi=np.nanmax(np.where(inside,M,0)); kmin=np.nanmin(np.where(inside,k,np.inf))
    th=np.linspace(0.01,np.pi-0.01,300); from scipy.interpolate import RegularGridInterpolator as RGI
    kb=RGI((R,Z),np.where(inside,k,1.0),bounds_error=False,fill_value=None)(np.c_[0.93*np.sin(th),0.93*np.cos(th)])
    Q=4*np.pi*np.trapezoid((kb-1)*0.5*np.sin(th),th)
    im=np.unravel_index(np.nanargmax(np.where(inside,M,0)),M.shape); last=(A,Mi,kmin,Q); print(f"   [max Mach at r={R[im[0]]:.3f}, z={Z[im[1]]:+.3f}]"); print(f"A={A:.3f}  max Mach={Mi:.3f}  min k={kmin:.3f}  charge Q/a={Q:+.4f}  (low-Mach -2.67*A^2 = {-2.67*A*A:+.4f})")
    A*=1.08
print("no subsonic solution beyond A ≈",A, " last:",last)
np.save("transonic_last.npy",np.array([*last]))
np.save("transonic_fields.npy",np.stack([chi,k,M]))
