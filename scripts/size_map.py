import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.interpolate import interp1d
rows=np.load("energy_spin_rows.npy"); A_,ew_,ef_,j_,Q_,band_,M_=rows.T
p0=5/7; sig=2*np.sqrt(2)/3
pw=np.polyfit(A_,ew_/A_**2,3); pf=np.polyfit(A_,ef_/A_**2,3); pj=np.polyfit(A_,j_/A_,3)   # smooth fits (grid kinks removed)
EW=lambda A: np.polyval(pw,A)*A**2; EF=lambda A: np.polyval(pf,A)*A**2; J=lambda A: np.polyval(pj,A)*A
Aband=0.274; As=np.linspace(A_[0],A_[-1],160); ewA,efA,jA=EW(As),EF(As),J(As)
def path(nu,Gam):
    # solve nu a^4 j = E(a)/2 for all A at once by bisection in log a
    lo=np.full(As.shape,-4.0); hi=np.full(As.shape,7.0)
    E=lambda a: nu*a**3*ewA+Gam*a*efA+p0*4*np.pi/3*a**3+4*np.pi*sig*a**2
    for _ in range(60):
        mid=(lo+hi)/2; a=10**mid; g=nu*a**4*jA-E(a)/2
        hi=np.where(g>0,mid,hi); lo=np.where(g>0,lo,mid)
    a=10**((lo+hi)/2); return np.c_[a,E(a)]
nus=np.logspace(-1,3,41); Gs=np.logspace(-3,4,36)
cls=np.zeros((len(Gs),len(nus))); asz=np.full_like(cls,np.nan)
for i,G in enumerate(Gs):
    for j,nu in enumerate(nus):
        P=path(nu,G); E=P[:,1]
        mins=[k for k in range(1,len(As)-1) if E[k]<E[k-1] and E[k]<E[k+1]]
        good=[k for k in mins if As[k]<Aband]
        if len(good)==1 and len(mins)==1: cls[i,j]=2; asz[i,j]=P[good[0],0]
        elif len(mins)==0: cls[i,j]=0 if np.argmin(E)==len(As)-1 else 1
        elif len(good)==0: cls[i,j]=3
        else: cls[i,j]=4
fig,ax=plt.subplots(1,2,figsize=(15,6))
cmap=matplotlib.colors.ListedColormap(["#d9d9d9","#9ecae1","#31a354","#fdae6b","#de2d26"])
ax[0].pcolormesh(nus,Gs,cls,cmap=cmap,vmin=-0.5,vmax=4.5,shading="nearest")
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("ν (weighted inertia / bulk)"); ax[0].set_ylabel("Γ (field stiffness)")
for c,l in zip(cmap.colors,["runs to the choke","runs to A → 0","ONE stable size (closed band)","minimum only after band lost","several minima"]):
    ax[0].plot([],[],"s",color=c,ms=10,label=l)
ax[0].legend(loc="lower left",fontsize=8); ax[0].set_title("energy at fixed spin ħ/2: number of stable sizes")
pc=ax[1].pcolormesh(nus,Gs,np.log10(asz),cmap="viridis",shading="nearest"); fig.colorbar(pc,ax=ax[1],label="log₁₀ a/ξ")
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("ν"); ax[1].set_ylabel("Γ"); ax[1].set_title("ring size where unique")
plt.tight_layout(); plt.savefig("size_map.png",dpi=110)
print("fit residual rms (e_w/A^2, e_f/A^2, j/A):",np.std(ew_/A_**2-np.polyval(pw,A_)).round(3),np.std(ef_/A_**2-np.polyval(pf,A_)).round(3),np.std(j_/A_-np.polyval(pj,A_)).round(4))
u=cls==2; print("fraction with one stable size:",u.mean().round(3))
print("nu range with unique size:",nus[u.any(0)].min().round(2),"to",nus[u.any(0)].max().round(1))
print("size range a/xi:",np.nanmin(asz).round(3),"to",np.nanmax(asz).round(2))
for nu_t in (3,10,100):
    j=np.argmin(abs(nus-nu_t)); print(f" nu≈{nus[j]:.1f}: Γ with unique size {Gs[cls[:,j]==2].min() if (cls[:,j]==2).any() else None} – {Gs[cls[:,j]==2].max() if (cls[:,j]==2).any() else None}, a/ξ {np.nanmin(asz[:,j]):.2f}–{np.nanmax(asz[:,j]):.2f}")
print("classes present:",{int(k):int((cls==k).sum()) for k in np.unique(cls)})
