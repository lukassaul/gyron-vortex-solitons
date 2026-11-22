"""Kinematic sketch of pair production showing what the pair pushes against:
(1) the nucleus's spin-pressure gradient (its Coulomb field, held by the collisional term C_F = c^2 grad K),
(2) the bulk medium (added-mass backflow of the moving rings), and a locking-failure zone where the field is strong.
Built from the exact static pair rings; not a dynamical solution."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import LinearSegmentedColormap
exec(open("pp_fields.py").read())
BKR=LinearSegmentedColormap.from_list("bkr",["#1f5bff","#111111","#e8261c"])
nx=230; X=np.linspace(-6,6,nx); Z=np.linspace(-4,3.5,int(nx*7.5/12)); XX,ZZ=np.meshgrid(X,Z); dx=X[1]-X[0]; dz=Z[1]-Z[0]
def ring_field(x0,z0,sign,amp):
    r=np.abs(XX-x0)+1e-6; sg=np.sign(XX-x0); sg[sg==0]=1; pts=np.c_[r.ravel(),(ZZ-z0).ravel()]
    chi=iC(pts).reshape(r.shape); I=iI(pts).reshape(r.shape); Kf=iK(pts).reshape(r.shape)
    cz,cx=np.gradient(chi,dz,dx)
    return sign*amp*(-cz/r*sg), sign*amp*(sg*cx/r), sign*amp*I/r*sg, sign*amp**2*Kf
def photon(xc,amp):
    k=2*np.pi/1.2; env=np.exp(-((XX-xc)/1.3)**2-((ZZ-0)/1.6)**2); psi=amp*0.6*env*np.sin(k*(XX-xc))
    pz,px=np.gradient(psi,dz,dx); return pz,-px
xN,zN=0.0,-2.2; Qn=0.6                          # nucleus: positive spin-pressure peak (Coulomb)
def nucleus_K(zn): d=np.sqrt((XX-xN)**2+(ZZ-zn)**2+0.25**2); return Qn/(4*np.pi*d)*4
def bulk_flow(objs):
    """potential (added-mass) backflow of the bulk around moving spheres: phi = -(U a^3/2) cos/r^2"""
    ux=np.zeros_like(XX); uz=np.zeros_like(XX)
    for (x0,z0,Ux,Uz,a) in objs:
        dxr=XX-x0; dzr=ZZ-z0; r2=dxr**2+dzr**2+1e-3; r=np.sqrt(r2); out=r>a
        Udotr=(Ux*dxr+Uz*dzr)/r
        fx=(a**3/(2*r**3))*(3*Udotr*dxr/r-Ux); fz=(a**3/(2*r**3))*(3*Udotr*dzr/r-Uz)
        ux+=np.where(out,fx,Ux); uz+=np.where(out,fz,Uz)
    return ux,uz
N=160; smooth=lambda u: 0.5-0.5*np.cos(np.pi*np.clip(u,0,1))
def state(t):
    xc=-6+6*smooth(t/0.32); pa=1-smooth((t-0.30)/0.12); ra=smooth((t-0.30)/0.14)
    u=np.clip((t-0.44)/0.56,0,None); s_=(u**2/0.6 if u<0.3 else u-0.15)/0.85; d=0.25*smooth((t-0.30)/0.12)+3.0*s_
    xe,ze=-d,-0.9*s_          # electron pulled toward the nucleus
    xp,zp=+d,+0.6*s_          # positron pushed away
    zn=zN-0.25*np.clip((t-0.30)/0.7,0,None)   # nucleus recoils (absorbs the photon's surplus momentum)
    return xc,pa,ra,xe,ze,xp,zp,zn
fig,axs=plt.subplots(1,2,figsize=(18,7.2)); writer=FFMpegWriter(fps=20,bitrate=3500)
dt=1/(N-1)
with writer.saving(fig,"pair_production.mp4",dpi=100):
    for f in range(N):
        t=f*dt; xc,pa,ra,xe,ze,xp,zp,zn=state(t)
        _,_,_,xe1,ze1,xp1,zp1,zn1=state(t+dt)
        Ue=((xe1-xe)/dt*0.08,(ze1-ze)/dt*0.08); Up=((xp1-xp)/dt*0.08,(zp1-zp)/dt*0.08); Un=(0,(zn1-zn)/dt*0.08)
        Fx=np.zeros_like(XX); Fz=np.zeros_like(XX); Fy=np.zeros_like(XX); Kt=nucleus_K(zn)
        if pa>1e-3: px,pz=photon(xc,pa); Fx+=px; Fz+=pz
        if ra>1e-3:
            for x0,z0,sgn in ((xe,ze,+1),(xp,zp,-1)):
                a,b_,c_,k_=ring_field(x0,z0,sgn,ra); Fx+=a; Fz+=b_; Fy+=c_; Kt+=k_
        Ex,Ez=[-g for g in np.gradient(Kt,dz,dx)[::-1]]
        # ---- left: weighted sector
        ax=axs[0]; ax.clear(); m=0.3
        ax.pcolormesh(XX,ZZ,Kt,cmap="PuOr_r",vmin=-m,vmax=m,shading="auto")
        Emag=np.hypot(Ex,Ez); w_=np.clip(pa*ra*4,0,1)+ra*(1-smooth((t-0.40)/0.12)); fail=(Emag>np.percentile(Emag,96))&(np.hypot(XX,ZZ+0.6)<1.8)&(w_>0.2)
        ax.contourf(XX,ZZ,fail.astype(float),levels=[0.5,1.5],colors=["#ffe66b"],alpha=0.55)
        ax.contour(XX,ZZ,Kt,levels=np.linspace(-0.6,0.9,16),colors="0.55",linewidths=0.4)
        spd=np.hypot(Fx,Fz)
        if spd.max()>1e-6:
            ax.streamplot(X,Z,Fx,Fz,color=np.clip(Fy,-0.15,0.15),cmap=BKR,norm=plt.Normalize(-0.15,0.15),density=1.5,linewidth=0.4+2*np.clip(spd/0.25,0,1),arrowsize=0.8)
        ax.add_patch(plt.Circle((xN,zn),0.18,color="#c0392b")); ax.text(xN+0.3,zn-0.15,"nucleus +Ze",fontsize=9)
        if ra>1e-3 and abs(xp-xe)>1.5:
            ax.text(xe,ze+1.25,"e⁻",ha="center",fontsize=14); ax.text(xp,zp+1.25,"e⁺",ha="center",fontsize=14)
        ax.set_title("weighted sector: net flux F (streamlines), spin pressure K (colour, grey equipotentials)\nyellow: strong-gradient zone where locking fails and F-lines can reconnect",fontsize=9.5)
        # ---- right: bulk sector
        ax=axs[1]; ax.clear()
        objs=[(xN,zn,Un[0],Un[1],0.3)]
        if ra>0.05: objs+=[(xe,ze,*Ue,ra*1.0),(xp,zp,*Up,ra*1.0)]
        ux,uz=bulk_flow(objs); us=np.hypot(ux,uz)
        nb=np.ones_like(XX)
        for (x0,z0,_,_,a) in objs[1:]: nb-=0.25*np.exp(-((XX-x0)**2+(ZZ-z0)**2)/(0.6*a**2+1e-3))
        ax.pcolormesh(XX,ZZ,nb,cmap="Greys_r",vmin=0.6,vmax=1.05,shading="auto")
        if us.max()>1e-4:
            ax.streamplot(X,Z,ux,uz,color=np.log10(us+1e-5),cmap="autumn",norm=plt.Normalize(-3,-0.5),density=1.5,linewidth=0.8,arrowsize=0.8)
        for (x0,z0,Ux,Uz,a),col in zip(objs,["#c0392b","#6a3d9a","#e67e22"]):
            ax.add_patch(plt.Circle((x0,z0),max(a,0.18),fill=False,ec=col,lw=1.5))
            if np.hypot(Ux,Uz)>1e-3: ax.arrow(x0,z0,Ux*4,Uz*4,width=0.04,color=col,length_includes_head=True)
        if pa>1e-3: ax.arrow(xc-1,2.6,1.4,0,width=0.05,color="k"); ax.text(xc-1,2.85,"photon momentum (carried by the weighted wave; bulk untouched)",fontsize=8)
        ax.set_title("bulk sector: added-mass backflow around the moving rings and the recoiling nucleus\n(grey: mild bulk compression in each ring, n = 0 picture); arrows: velocities",fontsize=9.5)
        stage=("photon arrives" if t<0.30 else "conversion near the nucleus (locking fails, F-lines reconnect)" if t<0.44
               else "pair separates: e⁻ pulled toward the nucleus, e⁺ pushed away; nucleus and bulk take the recoil")
        fig.suptitle("Pair production, kinematic sketch from the exact static rings: "+stage,fontsize=12)
        for a_ in axs: a_.set_xlim(-6,6); a_.set_ylim(-4,3.5); a_.set_aspect("equal"); a_.set_xlabel("x / a"); a_.set_ylabel("z / a")
        writer.grab_frame()
        if f in (25,62,100,159): plt.savefig(f"pp2_{f:03d}.png",dpi=70)
print("done")
