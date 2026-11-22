"""Kinematic (not dynamical) animation of pair production built from the exact static pair rings.
Electron: pair ring (lam+=5.4 member, lam-=-3.6 member, g = 2). Positron: roles swapped -> net flux and charge reversed (C image).
Meridional plane y = 0 through both ring axes (both axes along z). Lengths in units of the ring radius a."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import LinearSegmentedColormap
BKR=LinearSegmentedColormap.from_list("bkr",["#1f5bff","#111111","#e8261c"])
from scipy.special import spherical_jn as jn
from scipy.interpolate import RegularGridInterpolator as RGI
curve=np.load("q1curve.npy")
def member(l): i=np.argmin(abs(curve[:,0]-l)); return curve[i,0],curve[i,1],curve[i,2]
l1,b1,G1=member(5.4); l2,b2,G2=member(3.6); L2=-l2
A1=0.248; A2=-A1*G1/G2
def shape(lam,b): return lambda x: x*jn(1,abs(lam)*x)+b/lam**2*x**2
S1,S2=shape(l1,b1),shape(L2,b2)
# axisymmetric tables of the electron ring: net chi, net I, net K on (r,z)
rt=np.linspace(1e-3,8,500); zt=np.linspace(-8,8,1000); RT,ZT=np.meshgrid(rt,zt,indexing="ij")
rho=np.hypot(RT,ZT); s2=(RT/rho)**2; ins=rho<1
def chiI(S,lam,A,G):
    chi=A*np.where(ins,S(rho)*s2,G*s2/rho); return chi
c1=chiI(S1,l1,A1,G1); c2=chiI(S2,L2,A2,G2)
chi_net=c1-c2
I_in=l1*c1-L2*c2; D=l1*A1*G1-L2*A2*G2
I_net=np.where(ins,I_in,D*s2/rho)
# net spin pressure K: interior Bernoulli difference, exterior monopole from boundary average
def f2(chi,lam):
    gr=np.gradient(chi,rt,axis=0); gz=np.gradient(chi,zt,axis=1); return (gr**2+gz**2+(lam*chi)**2)/RT**2
kin=(b1*A1*c1-0.5*f2(c1,l1))-(b2*A2*c2-0.5*f2(c2,L2))
th=np.linspace(0.01,np.pi-0.01,200); kb=RGI((rt,zt),kin)(np.c_[0.995*np.sin(th),0.995*np.cos(th)])
k0=np.trapezoid(kb*0.5*np.sin(th),th)
K=np.where(ins,kin,k0/rho)
if k0>0: K=-K   # electron carries negative net spin pressure (charge -)
iC,iI,iK=[RGI((rt,zt),F,bounds_error=False,fill_value=0) for F in (chi_net,I_net,K)]
