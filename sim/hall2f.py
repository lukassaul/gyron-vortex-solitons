"""Dissipative two-class spin-pressure fluid, 3D periodic pseudo-spectral (box 2*pi, units with h = 1).
Classes V± = U ± W (mass drift U, correlation flux W). Field coupling by locking: A = (eta - 1) W, so the canonical
momentum of class ± is P± = V± ± (eta - 1) W and its vorticity Om± = curl U ± eta curl W. Incompressible.
    d_t U   = P[ (V+ x Om+ + V- x Om-)/2 ]           + nu_U lap U - r_Ub U
    eta d_t W = P[ (V+ x Om+ - V- x Om-)/2 ]         + nu_W lap W - r_cl W
P = Leray projection (removes all gradients: pressures, Bernoulli, electric potential).
Ideal invariants: energy E = int (U^2 + eta W^2) dV / V_box, canonical helicities K± = int P±.Om± dV / V_box."""
import numpy as np
try:
    from scipy import fft as sfft
    rfftn = lambda a: sfft.rfftn(a, axes=(-3, -2, -1), workers=-1); irfftn = lambda a, s: sfft.irfftn(a, s=s, axes=(-3, -2, -1), workers=-1)
except ImportError:
    rfftn = lambda a: np.fft.rfftn(a, axes=(-3, -2, -1)); irfftn = lambda a, s: np.fft.irfftn(a, s=s, axes=(-3, -2, -1))

class Solver:
    def __init__(self, N=32, eta=1.5, nu_U=0.0, nu_W=0.0, r_Ub=0.0, r_cl=0.0):
        self.N, self.eta = N, eta; self.co = dict(nu_U=nu_U, nu_W=nu_W, r_Ub=r_Ub, r_cl=r_cl)
        k = np.fft.fftfreq(N, 1 / N); kz = np.fft.rfftfreq(N, 1 / N)
        KX, KY, KZ = np.meshgrid(k, k, kz, indexing="ij"); self.K = np.stack([KX, KY, KZ]).astype(float)
        self.K2 = np.sum(self.K**2, 0); self.K2s = np.where(self.K2 == 0, 1, self.K2)
        kmax = N // 3; self.dealias = (abs(KX) < kmax) & (abs(KY) < kmax) & (abs(KZ) < kmax)
        x = np.arange(N) * 2 * np.pi / N; self.X = np.stack(np.meshgrid(x, x, x, indexing="ij"))
        self.s = (N, N, N)
    def fwd(self, a): return rfftn(a)
    def inv(self, a): return irfftn(a, self.s)
    def proj(self, ah):
        kd = np.sum(self.K * ah, 0) / self.K2s; return ah - self.K * kd
    def curl_h(self, ah):
        K = self.K; return 1j * np.stack([K[1] * ah[2] - K[2] * ah[1], K[2] * ah[0] - K[0] * ah[2], K[0] * ah[1] - K[1] * ah[0]])
    @staticmethod
    def cross(a, b): return np.stack([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]])
    def rhs(self, Uh, Wh):
        e = self.eta; U, W = self.inv(Uh), self.inv(Wh); wU, wW = self.inv(self.curl_h(Uh)), self.inv(self.curl_h(Wh))
        Vp, Vm = U + W, U - W; Op, Om = wU + e * wW, wU - e * wW
        Lp, Lm = self.cross(Vp, Op), self.cross(Vm, Om)
        dU = self.proj(self.fwd(0.5 * (Lp + Lm))) * self.dealias
        dW = self.proj(self.fwd(0.5 * (Lp - Lm))) * self.dealias
        c = self.co
        dU += -(c["nu_U"] * self.K2 + c["r_Ub"]) * Uh
        dW += -(c["nu_W"] * self.K2 + c["r_cl"]) * Wh
        return dU, dW / e
    def step(self, Uh, Wh, dt):
        k1 = self.rhs(Uh, Wh); k2 = self.rhs(Uh + 0.5 * dt * k1[0], Wh + 0.5 * dt * k1[1])
        k3 = self.rhs(Uh + 0.5 * dt * k2[0], Wh + 0.5 * dt * k2[1]); k4 = self.rhs(Uh + dt * k3[0], Wh + dt * k3[1])
        return (Uh + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), Wh + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]))
    def diag(self, Uh, Wh):
        e = self.eta; U, W = self.inv(Uh), self.inv(Wh); wU, wW = self.inv(self.curl_h(Uh)), self.inv(self.curl_h(Wh))
        m = lambda f: float(np.mean(f))
        Pp, Pm = U + e * W, U - e * W; Op, Om = wU + e * wW, wU - e * wW     # P± = U ± W ± (eta-1)W
        Vp, Vm = U + W, U - W
        bel = lambda V, O: m(np.sum(self.cross(V, O)**2, 0)) / max(m(np.sum(V**2, 0)) * m(np.sum(O**2, 0)), 1e-300)
        return dict(E=m(np.sum(U**2 + e * W**2, 0)), Kp=m(np.sum(Pp * Op, 0)), Km=m(np.sum(Pm * Om, 0)),
                    EU=m(np.sum(U**2, 0)), EW=m(np.sum(W**2, 0)), beltrami_p=bel(Vp, Op), beltrami_m=bel(Vm, Om))
    def cfl_dt(self, Uh, Wh, cfl=0.3):
        vmax = max(np.max(np.abs(self.inv(Uh))) + np.max(np.abs(self.inv(Wh))), 1e-12)
        dx = 2 * np.pi / self.N; d = cfl * dx / (vmax * max(self.eta, 1))
        nu = max(self.co["nu_U"], self.co["nu_W"]);
        if nu > 0: d = min(d, 0.2 * dx**2 / nu)
        return d

# ---------------- initial conditions ----------------
def abc(S, k, A=1.0, B=1.0, C=1.0, phase=(0, 0, 0)):
    """ABC Beltrami flow with curl F = k F (k integer)."""
    x, y, z = S.X[0] * k + phase[0], S.X[1] * k + phase[1], S.X[2] * k + phase[2]
    return np.stack([A * np.sin(z) + C * np.cos(y), B * np.sin(x) + A * np.cos(z), C * np.sin(y) + B * np.cos(x)])

def alphas(eta, m1, m2):
    s = 2 * eta * (m1 + m2) / (1 + eta); d = np.sqrt(s * s - 4 * eta * m1 * m2)
    return (2 * eta * m1 - s) / d, (2 * eta * m2 - s) / d

def double_beltrami(S, m1, m2, c2=0.5):
    a1, a2 = alphas(S.eta, m1, m2); W1, W2 = abc(S, m1), c2 * abc(S, m2, 1.0, 0.7, 0.4)
    return S.fwd(a1 * W1 + a2 * W2), S.fwd(W1 + W2)

def random_field(S, kmin=1, kmax=4, amp=1.0, seed=0):
    rng = np.random.default_rng(seed); out = []
    for _ in range(2):
        a = S.fwd(rng.normal(size=(3,) + S.s)); k = np.sqrt(S.K2)
        a = S.proj(a * ((k >= kmin) & (k <= kmax)) / np.maximum(k, 1)**1.5); f = S.inv(a)
        out.append(S.fwd(amp * f / np.sqrt(np.mean(np.sum(f**2, 0)))))
    return out
