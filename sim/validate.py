import sys, time, numpy as np
sys.path.insert(0, ".")
from sim.hall2f import Solver, double_beltrami, random_field, alphas
out = []
def log(s): print(s, flush=True); out.append(s)
# T1: ideal invariants with a random field
S = Solver(N=32, eta=1.5); Uh, Wh = random_field(S, 1, 4, amp=0.5, seed=1); d0 = S.diag(Uh, Wh); t = 0; t0 = time.time()
while t < 2.0:
    dt = S.cfl_dt(Uh, Wh); Uh, Wh = S.step(Uh, Wh, dt); t += dt
d1 = S.diag(Uh, Wh)
log(f"T1 ideal, random field (N=32, eta=1.5, t=2, {time.time()-t0:.0f}s): dE/E = {d1['E']/d0['E']-1:+.2e}, dK+/K+ = {d1['Kp']/d0['Kp']-1:+.2e}, dK-/K- = {d1['Km']/d0['Km']-1:+.2e}")
# T2: exact double-Beltrami steady states
for eta, m1, m2 in [(1.5, 1, 2), (1.05, 1, 3), (3.0, 2, 3)]:
    S = Solver(N=32, eta=eta); Uh, Wh = double_beltrami(S, m1, m2); dU, dW = S.rhs(Uh, Wh)
    nrm = lambda a: np.sqrt(np.sum(abs(a)**2)); r = (nrm(dU) + nrm(dW)) / (nrm(Uh) + nrm(Wh))
    a1, a2 = alphas(eta, m1, m2)
    log(f"T2 double Beltrami eta={eta}, mu=({m1},{m2}), alpha=({a1:+.3f},{a2:+.3f}): |rhs|/|field| = {r:.1e}")
# control: wrong alphas (not a canonical Beltrami state) should evolve
S = Solver(N=32, eta=1.5); from sim.hall2f import abc
Uh, Wh = S.fwd(0.3 * abc(S, 1) + 0.9 * 0.5 * abc(S, 2, 1, .7, .4)), S.fwd(abc(S, 1) + 0.5 * abc(S, 2, 1, .7, .4)); dU, dW = S.rhs(Uh, Wh)
nrm = lambda a: np.sqrt(np.sum(abs(a)**2)); log(f"T2 control (wrong alphas): |rhs|/|field| = {(nrm(dU)+nrm(dW))/(nrm(Uh)+nrm(Wh)):.2e}")
open("sim/validation_T12.txt", "w").write("\n".join(out) + "\n")
