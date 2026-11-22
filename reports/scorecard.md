# Scorecard (nu = 10.0, Gamma = 1.0)

Steady residuals R0, R1, R2 = relative curl of the Lamb terms of the three Euler orders (0 = exactly steady; R0, R1 need U); steady_canon = same test for the two-class fluid with field coupling (canonical vorticity), the right test for Hall-like targets). Electron targets: q = 1, impulse 0, orientation right, g = 2, mach_peak < 1, alpha = 1/137.

| name | q | steady_R2 | steady_R1 | steady_canon | impulse | spin_j | charge_Q | moment_D | orientation | rigid_fraction | a_xi | A | mach_peak | kappa | g | alpha |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| spheromak_l3.142_e0.5 | 0.825 | 3e-05 | 3e-05 | 3e-05 | +0.667 | +1.27 | -5.71 | +1 | wrong | 0.08 | 2.64 | 0.163 | 0.34 | 203 | 143 | 6.34e-05 |
| hill_l3.633_b2.077 | 1.000 | 4.8e-05 | 4.8e-05 | 4.8e-05 | +0.764 | +1.74 | -9.47 | +1.32 | wrong | 0.10 | 2.47 | 0.133 | 0.36 | 171 | 130 | 8.03e-05 |
| hall_eta1.5_3.850_8.814 | 0.864 | 0.69 | 0.48 | 0.00082 | -0.0121 | -1.96 | +0.743 | +1.03 | wrong | 0.15 | 2.1 | 0.151 | 0.82 | 115 | 966 | 8.79e-07 |
| hall_eta1.5_4.072_8.238 | 0.989 | 0.5 | 0.46 | 0.00065 | -0.0165 | -2.19 | +3.82 | +0.608 | wrong | 0.21 | 1.84 | 0.167 | 0.84 | 84.4 | 88.3 | 3.64e-05 |
| hall_eta1.5_4.722_13.454 | 1.000 | 1.5 | 0.82 | 0.0031 | -0.0164 | -1.68 | +3.67 | -0.481 | right | 0.17 | 2.31 | 0.149 | 1.06 | 142 | 102 | 1.96e-05 |
| hall_eta3_3.405_12.188 | 0.958 | 2.7 | 2.8 | 0.0021 | -0.00798 | -1.21 | -1.72 | +1.08 | right | 0.05 | 2.99 | 0.142 | 1.72 | 275 | 664 | 3.10e-06 |
| hall_eta1.5_4.722_13.454_rigidU | 1.000 | 1.5 | 1.4 | 1.6 | +0 | +1.68 | +2.69 | -0.481 | wrong | 1.00 | 2.13 | 0.167 | 1.19 | 116 | 114 | 1.77e-05 |
| spheromak_l3.142_e0_rigidU | 0.825 | 3e-05 | 0.59 | 0.62 | +0 | +1.68 | -8.5 | +1 | wrong | 1.00 | 1.82 | 0.225 | 0.47 | 81.9 | 47.8 | 5.91e-04 |
