# Infinite line vortex: scorecard, transport and stability

Fields depend on r only, so every target is an exact steady state of the ideal two-class fluid. Core radius r_c is selected by energy at fixed spin per length (ell = 1) and W flux (Phi_W = 1), with outer log cutoff 100 r_c. Units: xi, c.

## Shapes

| target | charge/len | spin j | axial P_z | W flux | W net circ | pitch mid | pitch edge | U sheet |
|---|---|---|---|---|---|---|---|---|
| line_RGH_x0.5 | -5.03 | +1.57 | +0 | +1.26 | +2.51 | 3.14 | 3.14 | 1.00 |
| line_RGH_x1 | -6.28 | +1.57 | +0 | +2.18 | +3.14 | 6.28 | 6.28 | 1.00 |
| line_RGH_x2 | -5.03 | +1.57 | +0 | +2.8 | +2.51 | 12.6 | 12.6 | 1.00 |
| line_RL_k2.405 | -6.52 | +1.57 | +0 | +1.36 | +3.26 | 4.22 | 0.00268 | 1.00 |
| line_RL_k3 | -4.26 | +1.57 | +0 | +0.71 | +2.13 | 2.88 | -4.81 | 1.00 |
| line_RL_k3.832 | +0.00149 | +1.57 | +0 | -0.000511 | -0.000744 | 1.47 | -9.46e+03 | 1.00 |
| line_AL_ri0.3_k2.405 | -1.96 | +0.858 | +0 | +1.36 | +3.26 | 4.22 | 0.00268 | 0.30 |
| line_AL_ri0.6_k2.405 | -3.91 | +1.21 | +0 | +1.36 | +3.26 | 4.22 | 0.00268 | 0.60 |
| line_AB_k3_e1 | -1.44 | +1.02 | +0.71 | +0.71 | +2.13 | 2.88 | -4.81 | 0.43 |
| line_AB_k2_e0.5 | -2.09 | +0.554 | +0.906 | +1.81 | +3.62 | 5.46 | 2.44 | 0.62 |
| line_hall_eta1.05_2_5 | +1.89 | -0.994 | -2.32 | +1.29 | +1.03 | 2.07 | -0.00873 | 0.94 |
| line_hall_eta1.5_2_5 | +1.14 | -0.67 | -2.39 | +1.29 | +1.03 | 2.07 | -0.00873 | 0.77 |
| line_hall_eta1.5_2_7 | -0.248 | -0.571 | -1.52 | +1.81 | +3.65 | 9.76 | 0.00319 | 0.35 |
| line_hall_eta3_2_5 | +nan | +nan | +nan | +1.29 | +1.03 | 2.07 | -0.00873 | nan |
| line_hall_eta1.5_3_6 | -2.44 | -1.7 | -1.66 | +0.21 | -0.871 | 0.173 | 0.0209 | 0.44 |

## Findings

- Every line family has exactly one stable core radius (one energy minimum at fixed spin per length and W flux, E'' > 0) in all nine media: r_c = 0.7 to 7 xi. It is subsonic for nu >= 10 or Gamma >= 1 in most families; supersonic at nu = 1, Gamma = 0.01. The Hall-like eta = 3, (2, 5) target has no real two-scale split (complex mu) and is dropped.
- Charge per length is nonzero for all but Rankine-Lundquist k = 3.832 (W_phi = 0 at the edge: no net W circulation, so no overlap at the sheet). Sign follows the W-U overlap at the edge, so both signs occur (Hall-like lines are positive with these parameters, the rest negative; reversing W flips it).
- Transport: same picture as the rings. Diffusive channels need suppression 10^18 to 10^20 below c*l (hollow annular cores 10^22-10^23: the thin hollow edge shears hardest); relaxation channels ~10^57 below nu_c. Bulk viscosity, heat conduction, inclination and axial-spin relaxation do not act.
- Stability (local criteria):
  - Kelvin-Helmholtz at the U sheet decides most cases: a sheet survives only if the jump in U is parallel to the field B = curl W there and below the tension limit. That holds for the aligned Beltrami lines (U parallel to W parallel to B) and the hollow annular cores (azimuthal jump, azimuthal B at the edge); every Rankine, Gold-Hoyle and Hall-like line fails (jump at 10-90 deg to B).
  - Centrifugal (Rayleigh, Michael) tests fail for at least one class in most families, usually the -w class (V- = U - W reverses its angular momentum profile). Passing both classes: Rankine-Gold-Hoyle x0 = 1, 2 and aligned Beltrami k = 2, eps = 0.5.
  - Only the aligned Beltrami line (k = 2, eps = 0.5) passes all axisymmetric tests, Howard-Gupta (Ri >= 0.72 > 1/4), KH and breathing.
  - Kink (m = 1) is the remaining threat: in the MHD-analogue criterion the field twist allows only periods L < L_max = 0.15 to 13 r_c. Bent into a ring of radius R (L = 2 pi R), only fat rings (R ~ r_c) escape it. Field-aligned flow can stabilize kinks (flow-shear and Alfvenic-flow stabilization), which this criterion omits.
- Status of the survey: local, MHD-analogue criteria with Alfven ratio Gamma/nu for each class. They flag instabilities but are not a normal-mode calculation. Next: a two-class normal-mode solver in (m, k) for the line, which settles the kink and the -w class centrifugal question, then the same for the rings.

## Sizing over media (nu, Gamma): r_c/xi, minima count, peak Mach

| target | nu=1,G=0.01 | nu=1,G=1 | nu=1,G=100 | nu=10,G=0.01 | nu=10,G=1 | nu=10,G=100 | nu=100,G=0.01 | nu=100,G=1 | nu=100,G=100 |
|---|---|---|---|---|---|---|---|---|---|
| line_RGH_x0.5 | 0.86 (1) M1.06 | 1.2 (1) M0.58 | 2.2 (1) M0.16 | 1.4 (1) M0.38 | 1.5 (1) M0.36 | 2.3 (1) M0.16 | 2.7 (1) M0.11 | 2.7 (1) M0.11 | 2.9 (1) M0.09 |
| line_RGH_x1 | 0.8 (1) M1.25 | 1.1 (1) M0.53 | 1.8 (1) M0.14 | 1.2 (1) M0.32 | 1.2 (1) M0.30 | 1.8 (1) M0.14 | 2.3 (1) M0.09 | 2.3 (1) M0.09 | 2.4 (1) M0.08 |
| line_RGH_x2 | 0.74 (1) M1.55 | 0.97 (1) M0.70 | 1.5 (1) M0.18 | 0.98 (1) M0.37 | 0.99 (1) M0.36 | 1.4 (1) M0.17 | 1.9 (1) M0.10 | 1.9 (1) M0.10 | 2 (1) M0.09 |
| line_RL_k2.405 | 0.91 (1) M0.89 | 1.2 (1) M0.49 | 2.3 (1) M0.14 | 1.6 (1) M0.29 | 1.6 (1) M0.28 | 2.4 (1) M0.13 | 3 (1) M0.08 | 3 (1) M0.08 | 3.1 (1) M0.07 |
| line_RL_k3 | 1 (1) M1.33 | 1.4 (1) M0.68 | 3 (1) M0.15 | 1.9 (1) M0.41 | 1.9 (1) M0.38 | 3.1 (1) M0.15 | 3.5 (1) M0.12 | 3.5 (1) M0.11 | 3.8 (1) M0.10 |
| line_RL_k3.832 | 0.7 (1) M1.88 | 1.2 (1) M0.86 | 4 (1) M0.25 | 0.46 (1) M2.17 | 1.1 (1) M0.90 | 4 (1) M0.25 | 0.33 (1) M3.02 | 1.1 (1) M0.90 | 4 (1) M0.25 |
| line_AL_ri0.3_k2.405 | 0.92 (1) M1.49 | 1.2 (1) M0.70 | 2.3 (1) M0.14 | 1.6 (1) M0.29 | 1.6 (1) M0.28 | 2.4 (1) M0.13 | 3 (1) M0.08 | 3 (1) M0.08 | 3.1 (1) M0.07 |
| line_AL_ri0.6_k2.405 | 0.92 (1) M1.08 | 1.2 (1) M0.51 | 2.3 (1) M0.14 | 1.6 (1) M0.29 | 1.6 (1) M0.28 | 2.4 (1) M0.13 | 3 (1) M0.08 | 3 (1) M0.08 | 3.1 (1) M0.07 |
| line_AB_k3_e1 | 1 (1) M1.28 | 1.4 (1) M0.69 | 3 (1) M0.15 | 1.9 (1) M0.41 | 1.9 (1) M0.38 | 3.1 (1) M0.15 | 3.5 (1) M0.12 | 3.5 (1) M0.11 | 3.8 (1) M0.10 |
| line_AB_k2_e0.5 | 0.92 (1) M1.16 | 1.1 (1) M0.61 | 2 (1) M0.13 | 1.4 (1) M0.27 | 1.4 (1) M0.26 | 2.1 (1) M0.12 | 2.7 (1) M0.08 | 2.7 (1) M0.08 | 2.8 (1) M0.07 |
| line_hall_eta1.05_2_5 | 0.99 (1) M1.79 | 1.4 (1) M0.90 | 3 (1) M0.19 | 1.3 (1) M1.04 | 1.5 (1) M0.74 | 3 (1) M0.19 | 2.5 (1) M0.29 | 2.5 (1) M0.28 | 3.3 (1) M0.16 |
| line_hall_eta1.5_2_5 | 1.2 (1) M1.47 | 1.5 (1) M0.79 | 3 (1) M0.19 | 1.3 (1) M1.01 | 1.5 (1) M0.74 | 3.1 (1) M0.19 | 2.5 (1) M0.29 | 2.5 (1) M0.28 | 3.3 (1) M0.16 |
| line_hall_eta1.5_2_7 | 1.1 (1) M2.95 | 1.3 (1) M1.92 | 2.5 (1) M0.24 | 1.5 (1) M0.29 | 1.5 (1) M0.27 | 2.6 (1) M0.09 | 2.7 (1) M0.08 | 2.7 (1) M0.08 | 3.1 (1) M0.07 |
| line_hall_eta3_2_5 | none | none | none | none | none | none | none | none | none |
| line_hall_eta1.5_3_6 | 1.9 (1) M3.49 | 3 (1) M1.44 | 6.5 (1) M0.31 | 3.5 (1) M1.07 | 3.7 (1) M0.93 | 6.6 (1) M0.30 | 6.4 (1) M0.32 | 6.4 (1) M0.31 | 7.4 (1) M0.23 |

Cell: r_c/xi (number of minima) and peak Mach. One minimum with Mach < 1 = single stable subsonic core radius.

## Stability survey (each target in its subsonic single-minimum medium closest to Mach 1)

| target | medium | S1 Rayleigh +/- | S2 Michael +/- | S3 HG Ri min +/- | S4 Suydam | S5 kink L_max/r_c | S6 KH (sheet) | S7 breathing | S8 sonic |
|---|---|---|---|---|---|---|---|---|---|
| line_RGH_x0.5 | nu=1, G=1, r_c=1.2 | True/False | True/True | 19/7.6 | stable (pressure gradient outward, shear term >= 0) | 3.17 | unstable (angle 27 deg, |dU|^2/(4 vA^2 B^2) = 1.2) | stable | subsonic |
| line_RGH_x1 | nu=1, G=1, r_c=1.1 | True/True | True/True | 39/8.8 | stable (pressure gradient outward, shear term >= 0) | 6.35 | unstable (angle 45 deg, |dU|^2/(4 vA^2 B^2) = 0.94) | stable | subsonic |
| line_RGH_x2 | nu=1, G=1, r_c=0.97 | True/True | True/True | 1.9e+02/87 | stable (pressure gradient outward, shear term >= 0) | 12.7 | unstable (angle 64 deg, |dU|^2/(4 vA^2 B^2) = 1.5) | stable | subsonic |
| line_RL_k2.405 | nu=1, G=0.01, r_c=0.91 | True/False | True/True | 3.6/0.015 | stable (pressure gradient outward, shear term >= 0) | 0.148 | unstable (angle 1 deg, |dU|^2/(4 vA^2 B^2) = 12) | stable | subsonic |
| line_RL_k3 | nu=1, G=1, r_c=1.4 | False/False | True/True | 3.5/2.7 | stable (pressure gradient outward, shear term >= 0) | 4.48 | unstable (angle 36 deg, |dU|^2/(4 vA^2 B^2) = 0.031) | stable | subsonic |
| line_RL_k3.832 | nu=100, G=1, r_c=1.1 | False/False | False/False | -43/-8.1 | stable (pressure gradient outward, shear term >= 0) | 166 | unstable (angle 89 deg, |dU|^2/(4 vA^2 B^2) = 0.00034) | stable | subsonic |
| line_AL_ri0.3_k2.405 | nu=1, G=1, r_c=1.2 | True/False | True/False | 2.8/-19 | stable (pressure gradient outward, shear term >= 0) | 0.148 | stable | stable | subsonic |
| line_AL_ri0.6_k2.405 | nu=1, G=1, r_c=1.2 | True/False | True/False | 2.8/-2.9e+02 | stable (pressure gradient outward, shear term >= 0) | 0.148 | stable | stable | subsonic |
| line_AB_k3_e1 | nu=1, G=1, r_c=1.4 | False/False | True/True | 1.3/10 | stable (pressure gradient outward, shear term >= 0) | 4.48 | stable | stable | subsonic |
| line_AB_k2_e0.5 | nu=1, G=1, r_c=1.1 | True/True | True/True | 0.72/10 | stable (pressure gradient outward, shear term >= 0) | 2.56 | stable | stable | subsonic |
| line_hall_eta1.05_2_5 | nu=1, G=1, r_c=1.4 | False/True | True/False | 0.035/-0.29 | stable (pressure gradient outward, shear term >= 0) | 5.41 | unstable (angle 14 deg, |dU|^2/(4 vA^2 B^2) = 0.4) | stable | subsonic |
| line_hall_eta1.5_2_5 | nu=1, G=1, r_c=1.5 | False/False | True/True | 0.18/0.24 | stable (pressure gradient outward, shear term >= 0) | 5.41 | unstable (angle 13 deg, |dU|^2/(4 vA^2 B^2) = 1.4) | stable | subsonic |
| line_hall_eta1.5_2_7 | nu=10, G=0.01, r_c=1.5 | False/False | False/False | -1.6/-0.28 | stable (pressure gradient outward, shear term >= 0) | 5.32 | unstable (angle 88 deg, |dU|^2/(4 vA^2 B^2) = 5.2) | stable | subsonic |
| line_hall_eta1.5_3_6 | nu=10, G=1, r_c=3.7 | False/False | False/False | -1.8e+06/-6.4e+05 | stable (pressure gradient outward, shear term >= 0) | 2.07 | unstable (angle 10 deg, |dU|^2/(4 vA^2 B^2) = 0.00024) | stable | subsonic |

Ri < 1/4 does not prove instability (Howard-Gupta is sufficient only). S5: the line is kink-stable only if its period (ring circumference) is shorter than L_max.

## Transport: suppression needed (ordinary collisional value / value allowed by persistence)

| target | bulk shear viscosity | bulk volume viscosity | U shear viscosity | W diffusivity (resistive) | spin-pressure diffusion (Fick) | weighted heat conduction | U-bulk drag (mutual friction) | inter-class drag (+w vs -w) | K-current relaxation (C_F ~ -F) | transverse spin relaxation | inclination relaxation | axial-spin relaxation |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| line_RGH_x0.5 | 0.0e+00 | 0.0e+00 | 5.9e+19 | 7.5e+18 | 7.4e+18 | 0.0e+00 | 5.0e+57 | 1.1e+57 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |
| line_RGH_x1 | 0.0e+00 | 0.0e+00 | 7.1e+19 | 4.4e+18 | 7.6e+18 | 0.0e+00 | 5.0e+57 | 1.2e+57 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |
| line_RGH_x2 | 0.0e+00 | 0.0e+00 | 8.5e+19 | 3.3e+18 | 1.5e+19 | 0.0e+00 | 5.0e+57 | 1.9e+57 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |
| line_RL_k2.405 | 0.0e+00 | 0.0e+00 | 9.7e+19 | 1.2e+19 | 1.3e+19 | 0.0e+00 | 5.0e+57 | 8.9e+56 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |
| line_RL_k3 | 0.0e+00 | 0.0e+00 | 3.9e+19 | 1.3e+19 | 1.2e+19 | 0.0e+00 | 5.0e+57 | 1.4e+57 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |
| line_RL_k3.832 | 0.0e+00 | 0.0e+00 | 6.5e+19 | 1.2e+20 | 7.8e+19 | 0.0e+00 | 5.0e+57 | 5.0e+57 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |
| line_AL_ri0.3_k2.405 | 0.0e+00 | 0.0e+00 | 3.9e+22 | 7.4e+18 | 5.8e+22 | 0.0e+00 | 5.0e+57 | 8.9e+56 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |
| line_AL_ri0.6_k2.405 | 0.0e+00 | 0.0e+00 | 4.5e+22 | 7.2e+18 | 6.8e+22 | 0.0e+00 | 5.0e+57 | 8.9e+56 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |
| line_AB_k3_e1 | 0.0e+00 | 0.0e+00 | 4.4e+19 | 1.3e+19 | 2.9e+19 | 0.0e+00 | 5.0e+57 | 1.4e+57 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |
| line_AB_k2_e0.5 | 0.0e+00 | 0.0e+00 | 3.1e+19 | 5.3e+18 | 1.0e+19 | 0.0e+00 | 5.0e+57 | 8.6e+56 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |
| line_hall_eta1.05_2_5 | 0.0e+00 | 0.0e+00 | 7.0e+19 | 5.1e+19 | 1.3e+20 | 0.0e+00 | 5.0e+57 | 4.1e+57 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |
| line_hall_eta1.5_2_5 | 0.0e+00 | 0.0e+00 | 9.6e+19 | 4.3e+19 | 6.4e+19 | 0.0e+00 | 5.0e+57 | 4.1e+57 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |
| line_hall_eta1.5_2_7 | 0.0e+00 | 0.0e+00 | 1.2e+20 | 1.1e+19 | 2.7e+20 | 0.0e+00 | 5.0e+57 | 1.1e+57 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |
| line_hall_eta1.5_3_6 | 0.0e+00 | 0.0e+00 | 2.5e+19 | 1.5e+19 | 1.5e+19 | 0.0e+00 | 5.0e+57 | 4.3e+57 | 5.0e+57 | 5.0e+57 | 0.0e+00 | 0.0e+00 |

| target | s_U | s_W | s_K | W energy inside | rotation rate (s^-1) |
|---|---|---|---|---|---|
| line_RGH_x0.5 | 8 | 1.02 | 1 | 0.21 | 2.65e+20 |
| line_RGH_x1 | 8 | 0.501 | 0.863 | 0.23 | 3.84e+20 |
| line_RGH_x2 | 8 | 0.304 | 1.41 | 0.38 | 5.58e+20 |
| line_RL_k2.405 | 8 | 1.03 | 1.07 | 0.18 | 7.13e+20 |
| line_RL_k3 | 8 | 2.57 | 2.47 | 0.29 | 1.16e+20 |
| line_RL_k3.832 | 8 | 14.7 | 9.52 | 1.00 | 3.24e+18 |
| line_AL_ri0.3_k2.405 | 5.54e+03 | 1.03 | 8.11e+03 | 0.18 | 4.55e+20 |
| line_AL_ri0.6_k2.405 | 6.52e+03 | 1.03 | 9.75e+03 | 0.18 | 3.09e+20 |
| line_AB_k3_e1 | 9 | 2.57 | 5.97 | 0.29 | 1.83e+20 |
| line_AB_k2_e0.5 | 4 | 0.688 | 1.32 | 0.17 | 8.36e+20 |
| line_hall_eta1.05_2_5 | 13.5 | 9.83 | 25.4 | 0.82 | 2.07e+20 |
| line_hall_eta1.5_2_5 | 21.8 | 9.83 | 14.6 | 0.82 | 2.24e+20 |
| line_hall_eta1.5_2_7 | 24.9 | 2.26 | 57.5 | 0.21 | 3.04e+19 |
| line_hall_eta1.5_3_6 | 35.2 | 21.4 | 20.7 | 0.86 | 2.32e+17 |