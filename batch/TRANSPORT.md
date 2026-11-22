# Transport scorecard

Every transport coefficient, evaluated to first order on each candidate (dissipative terms on the ideal solution) at its energy-selected size. Persistence requires every rate below 2e-18 s^-1. 'Ordinary' = the collisional value (diffusivities c*l with l = 1e-32 m; rates nu_c = 1e40 s^-1). 'Suppression' = ordinary value / allowed value: the factor by which the kernel must suppress that channel.

## Summary over candidates (suppression factor needed, min to max)

| channel | min | max |
|---|---|---|
| bulk shear viscosity | 0.0e+00 | 0.0e+00 |
| bulk volume viscosity | 0.0e+00 | 0.0e+00 |
| U shear viscosity | 8.5e+18 | 1.6e+21 |
| W diffusivity (resistive) | 2.3e+18 | 4.1e+20 |
| spin-pressure diffusion (Fick) | 7.1e+18 | 1.2e+21 |
| weighted heat conduction | 0.0e+00 | 0.0e+00 |
| U-bulk drag (mutual friction) | 5.0e+57 | 5.0e+57 |
| inter-class drag (+w vs -w) | 2.9e+57 | 4.7e+57 |
| K-current relaxation (C_F ~ -F) | 5.0e+57 | 5.0e+57 |
| transverse spin relaxation | 5.0e+57 | 5.0e+57 |
| inclination relaxation | 0.0e+00 | 0.0e+00 |
| axial-spin relaxation | 0.0e+00 | 0.0e+00 |

## Channels

| key | channel | kind | acts on | effect |
|---|---|---|---|---|
| nu_b | bulk shear viscosity | diffusivity | bulk | shears bulk vorticity; ring flow is irrotational (n = 0) |
| zeta_b | bulk volume viscosity | diffusivity | bulk | damps compression; zero for a static ring |
| nu_U | U shear viscosity | diffusivity | U | unwinds the spin-carrying drift |
| nu_W | W diffusivity (resistive) | diffusivity | W | diffuses the correlation flux, field and moment |
| D_K | spin-pressure diffusion (Fick) | diffusivity | K | spreads the charge density inside the core |
| chi_T | weighted heat conduction | diffusivity | T | isothermal closure: no steady gradient to relax |
| r_Ub | U-bulk drag (mutual friction) | rate | U | locks U to the bulk; converts it to bulk vorticity |
| r_cl | inter-class drag (+w vs -w) | rate | W_in | damps the class counterflow 2W inside the core |
| r_F | K-current relaxation (C_F ~ -F) | rate | W | relaxes the static current; moment and field decay |
| r_s | transverse spin relaxation | rate | B | relaxes the locked spin S = I2 B |
| r_H | inclination relaxation | rate | none | ring carries no inclination structure at this order |
| r_par | axial-spin relaxation | rate | none | gravity sector; eps < 1e-90 |

## hall_eta1.05_4.612_13.374  (nu = 30, Gamma = 0.01, a = 1.19 xi, g = 31.2)

Shape numbers: s_U = 55.7, s_W = 49.2, s_K = 30.8 (per a^2); W energy inside 0.95; U sheet jump at the edge 0.36 of peak; rotation rate 3.31e+20 s^-1 (locking must be faster).

| channel | ordinary value | ordinary lifetime (s) | allowed value | suppression needed |
|---|---|---|---|---|
| bulk shear viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| bulk volume viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| U shear viscosity | 3.00e-24 | 1.27e-03 | 7.61e-45 | 3.9e+20 |
| W diffusivity (resistive) | 3.00e-24 | 1.44e-03 | 8.61e-45 | 3.5e+20 |
| spin-pressure diffusion (Fick) | 3.00e-24 | 2.30e-03 | 1.38e-44 | 2.2e+20 |
| weighted heat conduction | 3.00e-24 | inf | inf | 0.0e+00 |
| U-bulk drag (mutual friction) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inter-class drag (+w vs -w) | 1.00e+40 | 1.05e-40 | 2.11e-18 | 4.7e+57 |
| K-current relaxation (C_F ~ -F) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| transverse spin relaxation | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inclination relaxation | 1.00e+40 | inf | inf | 0.0e+00 |
| axial-spin relaxation | 1.00e-50 | inf | inf | 0.0e+00 |

## hall_eta1.1_4.622_13.382  (nu = 30, Gamma = 0.01, a = 1.2 xi, g = 34)

Shape numbers: s_U = 61.5, s_W = 49, s_K = 40 (per a^2); W energy inside 0.95; U sheet jump at the edge 0.35 of peak; rotation rate 3.30e+20 s^-1 (locking must be faster).

| channel | ordinary value | ordinary lifetime (s) | allowed value | suppression needed |
|---|---|---|---|---|
| bulk shear viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| bulk volume viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| U shear viscosity | 3.00e-24 | 1.17e-03 | 7.04e-45 | 4.3e+20 |
| W diffusivity (resistive) | 3.00e-24 | 1.47e-03 | 8.83e-45 | 3.4e+20 |
| spin-pressure diffusion (Fick) | 3.00e-24 | 1.80e-03 | 1.08e-44 | 2.8e+20 |
| weighted heat conduction | 3.00e-24 | inf | inf | 0.0e+00 |
| U-bulk drag (mutual friction) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inter-class drag (+w vs -w) | 1.00e+40 | 1.06e-40 | 2.11e-18 | 4.7e+57 |
| K-current relaxation (C_F ~ -F) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| transverse spin relaxation | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inclination relaxation | 1.00e+40 | inf | inf | 0.0e+00 |
| axial-spin relaxation | 1.00e-50 | inf | inf | 0.0e+00 |

## hall_eta1.25_4.656_13.407  (nu = 30, Gamma = 0.01, a = 1.25 xi, g = 43.9)

Shape numbers: s_U = 81, s_W = 48.5, s_K = 74 (per a^2); W energy inside 0.94; U sheet jump at the edge 0.35 of peak; rotation rate 3.26e+20 s^-1 (locking must be faster).

| channel | ordinary value | ordinary lifetime (s) | allowed value | suppression needed |
|---|---|---|---|---|
| bulk shear viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| bulk volume viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| U shear viscosity | 3.00e-24 | 9.52e-04 | 5.71e-45 | 5.2e+20 |
| W diffusivity (resistive) | 3.00e-24 | 1.59e-03 | 9.54e-45 | 3.1e+20 |
| spin-pressure diffusion (Fick) | 3.00e-24 | 1.04e-03 | 6.25e-45 | 4.8e+20 |
| weighted heat conduction | 3.00e-24 | inf | inf | 0.0e+00 |
| U-bulk drag (mutual friction) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inter-class drag (+w vs -w) | 1.00e+40 | 1.06e-40 | 2.12e-18 | 4.7e+57 |
| K-current relaxation (C_F ~ -F) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| transverse spin relaxation | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inclination relaxation | 1.00e+40 | inf | inf | 0.0e+00 |
| axial-spin relaxation | 1.00e-50 | inf | inf | 0.0e+00 |

## hall_eta1.5_4.722_13.453  (nu = 30, Gamma = 0.01, a = 1.33 xi, g = 68.2)

Shape numbers: s_U = 116, s_W = 47.6, s_K = 118 (per a^2); W energy inside 0.93; U sheet jump at the edge 0.28 of peak; rotation rate 3.90e+20 s^-1 (locking must be faster).

| channel | ordinary value | ordinary lifetime (s) | allowed value | suppression needed |
|---|---|---|---|---|
| bulk shear viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| bulk volume viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| U shear viscosity | 3.00e-24 | 7.55e-04 | 4.53e-45 | 6.6e+20 |
| W diffusivity (resistive) | 3.00e-24 | 1.84e-03 | 1.11e-44 | 2.7e+20 |
| spin-pressure diffusion (Fick) | 3.00e-24 | 7.42e-04 | 4.45e-45 | 6.7e+20 |
| weighted heat conduction | 3.00e-24 | inf | inf | 0.0e+00 |
| U-bulk drag (mutual friction) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inter-class drag (+w vs -w) | 1.00e+40 | 1.07e-40 | 2.15e-18 | 4.7e+57 |
| K-current relaxation (C_F ~ -F) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| transverse spin relaxation | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inclination relaxation | 1.00e+40 | inf | inf | 0.0e+00 |
| axial-spin relaxation | 1.00e-50 | inf | inf | 0.0e+00 |

## hall_eta1.5_6.301_12.548  (nu = 10, Gamma = 0.01, a = 12.4 xi, g = 1.65e+03)

Shape numbers: s_U = 130, s_W = 34.7, s_K = 108 (per a^2); W energy inside 0.58; U sheet jump at the edge 0.09 of peak; rotation rate 5.19e+19 s^-1 (locking must be faster).

| channel | ordinary value | ordinary lifetime (s) | allowed value | suppression needed |
|---|---|---|---|---|
| bulk shear viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| bulk volume viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| U shear viscosity | 3.00e-24 | 5.87e-02 | 3.52e-43 | 8.5e+18 |
| W diffusivity (resistive) | 3.00e-24 | 2.21e-01 | 1.32e-42 | 2.3e+18 |
| spin-pressure diffusion (Fick) | 3.00e-24 | 7.09e-02 | 4.25e-43 | 7.1e+18 |
| weighted heat conduction | 3.00e-24 | inf | inf | 0.0e+00 |
| U-bulk drag (mutual friction) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inter-class drag (+w vs -w) | 1.00e+40 | 1.73e-40 | 3.47e-18 | 2.9e+57 |
| K-current relaxation (C_F ~ -F) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| transverse spin relaxation | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inclination relaxation | 1.00e+40 | inf | inf | 0.0e+00 |
| axial-spin relaxation | 1.00e-50 | inf | inf | 0.0e+00 |

## hall_eta1.75_4.811_13.508  (nu = 30, Gamma = 0.01, a = 1.44 xi, g = 113)

Shape numbers: s_U = 148, s_W = 46.4, s_K = 131 (per a^2); W energy inside 0.91; U sheet jump at the edge 0.20 of peak; rotation rate 4.88e+20 s^-1 (locking must be faster).

| channel | ordinary value | ordinary lifetime (s) | allowed value | suppression needed |
|---|---|---|---|---|
| bulk shear viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| bulk volume viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| U shear viscosity | 3.00e-24 | 6.98e-04 | 4.19e-45 | 7.2e+20 |
| W diffusivity (resistive) | 3.00e-24 | 2.23e-03 | 1.33e-44 | 2.2e+20 |
| spin-pressure diffusion (Fick) | 3.00e-24 | 7.90e-04 | 4.73e-45 | 6.3e+20 |
| weighted heat conduction | 3.00e-24 | inf | inf | 0.0e+00 |
| U-bulk drag (mutual friction) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inter-class drag (+w vs -w) | 1.00e+40 | 1.09e-40 | 2.19e-18 | 4.6e+57 |
| K-current relaxation (C_F ~ -F) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| transverse spin relaxation | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inclination relaxation | 1.00e+40 | inf | inf | 0.0e+00 |
| axial-spin relaxation | 1.00e-50 | inf | inf | 0.0e+00 |

## hall_eta2_4.949_13.578  (nu = 100, Gamma = 0.01, a = 1.04 xi, g = 231)

Shape numbers: s_U = 171, s_W = 44.5, s_K = 129 (per a^2); W energy inside 0.88; U sheet jump at the edge 0.16 of peak; rotation rate 5.29e+20 s^-1 (locking must be faster).

| channel | ordinary value | ordinary lifetime (s) | allowed value | suppression needed |
|---|---|---|---|---|
| bulk shear viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| bulk volume viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| U shear viscosity | 3.00e-24 | 3.16e-04 | 1.89e-45 | 1.6e+21 |
| W diffusivity (resistive) | 3.00e-24 | 1.21e-03 | 7.28e-45 | 4.1e+20 |
| spin-pressure diffusion (Fick) | 3.00e-24 | 4.20e-04 | 2.52e-45 | 1.2e+21 |
| weighted heat conduction | 3.00e-24 | inf | inf | 0.0e+00 |
| U-bulk drag (mutual friction) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inter-class drag (+w vs -w) | 1.00e+40 | 1.13e-40 | 2.26e-18 | 4.4e+57 |
| K-current relaxation (C_F ~ -F) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| transverse spin relaxation | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inclination relaxation | 1.00e+40 | inf | inf | 0.0e+00 |
| axial-spin relaxation | 1.00e-50 | inf | inf | 0.0e+00 |

## hall_eta4_3.396_15.243  (nu = 100, Gamma = 0.01, a = 1.36 xi, g = 1.16e+03)

Shape numbers: s_U = 215, s_W = 9.97, s_K = 181 (per a^2); W energy inside 0.68; U sheet jump at the edge 0.10 of peak; rotation rate 4.55e+20 s^-1 (locking must be faster).

| channel | ordinary value | ordinary lifetime (s) | allowed value | suppression needed |
|---|---|---|---|---|
| bulk shear viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| bulk volume viscosity | 3.00e-24 | inf | inf | 0.0e+00 |
| U shear viscosity | 3.00e-24 | 4.29e-04 | 2.57e-45 | 1.2e+21 |
| W diffusivity (resistive) | 3.00e-24 | 9.24e-03 | 5.54e-44 | 5.4e+19 |
| spin-pressure diffusion (Fick) | 3.00e-24 | 5.09e-04 | 3.05e-45 | 9.8e+20 |
| weighted heat conduction | 3.00e-24 | inf | inf | 0.0e+00 |
| U-bulk drag (mutual friction) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inter-class drag (+w vs -w) | 1.00e+40 | 1.48e-40 | 2.96e-18 | 3.4e+57 |
| K-current relaxation (C_F ~ -F) | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| transverse spin relaxation | 1.00e+40 | 1.00e-40 | 2.00e-18 | 5.0e+57 |
| inclination relaxation | 1.00e+40 | inf | inf | 0.0e+00 |
| axial-spin relaxation | 1.00e-50 | inf | inf | 0.0e+00 |