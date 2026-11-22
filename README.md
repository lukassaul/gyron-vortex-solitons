# Vortex solitons in a spin fluid

Working code and notes for Paper III of the gyron-fluid programme: steady, persistent vortex solutions of the transport equations of a fluid of spinning, inclined particles (gyrons), aimed at charged, spinning rings (the electron), the line vortex, quark rings and pair production.

The central picture: a particle is a knot of vorticity in the spin-pressure current (the weighted flux, whose curl is the magnetic field), with its spin carried by the weighted fluid's own mass drift, while the bulk medium flows through it nearly unperturbed.

The full reasoning record, step by step, is in [`notes/working-record.md`](notes/working-record.md) (sections 1–34). Section numbers below refer to it.

## Requirements

Python 3 with NumPy, SciPy and Matplotlib (`pip install -r requirements.txt`). Tested with NumPy 2.5, SciPy 1.18, Matplotlib 3.11; videos need `ffmpeg`.

Run every script from inside `scripts/`: several scripts load others or reuse intermediate `.npy` files written by earlier ones.

## Scripts

| Script | What it computes | Record § | Figure |
| --- | --- | --- | --- |
| `ring_exact.py` | Exact low-Mach Beltrami ring, λa = π: F-lines, B-lines, 3D helical lines | 10 | `ring_flow.png` |
| `fluxes.py` | Fluxes and winding ratio q of the λa = π ring | 12 | — |
| `hill_q.py` | Swirling Hill family, winding ratio q(λ) | 12 | — |
| `hill_plot.py` | q = 1 ring (λa = 3.633), both handednesses | 12 | `hill_q1_both.png` |
| `charge.py` | Charge at O(V²/c²) from the Bernoulli deficit | 13 | — |
| `balance2.py` | Energy coefficients of the q = 1 ring | 14 | — |
| `solve_pics2.py` | Radial balance: thin-core vs bubble pictures | 14 | — |
| `ring_profiles.py` | Streamlines, spin lines, densities, spin pressure, profiles | 14 | `ring_profiles.png` |
| `ring_full.py` | Full meridional sections and 3D streamlines, bulk and weighted | 15 | `ring_full.png` |
| `impulse.py` | Zero total impulse and ring size | 15 | — |
| `swirl_family.py` | q = 1 curve in the (λ, β) family | 20 | — |
| `pairs.py`, `pair_plot.py` | Counter-streaming population pair | 20 | `pair_ring.png` |
| `pair_spin.py` | Spin of the pair; g as a shape number | 21 | — |
| `millikan.py` | Charge at the luminal limit | 22 | — |
| `transonic.py` | Compressible (transonic) ring, choke and closed band | 24 | `transonic_ring.png` |
| `energy_spin.py`, `minimize_spin.py` | Energy at fixed spin ħ/2; predicted α | 25 | — |
| `size_map.py` | Regime with a single stable size over (ν, Γ) | 26 | `size_map.png` |
| `line_rgh.py` | Exact Rankine–Gold–Hoyle line vortex, stable core radius | 28 | — |
| `thin_ring.py` | Thin ring: tension against backreaction | 29 | — |
| `double_beltrami.py`, `db_points.py`, `db_plot.py` | Hall-like (double-Beltrami) rings: isolated q = 1, zero-impulse solutions | 33 | `double_beltrami.png` |
| `db_props.py` | Spin, charge, moment, g and size of the Hall-like rings | 34 | — |
| `pp_fields.py`, `pair_production.py` | Kinematic pair-production animation (nucleus field, bulk recoil) | — | `pair_production.mp4`, `pp2_frames.png` |

## Status

Exact: the line vortex families, the spheromak and Hill-type rings, the Beltrami pointwise charge density, the Hall-like vortex equations. Numerical: the compressible ring, the size maps, the pair and spin scans, the Hall-like solution search. Kinematic sketch only: the pair-production animation (prescribed motions, not a time evolution).

Open: the flux quantum of the W-lines, g (the Hall-like rings give g = 60–2000, not 2), α, the time evolution, mass ratios. See sections 31–34 of the record.
