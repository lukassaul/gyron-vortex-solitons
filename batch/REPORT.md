# Batch report: Hall-like (double-Beltrami) ring family

Scan: eta in {1.05, 1.1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4}; for each, all isolated solutions with smooth W matching to an exterior dipole, zero impulse and winding ratio q = 1 (scanner), then the pipeline scorecard over nu in {1, 3, 10, 30, 100} x Gamma in {0.01, 0.1, 1, 10, 100} (25 media).

## Funnel

| Stage | Count |
|---|---|
| Isolated solutions found by the scanner | 40 |
| Steady in the field-coupled two-class fluid (canonical residual < 1e-2) | 40 |
| q = 1 confirmed by the pipeline (within 0.03) | 13 |
| Right orientation (moment antiparallel to spin for negative charge) | 8 |
| ...and subsonic at the energy-selected size, in at least one medium | 8 |

All scanner solutions are steady (canonical residual 3e-4 to 5e-3, grid level). Many scanner hits fail the pipeline's q check: the scanner picked a different family of tori where W has two O-points.

## Ranked candidates (best medium for each, ranked by |ln(g/2)|)

| Target | eta | mu1 a | mu2 a | nu | Gamma | a/xi | A | peak Mach | kappa | g | alpha |
|---|---|---|---|---|---|---|---|---|---|---|---|
| hall_eta1.05_4.612_13.374 | 1.05 | 4.612 | 13.374 | 30 | 0.01 | 1.19 | 0.137 | 0.99 | 31.7 | 31.2 | 4.33e-07 |
| hall_eta1.1_4.622_13.382 | 1.1 | 4.622 | 13.382 | 30 | 0.01 | 1.2 | 0.137 | 0.99 | 32.4 | 34 | 3.80e-07 |
| hall_eta1.25_4.656_13.407 | 1.25 | 4.656 | 13.407 | 30 | 0.01 | 1.25 | 0.135 | 0.97 | 35 | 43.9 | 2.56e-07 |
| hall_eta1.5_4.722_13.453 | 1.5 | 4.722 | 13.453 | 30 | 0.01 | 1.33 | 0.129 | 0.91 | 40.5 | 68.2 | 1.28e-07 |
| hall_eta1.75_4.811_13.508 | 1.75 | 4.811 | 13.508 | 30 | 0.01 | 1.44 | 0.117 | 0.90 | 48.7 | 113 | 5.55e-08 |
| hall_eta2_4.949_13.578 | 2 | 4.949 | 13.578 | 100 | 0.01 | 1.04 | 0.0642 | 0.71 | 23.4 | 231 | 3.08e-09 |
| hall_eta4_3.396_15.243 | 4 | 3.396 | 15.243 | 100 | 0.01 | 1.36 | 0.0573 | 0.80 | 42.8 | 1.16e+03 | 4.01e-10 |
| hall_eta1.5_6.301_12.548 | 1.5 | 6.301 | 12.548 | 10 | 0.01 | 12.4 | 0.112 | 0.92 | 1.17e+04 | 1.65e+03 | 4.12e-08 |

## Findings

- The lowest g among passing cases is 31.2; no candidate reaches g = 2.
- One branch carries the right orientation with q = 1 at every eta from 1.05 to 2 (mu1 a = 4.61 to 4.95, mu2 a = 13.37 to 13.58); along it g rises with eta (about 52 at eta = 1.05 to 280 at eta = 2), so the least screening (eta -> 1) is best.
- g depends weakly on the medium ratios; within the subsonic range it moves by tens of percent, not by the factor of ~25 needed.
- alpha stays between 4e-10 and 4e-7 in the best subsonic cases, far below 1/137 (larger values occur only in supersonic media).
- The obstruction is the one in working record s.34: g = 16 pi (D/Q)(a/A); the energy balance at fixed spin gives a/A ~ 10, while g = 2 needs a/A ~ 0.4 on the best branch.

## What to try next

- Add the omitted field-energy share of the mass at the alpha level and relativistic inertia (both lower a/A).
- Higher harmonics (l = 2, 3) and three or more classes, to raise Q/D.
- Extend eta below 1.05 toward 1 along the best branch; check whether g has a floor.

Figures: `batch/g_vs_eta.png` (g at nu = 10, Gamma = 1 against eta for every steady, q = 1 solution, by orientation) and `batch/best_candidate.png`.
