# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

This is a thesis project ("Trabajo de grado") in R comparing three different state-space model
formulations, evaluating them both statistically and computationally to determine which is most
efficient computationally and which best defines/recovers the latent states. The codebase currently
implements a Dynamic Factor Model (DFM) and the McCausland et al. (2011) precision-based simulation
smoother (MMP — "Mean and Marginal Precision" style algorithm), in the `McCausland-MMP/` directory,
based on:

- `McCausland-MMP/mccausland2011.pdf` — McCausland, Miller, Pelletier (2011), "Simulation smoothing
  for state-space models", the primary algorithmic reference (see p. 205 for the computational
  experiment parameters mirrored in the test files).
- `Mixed effects state-space models with Student-t errors.pdf` — secondary reference.
- `McCausland-MMP/formulas_modelo_factorial_McCausland_2011.md` — worked derivation of the
  formulas below (Ω blocks, `c` vector) from the article's general model down to this project's
  dynamic factor model; consult it before changing `compute_Omega_c()`.
- `librerias_testing_simulacion.md` — survey of candidate R packages for unit testing
  (`testthat`), cross-validating the smoother against external references (`KFAS`, `FKF`, `dlm`,
  `MARSS`), and running the Monte Carlo comparison study (`SimDesign`, `bench`); notes strategies
  for testing the stochastic parts (`modelo_de_simulacion()`, `mmp_draw()`) without flakiness.
  Nothing in it is implemented yet — it's a reference for future testing work.

## Environment

R project managed via RStudio (`Trabajo de grado.Rproj`). No package management file (renv/packrat)
is present — dependencies are loaded ad hoc with `library()` inside individual scripts.

## Running the code

There is no test runner or build system; scripts are plain R files meant to be sourced/run directly,
either from the RStudio console or via Rscript, from the repo root:

```
Rscript McCausland-MMP/algoritmo_simulacion_test.R
Rscript McCausland-MMP/algoritmo_estimacion_test.R
```

Each `*_test.R` file sources its corresponding implementation file (`source("algoritmo_*.R")`),
runs the algorithm, and prints diagnostic output (dimensions, matrices, coverage of the smoothed
state estimates against the simulated truth) via `cat()`/`print()` — there is no assertion-based
test framework, so "passing" is judged by inspecting the printed output.

## Architecture

The core model is a linear Gaussian state-space / dynamic factor model:

- Observation equation: `y_t = Z * alpha_t + u_t`, `u_t ~ N(0, D)`
- State equation: `alpha_{t+1} = T_mat * alpha_t + v_t`, `v_t ~ N(0, Q)`
- Initial state: `alpha_1 ~ N(a, Q1)`

with dimensions `m` (number of latent factors/states), `p` (number of observed series), and `n`
(number of time periods).

### Two-file pipeline (`McCausland-MMP/`)

- **`algoritmo_simulacion.R`** — `modelo_de_simulacion(Q, Q1, T_mat, D, a, Z, m, p, n)` simulates
  synthetic states `alpha` and observations `Y` from the model above (via Cholesky factors of
  `Q`, `Q1`, `D`), and returns them bundled with all model matrices. `validar_entradas()` performs
  strict dimension/symmetry/positive-definiteness checks on all inputs and is called at the top of
  simulation — reuse this validator rather than re-deriving checks when adding new model variants.

- **`algoritmo_estimacion.R`** — implements the McCausland precision-based simulation smoother in
  three stages, each a separate function so they can be reused/benchmarked independently:
  1. `compute_Omega_c()` — builds the block-tridiagonal precision matrix `Omega` (as a list of
     diagonal blocks `Omega_diag` plus one constant off-diagonal block `Omega_off`, since it is
     invariant across `t`) and the co-vector list `c_list`, from `Z`, `T_mat`, `Q`, `Q1`, `D`, `a`.
  2. `mmp_precompute()` — forward pass computing Cholesky factors `U_list` (upper-triangular
     `Lambda_t^T`) of the sequential precision matrices, `LiO_list` (`Lambda_t^{-1} * Omega_{t,t+1}`),
     and conditional means `m_list`; then a backward pass producing the smoothed means
     `mu_list` (`E[alpha_t | y]`). All linear solves use `backsolve()`/`crossprod()` against the
     triangular Cholesky factors rather than explicit matrix inversion — preserve this pattern for
     numerical stability and performance when modifying.
  3. `mmp_draw()` — given the precomputed factors, draws one sample of `alpha | y` via a backward
     recursion using `rnorm()` innovations.

  `estimar_estados_mmp(Y, Z, T_mat, Q, Q1, D, a, m, n, Ns = 1)` is the entry point tying the three
  stages together and returning both the smoothed mean (`mu_mat`) and `Ns` posterior draws
  (`draws`, an `m x n x Ns` array).

### `referencias/`

Earlier/alternate implementations kept for reference, not part of the active pipeline:
- `MMP.R` / `MMP_original.R` — scalar (univariate, `m=p=1`) version of the same precision-based
  smoother (`mcstates()`), useful as a simplified cross-check of the matrix algorithm above.
- `ZhuMcStates.R` — a related mixed-effects state-space simulation (`sim.messmtvp`) using
  `MASS`, `Matrix`, `LearnBayes`, `mvtnorm`.
- `Graphcyd.R` — plotting helpers (`graphvar`, `graphD`, `graphpop`, `graphind`) for visualizing
  MCMC/parameter trace output against true values.

The project is now under git, so earlier `*_v0.R` snapshots (previously kept in a
`versionamiento/` directory) have been removed — use `git log` / `git show` for version history
instead of on-disk copies.

## Conventions

- Code, comments, and console output are in Spanish; keep new code consistent with this.
- Matrices are passed/returned as R `matrix` objects (not vectors) even for `p=1`/`m=1` cases in
  the active (non-`referencias`) implementation, so a `1x1` matrix is still a matrix, not a scalar.
- Functions doc their parameters with a `# --- Parámetros de entrada ---` comment block rather than
  roxygen; follow this style for new functions in `algoritmo_*.R`.
