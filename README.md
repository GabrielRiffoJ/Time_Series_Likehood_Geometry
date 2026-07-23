# Likelihood geometry of moving average and autoregressive processes

Code accompanying the paper

> C. Améndola and G. Riffo, *Likelihood geometry of moving average and autoregressive processes*.
> Technische Universität Berlin.

This repository contains the scripts used to produce the tables and figures of the paper: maximum likelihood degrees of MA(*q*) and AR(*p*) models (in both the parametric and the implicit/autocovariance parametrisation), composite likelihood degrees, the symbolic proofs of Theorem 3.1 and Proposition 4.1, and the simulation studies of parameter estimation.

Scripts are organised in one folder per section of the paper, and each script is named after the table or figure it produces.

---

## Requirements

**Julia** (≥ 1.9) with

```julia
using Pkg
Pkg.add(["HomotopyContinuation", "Symbolics", "Distributions",
         "Primes", "Statistics", "StatsBase", "Plots"])
Pkg.add("Pandora")   # only needed for Figures 1-3
```

`LinearAlgebra` and `Random` are part of the Julia standard library.

**R** (≥ 4.0) — base only (`arima.sim`, `optim`); used for the `optim`-based simulation studies.

**Macaulay2** — used only for the two symbolic proofs (`.m2` files).

---

## Computational cost and how to reproduce each entry

All ML and composite likelihood degrees are computed by tracking *every* complex solution of the score equations with [`HomotopyContinuation.jl`](https://www.juliahomotopycontinuation.org/). **The cost of these computations grows very quickly** with the sample size `n` (for MA models) and with the order `q` (for the composite likelihood): the number of paths to be tracked, and hence the running time, grows by orders of magnitude from one column of a table to the next. The smallest cases run in seconds to minutes on a laptop, while the largest ones take hours. The entries reported as `-` or `--` in Tables 7 and 8 of the paper are precisely the cases that did not terminate within a 24-hour limit.

For this reason the scripts are not driven by a loop over the whole table. Instead, **a single line near the top of each script selects which entry is computed**, for example

```julia
n = 3     # sample size   (MA models)
q = 2     # model order   (composite likelihood)
```

Changing that line and re-running reproduces any other column of the corresponding table; nothing else has to be modified. Running a script as committed reproduces the entry given by the value currently set in it. Concretely,

```bash
julia 3-Model_MA_2/Table_2_and_3_Parametric_ML_degree.jl
```

computes the `n = 3` entry of Table 2; setting `n = 8` in the same file computes the `n = 8` entry.

Two exceptions: for AR(*p*) models the ML degree does not depend on the sample size, so there is one script per order *p* and no parameter to change; and the parametric composite likelihood script is written out explicitly for MA(2).

---

## Repository structure

### `3-Model_MA_2/` — Section 3, MA(2) model

| File | Output |
| --- | --- |
| `Table_2_and_3_Parametric_ML_degree.jl` | Parametric ML degree (Table 2) and its decomposition into the five solution groups (Table 3) |
| `Table_2_and_3_Implicit_ML_Degree.jl` | Implicit ML degree in autocovariance coordinates (Table 2) |
| `Table_4_Maximum_Probability_calcule.R` | Monte Carlo simulation for the probability of a local maximum at the non-invertible point (0, −1) (Table 4) |
| `Table_4_Maximum_Probability_results.R` | Collects and formats the output of the previous script |
| `Figure_3.jl` | Solution maps for `n = 3` with one sample value fixed (Figure 3) |
| `Theorem_3_1_Proof.m2` | Symbolic proof of Theorem 3.1 (Macaulay2) |

### `3.1-Images_with_Pandora/` — Figures 1 and 2

| File | Output |
| --- | --- |
| `FIgure_1_Parametric_ML_Degree.jl` | Real solution maps, parametric formulation (Figure 1) |
| `Figure_2_Implicit_ML_Degree.jl` | Real solution maps, implicit formulation (Figure 2) |

Both use [`Pandora.jl`](https://github.com/PBrdng/Pandora.jl) to count real solutions as one sample value varies, and certify the counts with `certify`.

### `4-Model_MA_3/` — Section 4, MA(3) model

| File | Output |
| --- | --- |
| `Table_7_Parametric_ML_Degree.jl` | Parametric ML degree of MA(3) (Table 7) |
| `Table_7_Implicit_ML_Degree.jl` | Implicit ML degree of MA(3) (Table 7) |
| `Proposition_4.1_Proof.m2` | Symbolic proof of Proposition 4.1 (Macaulay2) |

The parametric computation is the most expensive one in the paper: for `n = 5` it was stopped before completion, so the value 1220 reported in Table 7 is a lower bound.

### `5-Composite_Likehood_Degree/` — Section 5, composite likelihood

| File | Output |
| --- | --- |
| `Table_8_Implicit_Composite.jl` | Implicit composite likelihood degree for MA(*q*) (Table 8) |
| `Table_8_Parametric_Composite.jl` | Parametric composite likelihood degree for MA(2) — the `q = 2` entry of Table 8, together with the decomposition of its solutions into five groups discussed in Section 5 |
| `Table_8_Implicit_Composite_Example_q_2.jl` | Implicit composite likelihood degree for MA(2), worked out explicitly |

### `6-Autoregressive_Processes/` — Section 6, AR(*p*) models

| File | Output |
| --- | --- |
| `Table_11_AR1_Degree.jl` … `Table_11_AR6_Degree.jl` | ML degree of AR(*p*) for *p* = 1, …, 6 (Table 11) |

The solution counts are verified with `HomotopyContinuation.certify`, and a fixed random seed (`Random.seed!(1234)`) is set in the larger cases.

### `3.4_and_5.2_and_6.3_Simulation_Study_of_Parameter_Estimation/`

Simulation studies comparing numerical optimisation (`optim` in R) with homotopy continuation. Each study draws 500 replications from the corresponding model.

| File | Output |
| --- | --- |
| `Table_6_MA_2_estimation.R` | MA(2) estimation with `optim`, `n = 10` (Table 5) |
| `Table_5_MA_2_estimation.jl` | MA(2) estimation with `HomotopyContinuation.jl`, `n = 10` (Table 6) |
| `Table_9_MA_2_composite_estimation.R` | Composite likelihood estimation with `optim`, `n = 500` (Table 9) |
| `Table_10_MA_2_composite_estimation.jl` | Composite likelihood estimation with `HomotopyContinuation.jl`, `n = 500` (Table 10) |
| `Table_12_AR_1_estimation.R` | AR(1) estimation with `optim`, `n = 15` (Table 12) |
| `Table_13_AR_1_estimation.jl` | AR(1) estimation with `HomotopyContinuation.jl`, `n = 15` (Table 13) |

---

## Reproducibility notes

- The homotopy continuation solvers are randomised: the paths tracked depend on random start systems, so individual runs may differ in the numerical values of the solutions. The **number** of solutions is the quantity reported in the tables, and it was found to be stable across repeated runs.
- Sample values in the ML degree computations are chosen to be generic (the scripts use fixed lists of primes and their reciprocals), so that the computed degree is the generic one.

## Citation

If you use this code, please cite the paper:

```bibtex
@article{amendola-riffo-likelihood-geometry,
  author  = {Am\'endola, Carlos and Riffo, Gabriel},
  title   = {Likelihood geometry of moving average and autoregressive processes},
  year    = {2026}
}
```

## Contact

Carlos Améndola — `amendola@math.tu-berlin.de`
Gabriel Riffo — `riffo@math.tu-berlin.de`
