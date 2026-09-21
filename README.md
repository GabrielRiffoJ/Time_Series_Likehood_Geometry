# Likelihood geometry of moving average and autoregressive processes

Code accompanying the paper

> C. Améndola and G. Riffo, *Likelihood geometry of moving average and autoregressive processes*.
> Technische Universität Berlin. arXiv:2607.04046

This repository contains the scripts used for the tables and figures of the paper: maximum likelihood degrees of MA(*q*) and AR(*p*) models (parametric and implicit/autocovariance parametrisations), composite likelihood degrees, the symbolic proofs of Theorem 3.1 and Proposition 4.1, and the simulation studies of parameter estimation.

Scripts are organised in one folder per section of the paper, and each script is named after the table or figure it produces.

---

## Requirements

**Julia** (>= 1.9). Install the packages with

```julia
using Pkg
Pkg.add(["HomotopyContinuation", "Symbolics", "Distributions", "Primes",
         "StatsBase", "Plots", "GR", "MultivariatePolynomials"])
Pkg.add(url="https://github.com/tbrysiewicz/Pandora")   # Figures 1-3 only
```

`LinearAlgebra`, `Random`, `Statistics` and `DelimitedFiles` are part of the Julia standard library.
[TODO: add `Project.toml` and `Manifest.toml` with the exact versions used, and mention them here.]

**R** (>= 4.0), base only (`arima.sim`, `optim`).

**Macaulay2**, only for the two symbolic proofs (`.m2` files).

---

## Computational cost

All ML and composite likelihood degrees are computed with [`HomotopyContinuation.jl`](https://www.juliahomotopycontinuation.org/) by tracking the solutions of the score equations. The running time grows by orders of magnitude with the sample size `n` (MA models) and the order `q` (composite likelihood). The smallest cases take seconds to minutes on a laptop; the largest take hours. Entries marked `-` in Tables 7 and 8 did not finish within 24 hours, and the value 1220 in Table 7 (`n = 5`) is a lower bound because that run was stopped.

The solution counts are computational results. Homotopy continuation does not by itself guarantee that all solutions were found, so the largest entries should be read as computed values, not proofs. [TODO: list which entries were computed more than once and which only once, e.g. MA(2) with `n = 9, 10` and MA(3).] For AR(*p*) the counts were verified with `certify`.

Scripts are not driven by a loop over the whole table. **One line near the top of each script selects the entry** (for example `n = 3` for the sample size, or `q = 2` for the order). Change that line to compute another column, for example

```bash
julia 3-Model_MA_2/Table_2_and_3_Parametric_ML_degree.jl
```

computes the `n = 3` entry of Table 2. For AR(*p*) the ML degree does not depend on `n`, so there is one script per order.

---

## Repository structure

### `3-Model_MA_2/`: MA(2)

| File | Output |
| --- | --- |
| `Table_2_and_3_Parametric_ML_degree.jl` | Parametric ML degree (Table 2) and its decomposition into five groups (Table 3) |
| `Table_2_and_3_Implicit_ML_Degree.jl` | Implicit ML degree (Table 2) |
| `Table_4_Maximum_Probability_calcule.R` | Monte Carlo for the probability of a local maximum at (0, -1) (Table 4) |
| `Table_4_Maximum_Probability_results.R` | Collects and formats the output of the previous script |
| `Figure_3.jl` | Solution maps for `n = 3` (Figure 3) |
| `Theorem_3_1_Proof.m2` | Symbolic proof of Theorem 3.1 (Macaulay2) |

### `3.1-Images_with_Pandora/`: Figures 1 and 2

| File | Output |
| --- | --- |
| `FIgure_1_Parametric_ML_Degree.jl` | Real solution maps, parametric (Figure 1) |
| `Figure_2_Implicit_ML_Degree.jl` | Real solution maps, implicit (Figure 2) |

### `4-Model_MA_3/`: MA(3)

| File | Output |
| --- | --- |
| `Table_7_Parametric_ML_Degree.jl` | Parametric ML degree (Table 7) |
| `Table_7_Implicit_ML_Degree.jl` | Implicit ML degree (Table 7) |
| `Proposition_4.1_Proof.m2` | Symbolic proof of Proposition 4.1 (Macaulay2) |

### `5-Composite_Likehood_Degree/`: composite likelihood

| File | Output |
| --- | --- |
| `Table_8_Implicit_Composite.jl` | Implicit composite likelihood degree, MA(*q*) (Table 8) |
| `Table_8_Parametric_Composite.jl` | Parametric composite likelihood degree, MA(2) (Table 8, `q = 2`) with the five solution groups of Section 5 |
| `Table_8_Implicit_Composite_Example_q_2.jl` | Implicit composite likelihood degree, MA(2), worked out explicitly |

### `6-Autoregressive_Processes/`: AR(*p*)

| File | Output |
| --- | --- |
| `Table_11_AR1_Degree.jl` ... `Table_11_AR6_Degree.jl` | ML degree of AR(*p*), *p* = 1, ..., 6 (Table 11) |

### `3.4_and_5.2_and_6.3_Simulation_Study_of_Parameter_Estimation/`

Each study uses 500 replications. **Run the R script first**, since it generates the data (with a fixed seed) and writes the sample file that the Julia script reads. Run both from this folder.

| Step | File | Output | Data file |
| --- | --- | --- | --- |
| 1 | `Table_5_MA_2_estimation.R` | MA(2) with `optim`, `n = 10` (Table 5) | writes `samples_MA2.txt` |
| 2 | `Table_6_MA_2_estimation.jl` | MA(2) with homotopy continuation, `n = 10` (Table 6) | reads `samples_MA2.txt` |
| 1 | `Table_9_MA_2_composite_estimation.R` | Composite likelihood with `optim`, `n = 500` (Table 9) | writes `samples_MA2_composite.txt` |
| 2 | `Table_10_MA_2_composite_estimation.jl` | Composite likelihood with homotopy continuation, `n = 500` (Table 10) | reads `samples_MA2_composite.txt` |
| 1 | `Table_12_AR_1_estimation.R` | AR(1) with `optim`, `n = 15` (Table 12) | writes `samples_AR1.txt` |
| 2 | `Table_13_AR_1_estimation.jl` | AR(1) with homotopy continuation, `n = 15` (Table 13) | reads `samples_AR1.txt` |

Both methods are compared on the same simulated samples.

**MA(2) parameters.** The map from `(a_0, a_1, a_2)` to the autocovariances is 8-to-1, and the sign of `a_0` gives a global symmetry. The homotopy scripts keep the real solutions with `a_0 >= 0` and, among the critical points with all roots of the characteristic polynomial outside the unit disk (the invertible representation), the one with the largest likelihood. This rule does not use the true parameter values. Replications with no such solution are counted and reported.

---

## Reproducibility notes

- The homotopy continuation solvers are randomised, so the numerical values of individual solutions may differ between runs. The **number** of solutions is the quantity reported in the tables.
- Sample values in the ML degree computations are generic (fixed lists of primes and their reciprocals), so that the degree computed is the generic one.
- The R simulation scripts call `set.seed(1234)`. The scripts for AR(5) and AR(6) also fix `Random.seed!(1234)`.

## Citation

```bibtex
@article{amendola-riffo-likelihood-geometry,
  author  = {Am\'endola, Carlos and Riffo, Gabriel},
  title   = {Likelihood geometry of moving average and autoregressive processes},
  year    = {2026},
  eprint  = {2607.04046},
  archivePrefix = {arXiv}
}
```

## Contact

Carlos Améndola: `amendola@math.tu-berlin.de`
Gabriel Riffo: `riffo@math.tu-berlin.de`
