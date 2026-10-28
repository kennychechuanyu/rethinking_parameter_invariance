# Code and results

R code for the formal examples, simulations and figures in *Rethinking
Parameter Invariance in Cognitive Modeling*. Run every command from the
repository root.

```sh
Rscript -e 'install.packages("renv")'    # once, if renv is not installed
Rscript --vanilla -e 'renv::restore()'   # once: install the recorded package versions
Rscript --vanilla code/run_all.R         # rebuild all seven figures (about 10 seconds)
```

The standard rebuild reruns Simulation 1 and the supplementary checks, draws
the remaining figures from the saved Monte Carlo results, and compares every
file in `code/results/` with the SHA-256 checksums in
`code/results/manifest.csv`. The figures are written to `figures/` and
`supplement_figures/`. They are visually identical to the saved ones, but each
PDF records when it was written, so Git reports them as changed;
`git checkout -- figures supplement_figures` restores the saved files.

## Software

The results were produced with R 4.5.2 on macOS (Apple silicon). `renv.lock`
records the package versions, and `code/results/software_versions.csv` the
versions used for the saved results; `run_all.R` stops if an installed package
differs from that record. The figures need an R build with Cairo support.

## Where each result comes from

| Paper | Script | Results | Seed |
|---|---|---|---:|
| Figures 1 and 2 | `R/conceptual_figures.R` | | |
| Simulation 1, Figure 3, Supplement S2 | `R/simulation1_gaussian.R` | `results/simulation1/` | 20240115 |
| Simulation 2A, Figure 4, Supplement S3 | `R/simulation2a_psychometric.R` | `results/simulation2a/` | 20240115 |
| Simulation 2B, Figure 5, Supplement S4 | `R/simulation2b_directional.R` | `results/simulation2b/` | 20260823 |
| Figure S1, Supplement S2 and S5 checks | `R/supplement_checks.R` | `results/supplement_checks/` | 20260824 |
| Figure S2, Supplement S6.1 | `R/finite_validation.R` | `results/finite_validation/` | 20260821 |

`R/common.R` holds the paths, colours, figure themes and the Wilson interval.
`R/psychometric_model.R` holds the psychometric model, the population
calculations and the model fits shared by Simulations 2A and 2B and the
finite-validation analysis. Every script can be run on its own, for example

```sh
Rscript --vanilla code/R/simulation2b_directional.R
```

`DATA_DICTIONARY.md` describes the result files and their columns.

## Rerunning the simulations

Simulation 1 and the supplementary checks take a few seconds and always run.
The other simulations are rerun only on request:

```sh
Rscript --vanilla code/run_all.R --resimulate             # Simulations 2A and 2B
Rscript --vanilla code/run_all.R --resimulate-validation  # finite-validation analysis
Rscript --vanilla code/run_all.R --resimulate-all         # all of them
```

On an Apple M5 laptop, Simulation 2A took 5.5 minutes, Simulation 2B one
minute, and the finite-validation analysis 1.5 minutes. Each simulation calls
`RNGkind("Mersenne-Twister", "Inversion", "Rejection")`, R's default, before
setting its seed. Simulation 2A draws its recovery benchmark, main grid,
sensitivity grid and focused conditions, in that order, from one stream. The
finite-validation analysis seeds training and validation separately for each
slope gap, so changing the validation samples does not change any training
fit.

A rerun overwrites the saved results and is then compared with the checksums
in the manifest. On the platform above, every file is reproduced byte for
byte. Elsewhere, last-digit floating-point differences can make the comparison
fail even when the results agree to reported precision. In a Git checkout,
`git checkout -- code/results` restores the saved files.

To check the saved results without rebuilding anything:

```sh
Rscript --vanilla code/run_all.R --verify
```

## Conventions

- `delta` is the true difference between the task optima: the difference in
  means in Simulation 1 and in slopes in Simulation 2, with
  `s_A = 2 - delta/2` and `s_B = 2 + delta/2`.
- Losses are population excess losses per observation: mean quadratic loss
  (Gaussian negative log likelihood with unit variance) in Simulation 1, and
  mean Bernoulli KL divergence per trial in Simulation 2. AIC, BIC and the
  likelihood-ratio test use the summed negative log likelihood instead.
- Directions run from source to target. `R_A_to_B` is the Task B loss when the
  Task A slope is carried into Task B and the Task B location is refitted.
- The predictive tolerance is 0.02 in Simulation 1 and 0.007 nats per trial in
  Simulation 2. The calibrated joint criterion uses one quarter of it.
- Monte Carlo intervals are 95% Wilson intervals.
