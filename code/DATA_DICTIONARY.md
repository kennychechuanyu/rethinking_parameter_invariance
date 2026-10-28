# Result files

Everything under `code/results/` is either saved Monte Carlo output, kept
because it is slow to regenerate, or a summary or calculation derived from it
or from the analysis code. The standard rebuild (`code/run_all.R`) rewrites
the Simulation 1 and supplementary-check files and the deterministic
Simulation 2B files, and checks every file against `manifest.csv`. The
manifest lists each file with its role, dimensions, size, SHA-256 checksum,
seed and producing script; `software_versions.csv` records the software used.

## Names used in the files

The column names predate some of the paper's terminology. In particular,
`conflation` is the paper's *joint adequacy without transfer*.

| Name | Meaning |
|---|---|
| `delta` | True difference between the task optima (means in Simulation 1, slopes in Simulation 2). |
| `q` | `delta` divided by the population bidirectional-transfer boundary. |
| `joint`, `raw_joint`, `J` | Joint criterion under the common tolerance: the jointly fitted curve has excess loss within the tolerance in both tasks (joint adequacy). |
| `calibrated_joint`, `J_cal` | The same criterion with one quarter of the tolerance. |
| `substitution`, `bidirectional_substitution`, `S` | Bidirectional transfer under parameter substitution: both source-to-target directions pass. |
| `A_to_B`, `B_to_A` | One transfer direction, from source to target. |
| `conflation`, `C` | Joint criterion passes while bidirectional transfer fails. |
| `R_...`, `risk` | Population excess loss per observation; in the finite-validation files, `..._validation` is its estimate from validation trials. `regret` (in `regret_curve` and the finite-validation score) means the same. |
| `..._hat` | Fitted value. `s` is a slope and `mu` a location. |
| `lrt` | Likelihood-ratio statistic: twice the shared-slope NLL minus the separate-slope NLL, set to zero when optimizer error makes it slightly negative. |
| `aic_shared`, `bic_shared`, `M_AIC` | The criterion prefers the shared-slope model. |
| `lrt_shared`, `lrt_shared_not_rejected` | The likelihood-ratio test at the .05 level does not reject the shared-slope model. |
| `..._rate` | Proportion of replications; `replications` (or `denominator`) is its denominator. |
| `..._wilson_half_width`, `..._lo`, `..._hi`, `mc95_low`, `mc95_high` | 95% Wilson interval for a Monte Carlo proportion. |
| `conv...`, `rerun...`, `retry...` | Optimizer convergence and use of the restart grid. |
| `mcse` | Monte Carlo standard error. |
| `sJ_between` | Proportion of replications in which the joint slope lies between the two separate slopes. |
| `bound...`, `bnd...`, `near_bound` | Optimizer-bound check: an estimate within 0.01% of its optimizer interval from a bound (log scale for slopes). Distances are given on the optimizer scale, as a fraction of the interval (`scaled`), and on the parameter's own scale (`native`). |

## Simulation 1

- `simulation1_curves.csv`: one row per mean gap (61 gaps, 20,000
  replications each). Estimator bias, SD and correlation; mean losses; counts
  of every pass/fail combination (`..._k`); exact pass probabilities
  (`P_..._exact`); and simulated pass rates with Wilson intervals. The
  suffixes `_A_S` and `_B_S` mark the two task-specific (source) estimates and
  `_J` the joint estimate. `conflation_k` and `joint_only_k` count the same
  event.
- `simulation1_diagnostics.csv`: design constants, the closed-form values
  checked by the script, and summary diagnostics. Here `a` is the acceptable
  radius (r in the paper) and `s` the standard error (tau).

## Simulation 2A

- `population_results.csv`: rows with `kind = crossing` give the population
  boundaries at the three tolerances; rows with `kind = regret_curve` give the
  population losses, slopes and locations across the gap grid.
- `reference_tolerance_anchor.csv`: the largest change in response
  probability on the tolerance contour around the reference curve (Supplement
  S3).
- `finite_sample_by_delta.csv`: the main grid of Figure 4, 2,000 replications
  at each of 31 gaps. Pass rates, model-selection rates, slope recovery,
  probability-scale errors and fit checks.
- `sensitivity.csv`: rows with `analysis = grid` repeat the main grid at 360,
  720 and 1,440 trials per task and at three tolerances (1,500 replications);
  rows with `analysis = focused_grid` are the 10,000-replication conditions at
  `delta` = 0.3 to 0.7.
- `focus_delta_050.csv.gz`: every replication of the focused condition at
  `delta = 0.50`, with losses, decisions, `lrt`, parameter estimates and bound
  checks.
- `focus_event_summary.csv`, `focus_contingency_tables.csv`,
  `focus_criterion_events.csv`: event rates and 2 x 2 tables derived from that
  file, including the AIC table in Supplement S3.
- `recovery_benchmark.csv`: separate-model recovery at `delta = 0.60`.
- `finite_sample_convergence_audit.csv`: restarts, convergence and bound
  checks for the main grid (Table S1).
- `parameter_boundary_checks.csv`: one row per parameter and simulation block,
  with the flag rule on both scales.

## Simulation 2B

- `design_curvatures.csv`, `candidate_design_pairs.csv`: profile curvature of
  the four candidate designs and the ranking of design pairs (Table S3).
- `population_boundaries.csv`: the population boundaries and their local
  quadratic predictions.
- `population_curves.csv`, `profile_curvature_by_delta.csv`: population
  losses and curvatures across the gap grid.
- `reversed_orientation_boundaries.csv`: boundaries with the generating slopes
  swapped.
- `finite_sample_replications.csv.gz`: all 34,000 simulated experiments
  (2,000 at each of 17 gaps).
- `finite_sample_by_delta.csv`: summary of that file by gap; the standard
  rebuild recomputes it and checks the match.

## Finite validation

- `finite_validation_replications.csv.gz`: one row per training experiment and
  validation size (10,000 experiments at each of five gaps, three nested
  validation sizes). The fitted parameters and exact losses repeat across the
  three validation sizes of an experiment. `k_A_level_..` and `k_B_level_..`
  are the validation counts at each stimulus level.
- `specification.csv`: design, seeds and population boundaries for each gap.
- `decision_summary.csv`, `direction_summary.csv`, `event_summary.csv`,
  `contingency_tables.csv`: agreement between validation and exact decisions.
- `score_audit.csv`: the validation score against the exact excess loss.
- `training_consistency.csv`: the new training experiments against the
  Simulation 2A main grid.
- `convergence_audit.csv`, `parameter_boundary_checks.csv`: fit checks.
- `self_checks.csv`: checks of the design and the score that must pass for the
  run to finish.

## Supplementary checks

- `gaussian_independent.csv`: Simulation 1 with three independent estimators
  (Supplement S2).
- `nuisance_policy_boundaries.csv`: the A to B boundary on a one-sided grid
  with the target location fixed or refitted (Supplement S5).
