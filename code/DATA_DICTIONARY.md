# Result files

Everything under `code/results/` is either saved Monte Carlo output, kept
because it is slow to regenerate, or a summary or calculation derived from it
or from the analysis code. The standard rebuild (`code/run_all.R`) rewrites
the Simulation 1, probability-scale, task-design, supplementary-check,
benchmark-correction and misspecified-lapse files and the deterministic
Simulation 2B files, and checks every file against `manifest.csv`. The
manifest lists each file with its role, dimensions, size, SHA-256 checksum,
seed and producing script; `software_versions.csv` records the software used.

## Names used in the files

The column names predate some of the paper's terminology. In particular,
`conflation` is the paper's *joint adequacy without bidirectional
transfer*.

| Name | Meaning |
|---|---|
| `delta` | True difference between the task optima (means in Simulation 1, slopes in Simulations 2A and 2B). |
| `q` | `delta` divided by the population bidirectional-transfer boundary. |
| `joint`, `raw_joint`, `J` | Joint criterion under the common tolerance: the jointly fitted curve has excess loss within the tolerance in both tasks (joint adequacy). |
| `calibrated_joint`, `J_cal` | The same criterion with one quarter of the tolerance. |
| `substitution`, `bidirectional_substitution`, `S` | Bidirectional transfer under parameter substitution: both source-to-target directions pass. |
| `A_to_B`, `B_to_A` | One transfer direction, from source to target. |
| `conflation`, `C` | Joint criterion passes while bidirectional transfer fails. |
| `R_...`, `risk` | Population excess loss per observation; in the finite-validation files, `..._validation` is its estimate from validation trials. `regret` (in `regret_curve` and the finite-validation score) means the same. |
| `..._hat` | Fitted value. `s` is a slope and `mu` a location. |
| `lrt` | Likelihood-ratio statistic: twice the difference between the shared-slope and separate-slope NLLs, set to zero when optimizer error makes it slightly negative. |
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
- `finite_sample_by_delta.csv`: the main grid of Figure 3, 2,000 replications
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
  the four candidate designs (Table S3) and the ranking of design pairs. The
  concentrated grid used for Task B is called `very_near`.
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

## Probability-scale criterion

The fits are those of the population analyses of Simulations 2A and 2B: the
joint slope and every location minimize the KL loss. `simulation2b_reversed`
swaps the generating slopes of Simulation 2B. The `pp` measures are the
largest change in a predicted response probability, in percentage points.

- `boundaries.csv`: one row per design, criterion and limit. Rows with
  `criterion = kl` use the tolerance, and rows with `criterion = pp` the limit
  in percentage points given in `limit`. `boundary_bidirectional` is the
  nearer of the two directional boundaries.
- `changes.csv`: both transfer directions at the KL boundaries and at the gap
  of 0.70 discussed with Figure 4 (`point`): the largest change in the
  target's predicted probabilities, the range of its generating probabilities,
  the number of levels with generating probabilities below .10 or above .90,
  and the excess loss per unit of mean squared probability change.
- `scale_comparison.csv`: Simulation 2B at its two KL boundaries.
  `squared_change_ratio` is the squared ratio of the largest changes in Tasks
  B and A, `kl_ratio` the ratio of the excess loss per unit of mean squared
  change in Tasks A and B, and `spread_factor` the first divided by the
  second.

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

## Benchmark correction

The target is Task B of Simulation 2A at `delta = 0.60`. Each target dataset
has 40 fitting and 40 scoring trials at each stimulus level. The transfer
model holds the slope fixed and refits the location; the free model, the
unrestricted target model of Supplement S6.2, fits both.

- `benchmark_correction.csv`: one row per transferred slope (`slope`): the
  Task B optimum and, on each side of it (`side`), the slopes whose excess
  loss is 0.5, 1 and 1.5 times the tolerance (`tolerance_multiple`), with
  10,000 datasets each. `excess_loss` is the population excess loss of the
  slope. `expected_difference` is the mean over datasets of the expected
  held-out difference given the two fits (the exact excess loss of the
  transfer fit minus that of the free fit), and `realized_difference` the
  mean difference on the scoring trials. `error_uncorrected` and
  `error_corrected` are `expected_difference` minus `excess_loss`, without and
  with the correction K/(2n) (`correction`). `cost_transfer_2n` and
  `cost_free_2n` are the mean estimation costs of the two fits in units of
  1/(2n), for which the approximation predicts 1 and 2; `sandwich_ratio` is
  the large-sample cost of the transfer fit in the same units. The
  `pass_..._rate` columns give how often the realized difference, without and
  with the correction, was within the tolerance.
  `max_probability_difference_pp` is the largest difference between the best
  curve for the slope and the generating curve.

## Misspecified lapse rates

The design is that of Simulation 2A, but Task B's responses are simulated with
lapse rate `lapse_B` while the analysis model fixes it at .02; Task A is
correctly specified. `optimum_B` is Task B's pseudo-true slope and
`divergence_B` the KL divergence per trial of its best fixed-lapse curve, from
which Task B's excess losses are measured. `max_change_B_pp` is that curve's
largest departure from the generating probabilities, in percentage points.

- `population_cases.csv`: equal generating slopes with `lapse_B` = .12, and
  generating slopes of 1.70 and 2.30 with Task B correctly specified and with
  the lapse rate at which Task B's pseudo-true slope equals 1.70. The
  `_over_tolerance` columns give losses in units of the tolerance.
- `lapse_boundaries.csv`: the Task B lapse rates at which the verdicts
  change, and the lapse rate at which Task B's pseudo-true slope equals Task
  A's slope.
- `detection.csv`: noncentral chi-square approximations to the probability
  that the likelihood-ratio test of the fixed lapse rate (`power_lrt`) and the
  deviance test against the saturated model (`power_deviance`) reject at the
  .05 level with 720 trials, and the number of trials at which the
  likelihood-ratio test reaches .80.
- `detection_monte_carlo.csv`: rejection rates in 10,000 simulated Task B
  datasets for the first two cases, next to the approximations.

## Task-design example

Tasks A and B use the concentrated grid of Simulation 2B and re-estimate their
locations; Tasks C and D fix the location at zero and place their stimuli for
target probabilities from .25 to .75 under the slope used to design them.
`D_alternative` is Task D designed with slope 1.70, and `C_broad` and
`D_broad` use the broad grid of Simulation 2A. The `pp_` columns give the
largest change in a predicted response probability, in percentage points.

- `designs.csv`: stimulus levels `x` and generating probabilities `p0` of each
  task, with the slope used to design it (`design_slope`; `NA` for Tasks A
  and B).
- `pairwise_tests.csv`: one row per pair. `joint_slope` is the shared slope,
  `R_joint` the larger of the two joint excess losses, and
  `R_first_to_second` and `R_second_to_first` the transfer losses.
  `calibrated_joint` applies one quarter of the tolerance to the joint fit,
  and `bidirectional_transfer` the tolerance to both directions.
  `pair_calibration_factor` is w1 w2 / (w1 + w2)^2, from the curvatures of
  the two tasks' excess losses at their own slopes, and
  `pair_calibrated_joint` applies that fraction of the tolerance to the joint
  fit.
- `four_task_joint.csv`: the joint fit to Tasks A to D and each task's excess
  loss.
- `task_D_regions.csv`: Task D's acceptable region for design slopes from 1.5
  to 2.8 and for the accepted joint slope, under the tolerance (`kl_`) and
  under a five-point limit on the change in any predicted probability (`pp_`),
  with the range of probabilities the design realizes (`p_min`, `p_max`).
- `task_C_regions.csv`: Task C's acceptable region under the tolerance when
  it is designed from the A-B joint slope, as in the example, and when it is
  designed from its own slope, with the range of probabilities each design
  realizes.

## Supplementary checks

- `gaussian_independent.csv`: Simulation 1 with three independent estimators
  (Supplement S2).
- `gaussian_dependence_bounds.csv`: bounds on the probability that the joint
  criterion passes while bidirectional transfer fails, from the exact pass
  probabilities of Simulation 1. They hold for any dependence between the
  estimators (Supplement S2).
- `nuisance_policy_boundaries.csv`: the A to B boundary on a one-sided grid
  with the target location fixed or refitted (Supplement S5).
