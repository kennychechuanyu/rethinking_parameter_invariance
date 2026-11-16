# Benchmark-correction check (Supplement S6.2).
#
# A held-out comparison of the transfer model with a freely fitted target
# model understates the transfer excess loss, because the free model also
# estimates the focal parameter. Supplement S6.2 adds K/(2n) to correct for
# this, an approximation derived for a transferred value close to the target
# optimum. This script checks the correction for transferred slopes near the
# tolerance boundary, with Task B of Simulation 2A at delta = 0.60 as the
# target and half of its 720 trials used for fitting.
#
#   Rscript --vanilla code/R/benchmark_correction.R

source("code/R/common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})
source("code/R/psychometric_model.R")

output_dir <- file.path(results_dir, "benchmark_correction")

# Settings ----

BC_SEED <- 20261109L
BC_NREP <- 10000L                   # target datasets per transferred slope
BC_DELTA <- 0.60
BC_SLOPE_B <- 2 + BC_DELTA / 2      # Task B generating slope
BC_FIT_M <- 40L                     # fitting trials per stimulus level
BC_SCORE_M <- 40L                   # scoring trials per stimulus level
BC_FIT_N <- BC_FIT_M * length(X_LEVELS)
BC_SCORE_N <- BC_SCORE_M * length(X_LEVELS)
BC_FOCAL_DIMENSION <- 1L
BC_CORRECTION <- BC_FOCAL_DIMENSION / (2 * BC_FIT_N)
BC_MULTIPLES <- c(0.5, 1, 1.5)      # target excess losses, in tolerances

# Population quantities ----

# The slope on one side of the Task B optimum whose excess loss is a given
# multiple of the tolerance.
slope_at_multiple <- function(multiple, side, p0) {
  interval <- if (side == "shallower") {
    c(S_LO, BC_SLOPE_B)
  } else {
    c(BC_SLOPE_B, S_HI)
  }
  uniroot(
    function(slope) best_mu_pop(slope, p0)$R - multiple * EPS_PRIMARY,
    interval = interval, tol = 1e-12
  )$root
}

# Large-sample estimation cost of the transfer model's location, in units of
# 1/(2n): the variance of the per-trial score divided by the curvature of the
# loss, both at the best location for the slope. The two are equal, and the
# ratio is 1, when the slope is the target optimum.
sandwich_ratio <- function(slope, p0, step = 1e-4) {
  mu <- best_mu_pop(slope, p0)$mu
  q <- pfun(slope, mu)
  psi <- plogis(slope * (X_LEVELS - mu))
  dq <- -(1 - LAMBDA) * slope * psi * (1 - psi)
  score_variance <- mean(p0 * (1 - p0) * (dq / (q * (1 - q)))^2)
  curvature <- (
    R_T(slope, mu + step, p0) - 2 * R_T(slope, mu, p0) +
      R_T(slope, mu - step, p0)
  ) / step^2
  score_variance / curvature
}

# Simulation ----

# One target dataset: fit the transfer model (slope fixed, location free) and
# the free model to the fitting trials, then score both on the scoring trials.
# The scoring trials are balanced across stimulus levels, so given the fits
# the expected held-out difference is the difference between the two exact
# excess losses.
simulate_dataset <- function(slope, p0) {
  fit_counts <- rbinom(length(X_LEVELS), BC_FIT_M, p0)
  score_counts <- rbinom(length(X_LEVELS), BC_SCORE_M, p0)
  transfer <- fit_target_location(
    slope, fit_counts, BC_FIT_M, tolerance = 1e-10
  )
  free <- fit_sep(fit_counts, BC_FIT_M)
  transfer_risk <- R_T(slope, transfer$mu, p0)
  free_risk <- R_T(free$s, free$mu, p0)
  realized <- (
    nll_task(slope, transfer$mu, score_counts, BC_SCORE_M) -
      nll_task(free$s, free$mu, score_counts, BC_SCORE_M)
  ) / BC_SCORE_N
  c(
    expected_difference = transfer_risk - free_risk,
    realized_difference = realized,
    transfer_risk = transfer_risk,
    free_risk = free_risk,
    converged = free$conv,
    reruns = free$reruns,
    bound = free$bnd || transfer$bnd
  )
}

# Summary of the datasets for one transferred slope. Errors are mean
# differences minus the population excess loss; costs are in units of
# 1/(2n), so the K/(2n) approximation predicts 1 for the transfer model and
# 2 for the free model.
summarize_condition <- function(side, multiple, slope, p0, draws) {
  population <- best_mu_pop(slope, p0)
  excess_loss <- population$R
  n <- nrow(draws)
  expected <- draws[, "expected_difference"]
  realized <- draws[, "realized_difference"]
  pass_uncorrected <- sum(realized <= EPS_PRIMARY)
  pass_corrected <- sum(realized + BC_CORRECTION <= EPS_PRIMARY)
  tibble(
    side = side,
    tolerance_multiple = multiple,
    slope = slope,
    excess_loss = excess_loss,
    max_probability_difference_pp =
      100 * max(abs(pfun(slope, population$mu) - p0)),
    replications = n,
    expected_difference = mean(expected),
    expected_difference_mcse = sd(expected) / sqrt(n),
    realized_difference = mean(realized),
    realized_difference_mcse = sd(realized) / sqrt(n),
    error_uncorrected = mean(expected) - excess_loss,
    error_corrected = mean(expected) + BC_CORRECTION - excess_loss,
    error_uncorrected_over_tolerance =
      (mean(expected) - excess_loss) / EPS_PRIMARY,
    error_corrected_over_tolerance =
      (mean(expected) + BC_CORRECTION - excess_loss) / EPS_PRIMARY,
    cost_transfer_2n =
      mean(draws[, "transfer_risk"] - excess_loss) * 2 * BC_FIT_N,
    cost_free_2n = mean(draws[, "free_risk"]) * 2 * BC_FIT_N,
    sandwich_ratio = sandwich_ratio(slope, p0),
    pass_uncorrected_rate = pass_uncorrected / n,
    pass_uncorrected_wilson_half_width =
      wilson_interval(pass_uncorrected, n)[, "half_width"],
    pass_corrected_rate = pass_corrected / n,
    pass_corrected_wilson_half_width =
      wilson_interval(pass_corrected, n)[, "half_width"],
    reruns = sum(draws[, "reruns"])
  )
}

run_benchmark_correction <- function() {
  p0 <- pfun(BC_SLOPE_B, MU_B)
  sides <- rep(c("shallower", "steeper"), each = length(BC_MULTIPLES))
  multiples <- rep(BC_MULTIPLES, times = 2)
  conditions <- tibble(
    side = c("optimum", sides),
    tolerance_multiple = c(0, multiples),
    slope = c(
      BC_SLOPE_B,
      map2_dbl(multiples, sides, slope_at_multiple, p0 = p0)
    )
  )

  # No dataset is discarded, so every fit must converge and stay away from
  # the optimizer bounds. Otherwise the run stops.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(BC_SEED)
  results <- pmap_dfr(conditions, function(side, tolerance_multiple, slope) {
    draws <- replicate_rows(BC_NREP, function() simulate_dataset(slope, p0))
    if (any(draws[, "converged"] == 0)) {
      stop("A free fit did not converge after the restarts.")
    }
    if (any(draws[, "bound"] > 0)) {
      stop("A parameter estimate is at an optimizer bound.")
    }
    summarize_condition(side, tolerance_multiple, slope, p0, draws)
  })

  results %>%
    mutate(
      delta = BC_DELTA,
      slope_B = BC_SLOPE_B,
      mu_B = MU_B,
      lapse = LAMBDA,
      fitting_trials_per_level = BC_FIT_M,
      scoring_trials_per_level = BC_SCORE_M,
      fitting_trials = BC_FIT_N,
      focal_dimension = BC_FOCAL_DIMENSION,
      correction = BC_CORRECTION,
      tolerance = EPS_PRIMARY
    )
}

if (sys.nframe() == 0) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(
    run_benchmark_correction(),
    file.path(output_dir, "benchmark_correction.csv")
  )
  message("Saved ", file.path(output_dir, "benchmark_correction.csv"))
}
