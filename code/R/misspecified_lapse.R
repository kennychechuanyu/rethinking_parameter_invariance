# Misspecified lapse rates (Supplement S8).
#
# The design is that of Simulation 2A, and the analysis model keeps the lapse
# rate at .02 in both tasks, but Task B's responses are simulated with a lapse
# rate lambda_B that can be larger. Task A stays correctly specified. Task B's
# optimum is then a pseudo-true slope, the slope of the best-fitting
# fixed-lapse curve, and its excess losses are measured from that curve. The
# population results are exact. How often two checks of the Task B model would
# detect the misfit is approximated with noncentral chi-square distributions
# and checked by Monte Carlo simulation.
#
#   Rscript --vanilla code/R/misspecified_lapse.R

source("code/R/common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})
source("code/R/psychometric_model.R")

output_dir <- file.path(results_dir, "misspecified_lapse")

# Settings ----

SEED <- 20261116L
N_REPLICATIONS <- 10000L                 # Task B datasets per case
CHECK_ALPHA <- 0.05                      # level of both model checks
TARGET_POWER <- 0.80                     # power for the trials needed
EQUAL_SLOPES <- c(A = 2, B = 2)
UNEQUAL_SLOPES <- c(A = 1.70, B = 2.30)  # the slope gap of 0.60
EXAMPLE_LAPSE <- 0.12
FREE_LAPSE_STARTS <- c(0.02, 0.10, 0.20)
FREE_LAPSE_HI <- 0.5                     # optimizer bound for a free lapse rate

# Population analysis ----

# Effective loss of a slope under the fixed-lapse model: the mean KL per trial
# from the generating probabilities p0, minimized over the location.
effective_loss <- function(slope, p0) best_mu_pop(slope, p0)$R

# Best fixed-lapse curve for a task. Its loss is zero when the task is
# correctly specified and is otherwise the divergence that no slope and
# location can remove.
best_curve <- function(p0) {
  slope <- optimize(function(s) effective_loss(s, p0),
                    interval = c(S_LO, S_HI), tol = 1e-10)$minimum
  fit <- best_mu_pop(slope, p0)
  list(slope = slope, location = fit$mu, loss = fit$R,
       max_change_pp = 100 * max(abs(pfun(slope, fit$mu) - p0)))
}

# Population results for generating slopes `slopes` and Task B lapse rate
# lapse_B. Each excess loss is measured from the best curve of its task.
population_case <- function(slopes, lapse_B) {
  p0_A <- pfun(slopes[["A"]], MU_A)
  p0_B <- pfun(slopes[["B"]], MU_B, lapse = lapse_B)
  best_A <- best_curve(p0_A)
  best_B <- best_curve(p0_B)
  joint_slope <- optimize(
    function(slope) effective_loss(slope, p0_A) + effective_loss(slope, p0_B),
    interval = c(S_LO, S_HI), tol = 1e-10
  )$minimum

  tibble(
    slope_A = slopes[["A"]],
    slope_B = slopes[["B"]],
    lapse_B = lapse_B,
    optimum_A = best_A$slope,
    optimum_B = best_B$slope,
    divergence_B = best_B$loss,
    max_change_B_pp = best_B$max_change_pp,
    R_A_to_B = effective_loss(best_A$slope, p0_B) - best_B$loss,
    R_B_to_A = effective_loss(best_B$slope, p0_A) - best_A$loss,
    joint_slope = joint_slope,
    R_joint_A = effective_loss(joint_slope, p0_A) - best_A$loss,
    R_joint_B = effective_loss(joint_slope, p0_B) - best_B$loss
  ) %>%
    mutate(
      # Values a little below zero are rounding error.
      across(c(divergence_B, R_A_to_B, R_B_to_A, R_joint_A, R_joint_B),
             ~ pmax(.x, 0)),
      R_transfer = pmax(R_A_to_B, R_B_to_A),
      R_joint = pmax(R_joint_A, R_joint_B)
    )
}

# Task B lapse rate at which a population measure reaches a value.
lapse_root <- function(slopes, measure, value, interval) {
  uniroot(
    function(lapse) population_case(slopes, lapse)[[measure]] - value,
    interval = interval, tol = 1e-10
  )$root
}

# Model checks ----

# Detection probabilities of two checks of the Task B model at level
# CHECK_ALPHA: the likelihood-ratio test of the fixed lapse rate against a free
# one (1 df) and the deviance test against the saturated model (nine levels
# minus two parameters, 7 df). The free-lapse and saturated models both
# contain the generating curve, so both noncentralities are 2 n D, with D the
# divergence of the best fixed-lapse curve.
detection_power <- function(divergence, trials = N_TASK) {
  noncentrality <- 2 * trials * divergence
  power <- function(df) {
    pchisq(qchisq(1 - CHECK_ALPHA, df), df, ncp = noncentrality,
           lower.tail = FALSE)
  }
  c(noncentrality = noncentrality,
    lrt = power(1),
    deviance = power(length(X_LEVELS) - 2))
}

# Trials per task at which the likelihood-ratio test detects the misfit with
# probability TARGET_POWER.
trials_for_power <- function(divergence) {
  uniroot(
    function(trials) {
      detection_power(divergence, trials)[["lrt"]] - TARGET_POWER
    },
    interval = c(10, 1e6), tol = 1e-8
  )$root
}

# Binomial negative log likelihood with a free lapse rate, without the
# constant term, as in nll_task().
nll_free_lapse <- function(par, k, m) {
  p <- clip_probability(pfun(exp(par[1]), par[2], lapse = par[3]))
  -sum(k * log(p) + (m - k) * log(1 - p))
}

# Smallest negative log likelihood with the lapse rate free: the best of
# several starts from the fixed-lapse estimates, one of them at the fixed lapse
# rate, so the result can only improve on the nested fixed-lapse fit. The
# fitted lapse rate may sit at its lower bound of zero.
min_nll_free_lapse <- function(k, m, fixed) {
  fits <- lapply(FREE_LAPSE_STARTS, function(lapse) {
    optim(
      c(log(fixed$s), fixed$mu, lapse), nll_free_lapse, k = k, m = m,
      method = "L-BFGS-B",
      lower = c(log(S_LO), MU_LO, 0),
      upper = c(log(S_HI), MU_HI, FREE_LAPSE_HI)
    )
  })
  min(vapply(fits, `[[`, numeric(1), "value"))
}

# Negative log likelihood of the saturated model, which fits every level's
# proportion exactly.
nll_saturated <- function(k, m) {
  p <- k / m
  -sum(ifelse(k > 0, k * log(p), 0) + ifelse(k < m, (m - k) * log(1 - p), 0))
}

# One Task B dataset of Simulation 2A size and the decisions of both checks.
simulate_checks <- function(slope_B, lapse_B, m = M_TRAIN) {
  k <- rbinom(length(X_LEVELS), m, pfun(slope_B, MU_B, lapse = lapse_B))
  fixed <- fit_sep(k, m)
  lrt <- 2 * (fixed$nll - min_nll_free_lapse(k, m, fixed))
  deviance <- 2 * (fixed$nll - nll_saturated(k, m))
  c(lrt_reject = lrt > qchisq(1 - CHECK_ALPHA, 1),
    deviance_reject = deviance > qchisq(1 - CHECK_ALPHA, length(X_LEVELS) - 2),
    converged = fixed$conv,
    bound = fixed$bnd)
}

# Analysis ----

run_misspecified_lapse <- function() {
  # Population cases and the lapse rates at which the verdicts change.
  compensating_lapse <- lapse_root(
    UNEQUAL_SLOPES, "optimum_B", UNEQUAL_SLOPES[["A"]], c(LAMBDA, 0.4)
  )
  boundaries <- tibble(
    slope_A = c(EQUAL_SLOPES[["A"]], EQUAL_SLOPES[["A"]],
                rep(UNEQUAL_SLOPES[["A"]], 3)),
    slope_B = c(EQUAL_SLOPES[["B"]], EQUAL_SLOPES[["B"]],
                rep(UNEQUAL_SLOPES[["B"]], 3)),
    boundary = c(
      "bidirectional transfer fails above",
      "joint adequacy fails above",
      "bidirectional transfer succeeds from",
      "Task B's pseudo-true slope equals Task A's slope at",
      "bidirectional transfer succeeds up to"
    ),
    lapse_B = c(
      lapse_root(EQUAL_SLOPES, "R_transfer", EPS_PRIMARY, c(LAMBDA, 0.4)),
      lapse_root(EQUAL_SLOPES, "R_joint", EPS_PRIMARY, c(LAMBDA, 0.4)),
      lapse_root(UNEQUAL_SLOPES, "R_transfer", EPS_PRIMARY,
                 c(LAMBDA, compensating_lapse)),
      compensating_lapse,
      lapse_root(UNEQUAL_SLOPES, "R_transfer", EPS_PRIMARY,
                 c(compensating_lapse, 0.4))
    )
  )

  cases <- bind_rows(
    population_case(EQUAL_SLOPES, EXAMPLE_LAPSE),
    population_case(UNEQUAL_SLOPES, LAMBDA),
    population_case(UNEQUAL_SLOPES, compensating_lapse)
  ) %>%
    mutate(
      case = c("equal slopes", "unequal slopes, correctly specified",
               "unequal slopes, compensating lapse"),
      .before = 1
    ) %>%
    mutate(
      tolerance = EPS_PRIMARY,
      R_A_to_B_over_tolerance = R_A_to_B / EPS_PRIMARY,
      R_B_to_A_over_tolerance = R_B_to_A / EPS_PRIMARY,
      R_joint_over_tolerance = R_joint / EPS_PRIMARY
    )

  # Detection by the model checks, approximated for the two example cases and
  # for the lower end of the range in which transfer succeeds.
  lower_end <- boundaries$lapse_B[
    boundaries$boundary == "bidirectional transfer succeeds from"
  ]
  detection <- bind_rows(
    population_case(EQUAL_SLOPES, EXAMPLE_LAPSE),
    population_case(UNEQUAL_SLOPES, compensating_lapse),
    population_case(UNEQUAL_SLOPES, lower_end)
  ) %>%
    transmute(
      case = c("equal slopes", "unequal slopes, compensating lapse",
               "unequal slopes, lower end of transfer success"),
      slope_B, lapse_B, divergence_B, max_change_B_pp,
      trials = N_TASK,
      noncentrality = map_dbl(divergence_B,
                              ~ detection_power(.x)[["noncentrality"]]),
      power_lrt = map_dbl(divergence_B, ~ detection_power(.x)[["lrt"]]),
      power_deviance = map_dbl(divergence_B,
                               ~ detection_power(.x)[["deviance"]]),
      trials_for_lrt_power_80 = map_dbl(divergence_B, trials_for_power)
    )

  # Monte Carlo check of the approximation for the two example cases. No
  # dataset is discarded, so every fixed-lapse fit must converge and stay away
  # from the optimizer bounds. Otherwise the run stops.
  checked <- filter(
    detection, case %in% c("equal slopes", "unequal slopes, compensating lapse")
  )
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(SEED)
  monte_carlo <- pmap_dfr(checked, function(case, slope_B, lapse_B, power_lrt,
                                            power_deviance, ...) {
    draws <- replicate_rows(N_REPLICATIONS, function() {
      simulate_checks(slope_B, lapse_B)
    })
    if (any(draws[, "converged"] == 0)) {
      stop("A fixed-lapse fit did not converge after the restarts.")
    }
    if (any(draws[, "bound"] > 0)) {
      stop("A parameter estimate is at an optimizer bound.")
    }
    lrt_k <- sum(draws[, "lrt_reject"])
    deviance_k <- sum(draws[, "deviance_reject"])
    tibble(
      case = case,
      slope_B = slope_B,
      lapse_B = lapse_B,
      replications = N_REPLICATIONS,
      lrt_reject_rate = lrt_k / N_REPLICATIONS,
      lrt_reject_wilson_half_width =
        wilson_interval(lrt_k, N_REPLICATIONS)[, "half_width"],
      power_lrt = power_lrt,
      deviance_reject_rate = deviance_k / N_REPLICATIONS,
      deviance_reject_wilson_half_width =
        wilson_interval(deviance_k, N_REPLICATIONS)[, "half_width"],
      power_deviance = power_deviance
    )
  })

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(cases, file.path(output_dir, "population_cases.csv"))
  write_csv(boundaries, file.path(output_dir, "lapse_boundaries.csv"))
  write_csv(detection, file.path(output_dir, "detection.csv"))
  write_csv(monte_carlo, file.path(output_dir, "detection_monte_carlo.csv"))
  message("Saved ", output_dir)
}

if (sys.nframe() == 0) {
  run_misspecified_lapse()
}
