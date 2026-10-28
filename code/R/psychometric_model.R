# Psychometric model shared by Simulations 2A and 2B and the finite-validation
# analysis: design settings, population calculations, and model fits.
#
# Sourcing this file only defines objects. It draws no random numbers and
# writes no files.

# Settings ----

X_LEVELS <- seq(-2, 2, by = 0.5)  # nine stimulus levels
LAMBDA <- 0.02                    # lapse rate, fixed
MU_A <- -0.15                     # task locations
MU_B <- 0.15
DELTAS <- seq(0, 1.50, by = 0.05) # slope gaps; s_A = 2 - delta/2, s_B = 2 + delta/2

M_TRAIN <- 80                        # trials per stimulus level
N_TASK <- M_TRAIN * length(X_LEVELS) # 720 trials per task
N_TOTAL <- 2L * N_TASK               # 1,440 trials in the two tasks

# Shared- versus separate-slope model selection. With three and four free
# parameters, each criterion is a cutoff on the likelihood-ratio statistic
# (natural logs throughout).
K_SHARED <- 3L
K_SEPARATE <- 4L
AIC_LRT_CUTOFF <- 2 * (K_SEPARATE - K_SHARED)
LRT_05_CUTOFF <- qchisq(0.95, df = K_SEPARATE - K_SHARED)
BIC_LRT_CUTOFF <- log(N_TOTAL) * (K_SEPARATE - K_SHARED)

# The separate model nests the shared model, so a negative likelihood-ratio
# statistic can only be optimizer error. One below -LRT_NESTING_TOL stops the
# run.
LRT_NESTING_TOL <- 1e-6

EPS_PRIMARY <- 0.007                # predictive tolerance, nats per trial
EPS_GRID <- c(0.0035, 0.007, 0.014) # tolerance sensitivity

# The factor of 4 aligns the joint and transfer boundaries exactly for
# equal-curvature quadratic losses. It is applied here without retuning.
JOINT_CALIBRATION_FACTOR <- 4
EPS_JOINT_CALIBRATED <- EPS_PRIMARY / JOINT_CALIBRATION_FACTOR

# Optimizer bounds and restart grid (Supplement S3).
S_LO <- 0.1
S_HI <- 8
MU_LO <- -1.5
MU_HI <- 1.5
LOGS_STARTS <- log(c(0.5, 1, 2, 4))
MU_STARTS <- c(-0.5, 0, 0.5)

# A fitted parameter is flagged as near a bound when its distance from the
# nearer bound is at most 0.01% of the optimizer interval. Slopes are
# optimized, and checked, on the log scale.
BOUND_SCALED_TOL <- 1e-4
LOG_S_LO <- log(S_LO)
LOG_S_HI <- log(S_HI)

# The nine free parameters checked in every replication.
BOUND_PARAMETER_IDS <- c(
  "sep_A_s", "sep_A_mu", "sep_B_s", "sep_B_mu",
  "joint_s", "joint_mu_A", "joint_mu_B",
  "substitution_A_to_B_mu", "substitution_B_to_A_mu"
)
BOUND_FLAG_COLUMNS <- paste0("bound_", BOUND_PARAMETER_IDS)
BOUND_SCALED_DISTANCE_COLUMNS <- paste0(
  "bound_scaled_distance_", BOUND_PARAMETER_IDS
)
BOUND_OPTIMIZER_DISTANCE_COLUMNS <- paste0(
  "bound_optimizer_distance_", BOUND_PARAMETER_IDS
)
BOUND_NATIVE_DISTANCE_COLUMNS <- paste0(
  "bound_native_distance_", BOUND_PARAMETER_IDS
)

# Model and loss ----

# Response probability at stimulus levels x: a logistic curve with slope s and
# location mu, mixed with lapses.
pfun <- function(s, mu, x = X_LEVELS, lapse = LAMBDA) {
  lapse / 2 + (1 - lapse) * plogis(s * (x - mu))
}

clip_probability <- function(p) pmin(pmax(p, 1e-12), 1 - 1e-12)

# Elementwise KL(Bernoulli(p0) || Bernoulli(q)).
kl_bern <- function(p0, q) {
  q <- clip_probability(q)
  p0 <- clip_probability(p0)
  p0 * log(p0 / q) + (1 - p0) * log((1 - p0) / (1 - q))
}

# Population excess loss of the curve (s, mu): the mean KL per trial from the
# generating probabilities p0. The generating curve is in the model family, so
# the minimum over (s, mu) is zero.
R_T <- function(s, mu, p0, x = X_LEVELS) {
  mean(kl_bern(p0, pfun(s, mu, x)))
}

# Effective loss of slope s, min over mu of R_T(s, mu), and the location that
# attains it.
best_mu_pop <- function(s, p0, x = X_LEVELS) {
  o <- optimize(
    function(mu) R_T(s, mu, p0, x),
    interval = c(MU_LO, MU_HI), tol = 1e-10
  )
  list(mu = o$minimum, R = o$objective)
}

# Binomial negative log likelihood, without the constant term.
nll_task <- function(s, mu, k, m, x = X_LEVELS) {
  p <- clip_probability(pfun(s, mu, x))
  -sum(k * log(p) + (m - k) * log(1 - p))
}

# Distance of a fitted parameter from the nearer bound: on the optimizer scale
# (raw and as a fraction of the interval) and on the parameter's own scale.
bound_diagnostic <- function(value, lower, upper,
                             native_value = value,
                             native_lower = lower,
                             native_upper = upper) {
  optimizer_distance <- max(0, min(value - lower, upper - value))
  native_distance <- max(
    0, min(native_value - native_lower, native_upper - native_value)
  )
  scaled_distance <- optimizer_distance / (upper - lower)
  list(
    flag = scaled_distance <= BOUND_SCALED_TOL,
    scaled_distance = scaled_distance,
    optimizer_distance = optimizer_distance,
    native_distance = native_distance
  )
}

# Population analysis ----

# Population results at slope gap delta. The joint value shares the slope and
# gives each task its own best location. With equal designs and trial counts,
# the sum of the two mean losses has the same minimizer as the expected joint
# NLL.
population_at_delta <- function(delta) {
  slope_A <- 2 - delta / 2
  slope_B <- 2 + delta / 2
  generating_A <- pfun(slope_A, MU_A)
  generating_B <- pfun(slope_B, MU_B)

  joint_objective <- function(slope) {
    best_mu_pop(slope, generating_A)$R +
      best_mu_pop(slope, generating_B)$R
  }
  joint_solution <- optimize(
    joint_objective,
    interval = c(S_LO, S_HI),
    tol = 1e-10
  )
  joint_slope <- joint_solution$minimum
  joint_A <- best_mu_pop(joint_slope, generating_A)
  joint_B <- best_mu_pop(joint_slope, generating_B)
  joint_risk <- max(joint_A$R, joint_B$R)

  transfer_A_to_B <- best_mu_pop(slope_A, generating_B)
  transfer_B_to_A <- best_mu_pop(slope_B, generating_A)
  substitution_risk <- max(transfer_A_to_B$R, transfer_B_to_A$R)

  # Largest difference in response probability across the nine levels, in
  # percentage points.
  joint_max_probability_difference <- max(
    max(abs(pfun(joint_slope, joint_A$mu) - generating_A)),
    max(abs(pfun(joint_slope, joint_B$mu) - generating_B))
  ) * 100
  substitution_max_probability_difference <- max(
    max(abs(pfun(slope_A, transfer_A_to_B$mu) - generating_B)),
    max(abs(pfun(slope_B, transfer_B_to_A$mu) - generating_A))
  ) * 100

  tibble(
    delta = delta,
    sA = slope_A,
    sB = slope_B,
    sJ = joint_slope,
    mu_joint_A = joint_A$mu,
    mu_joint_B = joint_B$mu,
    mu_A_to_B = transfer_A_to_B$mu,
    mu_B_to_A = transfer_B_to_A$mu,
    R_J = joint_risk,
    R_S = substitution_risk,
    R_joint_A = joint_A$R,
    R_joint_B = joint_B$R,
    R_B_to_A = transfer_B_to_A$R,
    R_A_to_B = transfer_A_to_B$R,
    joint_max_probability_difference_pp = joint_max_probability_difference,
    substitution_max_probability_difference_pp =
      substitution_max_probability_difference
  )
}

# Population boundary: the gap at which a criterion's excess loss reaches eps,
# found with uniroot().
population_crossing <- function(eps, risk_name) {
  objective <- function(delta) population_at_delta(delta)[[risk_name]] - eps
  uniroot(objective, interval = range(DELTAS), tol = 1e-10)$root
}

run_population <- function() {
  pop <- map_dfr(DELTAS, population_at_delta)
  cross <- map_dfr(EPS_GRID, function(eps) {
    substitution_crossing <- population_crossing(eps, "R_S")
    joint_crossing <- population_crossing(eps, "R_J")
    calibrated_joint_crossing <- population_crossing(
      eps / JOINT_CALIBRATION_FACTOR,
      "R_J"
    )
    tibble(
      eps = eps,
      substitution_crossing = substitution_crossing,
      joint_crossing = joint_crossing,
      calibrated_joint_crossing = calibrated_joint_crossing,
      ratio = joint_crossing / substitution_crossing,
      calibrated_minus_substitution =
        calibrated_joint_crossing - substitution_crossing,
      calibration_factor = JOINT_CALIBRATION_FACTOR
    )
  })
  list(pop = pop, cross = cross)
}

# Behavioural meaning of the tolerance (Supplement S3): on the reference curve
# s = 2, mu = 0, move the slope down and up until the mean KL reaches
# EPS_PRIMARY, and record the largest change in response probability.
reference_tolerance_anchor <- function() {
  reference_slope <- 2
  reference_location <- 0
  generating <- pfun(reference_slope, reference_location)
  risk_at_slope <- function(slope) {
    R_T(slope, reference_location, generating)
  }
  lower_slope <- uniroot(
    function(slope) risk_at_slope(slope) - EPS_PRIMARY,
    interval = c(S_LO, reference_slope), tol = 1e-12
  )$root
  upper_slope <- uniroot(
    function(slope) risk_at_slope(slope) - EPS_PRIMARY,
    interval = c(reference_slope, S_HI), tol = 1e-12
  )$root

  slopes <- c(lower_slope, upper_slope)
  maximum_probability_change <- vapply(slopes, function(slope) {
    max(abs(pfun(slope, reference_location) - generating)) * 100
  }, numeric(1))
  risks <- vapply(slopes, risk_at_slope, numeric(1))

  tibble(
    side = c("lower slope", "higher slope"),
    reference_slope = reference_slope,
    contour_slope = slopes,
    reference_location = reference_location,
    tolerance = EPS_PRIMARY,
    mean_kl = risks,
    maximum_probability_change_percentage_points =
      maximum_probability_change
  )
}

# One table with the boundary rows (kind = "crossing") and the population
# curves (kind = "regret_curve").
population_output_table <- function(pop, cross) {
  pop_curve <- pop %>%
    transmute(
      kind = "regret_curve", eps = NA_real_,
      substitution = NA_real_, joint = NA_real_, joint_calibrated = NA_real_,
      ratio = NA_real_, calibrated_minus_substitution = NA_real_,
      calibration_factor = JOINT_CALIBRATION_FACTOR,
      delta, s_A = sA, s_B = sB, s_joint = sJ,
      R_joint = pmax(R_J, 0),
      R_substitution = pmax(R_S, 0),
      R_joint_A = pmax(R_joint_A, 0),
      R_joint_B = pmax(R_joint_B, 0),
      R_B_to_A = pmax(R_B_to_A, 0),
      R_A_to_B = pmax(R_A_to_B, 0),
      joint_max_probability_difference_pp,
      substitution_max_probability_difference_pp,
      mu_joint_A, mu_joint_B, mu_A_to_B, mu_B_to_A
    )
  pop_cross <- cross %>%
    transmute(
      kind = "crossing", eps, substitution = substitution_crossing,
      joint = joint_crossing,
      joint_calibrated = calibrated_joint_crossing,
      ratio, calibrated_minus_substitution, calibration_factor,
      delta = NA_real_, s_A = NA_real_, s_B = NA_real_,
      s_joint = NA_real_, R_joint = NA_real_, R_substitution = NA_real_,
      R_joint_A = NA_real_, R_joint_B = NA_real_,
      R_B_to_A = NA_real_, R_A_to_B = NA_real_,
      joint_max_probability_difference_pp = NA_real_,
      substitution_max_probability_difference_pp = NA_real_,
      mu_joint_A = NA_real_, mu_joint_B = NA_real_,
      mu_A_to_B = NA_real_, mu_B_to_A = NA_real_
    )
  bind_rows(pop_cross, pop_curve)
}

# Fitting ----

# Choose among the restarts that follow a failed fit: the converged restart
# with the lowest NLL or, if none converged, the lowest finite NLL overall, in
# which case the run stops at the convergence check.
select_retry_fit <- function(initial, retries) {
  converged <- vapply(
    retries,
    function(fit) fit$convergence == 0 && is.finite(fit$value),
    logical(1)
  )
  candidates <- if (any(converged)) {
    retries[converged]
  } else {
    Filter(function(fit) is.finite(fit$value), c(list(initial), retries))
  }
  if (length(candidates) == 0) {
    stop("Every optimizer start returned a non-finite objective.")
  }
  candidates[[which.min(vapply(candidates, `[[`, numeric(1), "value"))]]
}

# Separate fit of slope and location to one task, in (log slope, location).
# The restart grid is used only when the first fit does not converge.
fit_sep <- function(k, m, x = X_LEVELS) {
  f <- function(par) nll_task(exp(par[1]), par[2], k, m, x)
  o <- optim(c(log(2), 0), f, method = "L-BFGS-B",
             lower = c(log(S_LO), MU_LO), upper = c(log(S_HI), MU_HI))
  reruns <- 0L
  if (o$convergence != 0) {
    retries <- list()
    for (ls in LOGS_STARTS) for (mus in MU_STARTS) {
      reruns <- reruns + 1L
      retries[[reruns]] <- optim(
        c(ls, mus), f, method = "L-BFGS-B",
        lower = c(log(S_LO), MU_LO), upper = c(log(S_HI), MU_HI)
      )
    }
    o <- select_retry_fit(o, retries)
  }
  s <- exp(o$par[1])
  mu <- o$par[2]
  s_bound <- bound_diagnostic(
    o$par[1], LOG_S_LO, LOG_S_HI, s, S_LO, S_HI
  )
  mu_bound <- bound_diagnostic(mu, MU_LO, MU_HI)
  bnd <- s_bound$flag || mu_bound$flag
  list(s = s, mu = mu, nll = o$value,
       conv = (o$convergence == 0), bnd = bnd, reruns = reruns,
       bound_s = s_bound, bound_mu = mu_bound)
}

# Joint fit: one slope, a location for each task. Its summed NLL is the
# likelihood used for AIC, BIC and the likelihood-ratio test.
fit_joint <- function(kA, kB, m, x_A = X_LEVELS, x_B = X_LEVELS) {
  f <- function(par) nll_task(exp(par[1]), par[2], kA, m, x_A) +
                     nll_task(exp(par[1]), par[3], kB, m, x_B)
  o <- optim(c(log(2), 0, 0), f, method = "L-BFGS-B",
             lower = c(log(S_LO), MU_LO, MU_LO),
             upper = c(log(S_HI), MU_HI, MU_HI))
  reruns <- 0L
  if (o$convergence != 0) {
    retries <- list()
    for (ls in LOGS_STARTS) {
      reruns <- reruns + 1L
      retries[[reruns]] <- optim(
        c(ls, 0, 0), f, method = "L-BFGS-B",
        lower = c(log(S_LO), MU_LO, MU_LO),
        upper = c(log(S_HI), MU_HI, MU_HI)
      )
    }
    o <- select_retry_fit(o, retries)
  }
  s <- exp(o$par[1])
  s_bound <- bound_diagnostic(
    o$par[1], LOG_S_LO, LOG_S_HI, s, S_LO, S_HI
  )
  muA_bound <- bound_diagnostic(o$par[2], MU_LO, MU_HI)
  muB_bound <- bound_diagnostic(o$par[3], MU_LO, MU_HI)
  bnd <- s_bound$flag || muA_bound$flag || muB_bound$flag
  list(s = s, muA = o$par[2], muB = o$par[3], nll = o$value,
       conv = (o$convergence == 0), bnd = bnd, reruns = reruns,
       bound_s = s_bound, bound_muA = muA_bound, bound_muB = muB_bound)
}

# Substitution: hold the source slope fixed and fit the target location. The
# default tolerance is optimize()'s default; Simulation 2B passes 1e-10.
fit_target_location <- function(
    s_fixed, k, m, x = X_LEVELS,
    tolerance = .Machine$double.eps^0.25) {
  o <- optimize(
    function(mu) nll_task(s_fixed, mu, k, m, x),
    interval = c(MU_LO, MU_HI), tol = tolerance
  )
  mu_bound <- bound_diagnostic(o$minimum, MU_LO, MU_HI)
  list(mu = o$minimum, nll = o$objective,
       bnd = mu_bound$flag, bound_mu = mu_bound)
}

# Simulation ----

# One simulated experiment at slope gap delta: generate both tasks, fit the
# separate, joint and substitution models, and score every fitted curve by its
# population excess loss.
run_rep <- function(delta, m = M_TRAIN) {
  sA <- 2 - delta / 2
  sB <- 2 + delta / 2
  p0A <- pfun(sA, MU_A)
  p0B <- pfun(sB, MU_B)
  kA <- rbinom(length(X_LEVELS), m, p0A)
  kB <- rbinom(length(X_LEVELS), m, p0B)

  fA <- fit_sep(kA, m)
  fB <- fit_sep(kB, m)
  fJ <- fit_joint(kA, kB, m)
  target_B_fit <- fit_target_location(fA$s, kB, m) # A -> B
  target_A_fit <- fit_target_location(fB$s, kA, m) # B -> A

  RJ_A <- R_T(fJ$s, fJ$muA, p0A)
  RJ_B <- R_T(fJ$s, fJ$muB, p0B)
  RS_B <- R_T(fA$s, target_B_fit$mu, p0B)
  RS_A <- R_T(fB$s, target_A_fit$mu, p0A)
  RJ <- max(RJ_A, RJ_B)
  RS <- max(RS_B, RS_A)

  # Errors on the response-probability scale.
  joint_probability_rmse <- sqrt(mean(c(
    (pfun(fJ$s, fJ$muA) - p0A)^2,
    (pfun(fJ$s, fJ$muB) - p0B)^2
  )))
  joint_probability_max_error <- max(abs(c(
    pfun(fJ$s, fJ$muA) - p0A,
    pfun(fJ$s, fJ$muB) - p0B
  )))
  substitution_probability_rmse <- sqrt(mean(c(
    (pfun(fA$s, target_B_fit$mu) - p0B)^2,
    (pfun(fB$s, target_A_fit$mu) - p0A)^2
  )))
  substitution_probability_max_error <- max(abs(c(
    pfun(fA$s, target_B_fit$mu) - p0B,
    pfun(fB$s, target_A_fit$mu) - p0A
  )))

  nll_shared <- fJ$nll
  nll_sep <- fA$nll + fB$nll
  lrt_raw <- 2 * (nll_shared - nll_sep)
  if (lrt_raw < -LRT_NESTING_TOL) {
    stop(sprintf(
      "Separate fit is worse than the nested shared fit (LRT = %.9g).",
      lrt_raw
    ))
  }
  lrt <- max(0, lrt_raw)

  parameter_bounds <- list(
    sep_A_s = fA$bound_s,
    sep_A_mu = fA$bound_mu,
    sep_B_s = fB$bound_s,
    sep_B_mu = fB$bound_mu,
    joint_s = fJ$bound_s,
    joint_mu_A = fJ$bound_muA,
    joint_mu_B = fJ$bound_muB,
    substitution_A_to_B_mu = target_B_fit$bound_mu,
    substitution_B_to_A_mu = target_A_fit$bound_mu
  )
  bound_flags <- vapply(parameter_bounds, `[[`, logical(1), "flag")
  bound_scaled_distances <- vapply(
    parameter_bounds, `[[`, numeric(1), "scaled_distance"
  )
  bound_optimizer_distances <- vapply(
    parameter_bounds, `[[`, numeric(1), "optimizer_distance"
  )
  bound_native_distances <- vapply(
    parameter_bounds, `[[`, numeric(1), "native_distance"
  )
  names(bound_flags) <- BOUND_FLAG_COLUMNS
  names(bound_scaled_distances) <- BOUND_SCALED_DISTANCE_COLUMNS
  names(bound_optimizer_distances) <- BOUND_OPTIMIZER_DISTANCE_COLUMNS
  names(bound_native_distances) <- BOUND_NATIVE_DISTANCE_COLUMNS

  fit_bound_flags <- c(
    fA$bnd, fB$bnd, fJ$bnd, target_B_fit$bnd, target_A_fit$bnd
  )

  c(RJ = RJ, RS = RS,
    R_joint_A = RJ_A, R_joint_B = RJ_B,
    R_B_to_A = RS_A, R_A_to_B = RS_B,
    sA_hat = fA$s, sB_hat = fB$s, sJ_hat = fJ$s,
    muA_hat = fA$mu, muB_hat = fB$mu,
    mu_joint_A_hat = fJ$muA, mu_joint_B_hat = fJ$muB,
    mu_A_to_B_hat = target_B_fit$mu, mu_B_to_A_hat = target_A_fit$mu,
    nll_shared = nll_shared, nll_sep = nll_sep, lrt = lrt,
    joint_probability_rmse = joint_probability_rmse,
    joint_probability_max_error = joint_probability_max_error,
    substitution_probability_rmse = substitution_probability_rmse,
    substitution_probability_max_error = substitution_probability_max_error,
    convA = fA$conv, convB = fB$conv, convJ = fJ$conv,
    bndA = fA$bnd, bndB = fB$bnd, bndJ = fJ$bnd,
    bound_A_to_B_fit = target_B_fit$bnd,
    bound_B_to_A_fit = target_A_fit$bnd,
    reruns = fA$reruns + fB$reruns + fJ$reruns,
    sJ_between = as.numeric(fJ$s >= min(fA$s, fB$s) &&
                            fJ$s <= max(fA$s, fB$s)),
    bound_flags,
    bound_scaled_distances,
    bound_optimizer_distances,
    bound_native_distances,
    bound_parameter_count = sum(bound_flags),
    bound_fit_count = sum(fit_bound_flags),
    bound_replication = any(bound_flags))
}

# nrep experiments at one slope gap, one row each.
reps_matrix <- function(delta, nrep, m = M_TRAIN) {
  replicate_rows(nrep, function() run_rep(delta, m))
}

# No replication is discarded, so every fit must converge (after restarts) and
# stay away from the optimizer bounds. Otherwise the run stops.
check_fit_matrix <- function(replication_matrix, label = "simulation block") {
  if (any(replication_matrix[, c("convA", "convB", "convJ")] == 0)) {
    stop(label, ": a fit did not converge after the restarts.")
  }
  if (any(replication_matrix[, "bound_replication"] > 0)) {
    stop(label, ": a parameter estimate is at an optimizer bound.")
  }
  invisible(TRUE)
}

# Bound checks for a block of replications: overall counts and rates, and the
# count, rate and smallest distance for each of the nine parameters. The five
# fitting operations are the two separate fits, the joint fit, and the two
# target-location fits.
boundary_summary <- function(replication_matrix) {
  flags <- replication_matrix[, BOUND_FLAG_COLUMNS, drop = FALSE] > 0
  scaled_distances <- replication_matrix[, BOUND_SCALED_DISTANCE_COLUMNS, drop = FALSE]
  fit_flags <- replication_matrix[, c(
    "bndA", "bndB", "bndJ", "bound_A_to_B_fit", "bound_B_to_A_fit"
  ), drop = FALSE] > 0

  values <- list(
    boundary_rate = mean(fit_flags),
    boundary_fit_count = sum(fit_flags),
    boundary_parameter_rate = mean(flags),
    boundary_parameter_count = sum(flags),
    boundary_replication_rate = mean(rowSums(flags) > 0),
    boundary_replication_count = sum(rowSums(flags) > 0),
    min_scaled_bound_distance = min(scaled_distances)
  )

  for (i in seq_along(BOUND_PARAMETER_IDS)) {
    id <- BOUND_PARAMETER_IDS[i]
    values[[paste0("bound_count_", id)]] <- sum(flags[, i])
    values[[paste0("bound_rate_", id)]] <- mean(flags[, i])
    values[[paste0("min_scaled_bound_distance_", id)]] <-
      min(replication_matrix[, BOUND_SCALED_DISTANCE_COLUMNS[i]])
    values[[paste0("min_optimizer_bound_distance_", id)]] <-
      min(replication_matrix[, BOUND_OPTIMIZER_DISTANCE_COLUMNS[i]])
    values[[paste0("min_native_bound_distance_", id)]] <-
      min(replication_matrix[, BOUND_NATIVE_DISTANCE_COLUMNS[i]])
  }

  tibble::as_tibble_row(values)
}
