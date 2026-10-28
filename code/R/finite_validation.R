# Finite-validation analysis (Figure S2; Supplement S6.1).
#
# New Simulation 2A training experiments are fitted, and each fitted curve is
# scored both by its exact excess loss and on independent validation trials.
# The two pass/fail decisions are then compared for the same fits. The
# validation score needs the generating probabilities, so it can be computed
# only in a simulation.
#
#   Rscript --vanilla code/R/finite_validation.R               Figure S2 from saved results
#   Rscript --vanilla code/R/finite_validation.R --resimulate  rerun the simulation (a few minutes)

source("code/R/common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})
source("code/R/psychometric_model.R")

# Settings ----

FV_BASE_SEED <- 20260821L
FV_NREP <- 10000L
FV_DELTAS <- c(0, 0.40, 0.50, 0.60, 0.90)
FV_VALIDATION_M <- c(20L, 80L, 320L) # trials per level
FV_VALIDATION_N <- FV_VALIDATION_M * length(X_LEVELS)
FV_TRAIN_N <- M_TRAIN * length(X_LEVELS)

# Training and validation draws have their own seed for each gap, so changing
# the validation samples leaves the training fits unchanged.
FV_TRAIN_SEEDS <- setNames(
  FV_BASE_SEED + c(101L, 202L, 303L, 404L, 505L),
  sprintf("%.2f", FV_DELTAS)
)
FV_VALIDATION_SEEDS <- setNames(
  FV_BASE_SEED + 10000L + c(101L, 202L, 303L, 404L, 505L),
  sprintf("%.2f", FV_DELTAS)
)

FV_RESULTS_DIR <- file.path(results_dir, "finite_validation")
FV_REPLICATION_PATH <- file.path(
  FV_RESULTS_DIR, "finite_validation_replications.csv.gz"
)

delta_label <- function(delta) sprintf("%.2f", delta)

# Where each gap lies relative to the Simulation 2A population boundaries.
gap_role <- function(delta) {
  roles <- c(
    "0.00" = "true slope equality",
    "0.40" = "below the substitution population boundary",
    "0.50" = "focused condition just beyond the substitution boundary",
    "0.60" = "interior of the raw common-tolerance conflation region",
    "0.90" = "near the raw-joint population boundary"
  )
  unname(roles[delta_label(delta)])
}

# One row per gap and validation size: design, seeds, tolerances and the
# population boundaries.
specification_table <- function(
    nrep = FV_NREP,
    deltas = FV_DELTAS,
    validation_m = FV_VALIDATION_M) {
  crossing <- run_population()$cross %>%
    filter(abs(eps - EPS_PRIMARY) < 1e-12)

  crossing_substitution <- crossing$substitution_crossing[[1]]
  crossing_joint <- crossing$joint_crossing[[1]]
  crossing_calibrated <- crossing$calibrated_joint_crossing[[1]]

  tidyr::crossing(
    delta = deltas,
    validation_trials_per_level = validation_m
  ) %>%
    mutate(
      condition_role = vapply(delta, gap_role, character(1)),
      replications = as.integer(nrep),
      training_trials_per_level = M_TRAIN,
      training_trials_per_task = FV_TRAIN_N,
      validation_trials_per_task =
        validation_trials_per_level * length(X_LEVELS),
      validation_samples_nested = TRUE,
      train_seed = unname(FV_TRAIN_SEEDS[delta_label(delta)]),
      validation_seed = unname(FV_VALIDATION_SEEDS[delta_label(delta)]),
      stimulus_levels = paste(X_LEVELS, collapse = ";"),
      lapse = LAMBDA,
      mu_A = MU_A,
      mu_B = MU_B,
      slope_A = 2 - delta / 2,
      slope_B = 2 + delta / 2,
      epsilon_raw_joint = EPS_PRIMARY,
      epsilon_substitution = EPS_PRIMARY,
      joint_calibration_factor = JOINT_CALIBRATION_FACTOR,
      epsilon_calibrated_joint = EPS_JOINT_CALIBRATED,
      substitution_population_boundary = crossing_substitution,
      raw_joint_population_boundary = crossing_joint,
      calibrated_joint_population_boundary = crossing_calibrated,
      score = "oracle-referenced paired log-score regret",
      comparison = paste(
        "finite-validation decision versus exact-risk decision for the",
        "same fitted parameter estimates"
      )
    ) %>%
    arrange(delta, validation_trials_per_level)
}

# Validation score ----

# Mean log-score difference between the generating probabilities p0 and the
# fitted curve q, from k successes in m validation trials at each level. Given
# the fit, its expectation is the exact mean KL excess loss. It can be
# negative and is not truncated.
log_score_difference <- function(k, m, p0, q) {
  p0 <- clip_probability(p0)
  q <- clip_probability(q)
  sum(
    k * (log(p0) - log(q)) +
      (m - k) * (log1p(-p0) - log1p(-q))
  ) / (length(X_LEVELS) * m)
}

# Variance of log_score_difference() for a fixed fit.
log_score_variance <- function(m, p0, q) {
  p0 <- clip_probability(p0)
  q <- clip_probability(q)
  log_odds_ratio <- log(p0 * (1 - q) / (q * (1 - p0)))
  sum(p0 * (1 - p0) * log_odds_ratio^2) /
    (length(X_LEVELS)^2 * m)
}

# Pass/fail decisions from the exact excess losses and from the validation
# scores, for the same fits.
add_decisions <- function(x) {
  x %>%
    mutate(
      exact_joint_raw = pmax(
        R_joint_A_exact, R_joint_B_exact
      ) <= EPS_PRIMARY,
      validation_joint_raw = pmax(
        R_joint_A_validation, R_joint_B_validation
      ) <= EPS_PRIMARY,
      exact_joint_calibrated = pmax(
        R_joint_A_exact, R_joint_B_exact
      ) <= EPS_JOINT_CALIBRATED,
      validation_joint_calibrated = pmax(
        R_joint_A_validation, R_joint_B_validation
      ) <= EPS_JOINT_CALIBRATED,
      exact_B_to_A = R_B_to_A_exact <= EPS_PRIMARY,
      validation_B_to_A = R_B_to_A_validation <= EPS_PRIMARY,
      exact_A_to_B = R_A_to_B_exact <= EPS_PRIMARY,
      validation_A_to_B = R_A_to_B_validation <= EPS_PRIMARY,
      exact_substitution = exact_B_to_A & exact_A_to_B,
      validation_substitution =
        validation_B_to_A & validation_A_to_B,
      exact_conflation = exact_joint_raw & !exact_substitution,
      validation_conflation =
        validation_joint_raw & !validation_substitution
    )
}

# Simulation ----

# One slope gap: nrep training experiments, each scored on three nested
# validation samples. The smaller samples are the first trials of the larger
# ones, so validation_m must be increasing.
simulate_condition <- function(
    delta,
    nrep,
    validation_m,
    train_seed,
    validation_seed) {
  message(
    sprintf(
      "[finite validation] delta=%.2f: %s training experiments",
      delta, format(nrep, big.mark = ",")
    )
  )

  set.seed(train_seed)
  replication_matrix <- reps_matrix(delta, nrep, m = M_TRAIN)
  check_fit_matrix(
    replication_matrix,
    sprintf("finite-validation training at delta %.2f", delta)
  )

  boundary <- bind_cols(
    tibble(delta = delta, replications = nrep),
    boundary_summary(replication_matrix)
  )
  convergence <- tibble(
    delta = delta,
    replications = nrep,
    component_fits = 3L * nrep,
    converged_component_fits = sum(
      replication_matrix[, c("convA", "convB", "convJ"), drop = FALSE]
    ),
    nonconverged_component_fits =
      component_fits - converged_component_fits,
    convergence_rate = converged_component_fits / component_fits,
    replications_invoking_retry = sum(replication_matrix[, "reruns"] > 0),
    retry_rate = replications_invoking_retry / replications,
    retained_replications = nrep,
    discarded_replications = 0L
  )

  n_sizes <- length(validation_m)
  n_rows <- nrep * n_sizes
  increments <- diff(c(0L, validation_m))

  validation_counts_A <- matrix(
    0L, nrow = n_rows, ncol = length(X_LEVELS)
  )
  validation_counts_B <- matrix(
    0L, nrow = n_rows, ncol = length(X_LEVELS)
  )
  colnames(validation_counts_A) <- sprintf("k_A_level_%02d", seq_along(X_LEVELS))
  colnames(validation_counts_B) <- sprintf("k_B_level_%02d", seq_along(X_LEVELS))

  score_names <- c(
    "R_joint_A_validation", "R_joint_B_validation",
    "R_B_to_A_validation", "R_A_to_B_validation"
  )
  variance_names <- sub("_validation$", "_conditional_variance", score_names)
  validation_scores <- matrix(NA_real_, nrow = n_rows, ncol = 4L,
                              dimnames = list(NULL, score_names))
  conditional_variances <- matrix(
    NA_real_, nrow = n_rows, ncol = 4L,
    dimnames = list(NULL, variance_names)
  )

  sA_true <- 2 - delta / 2
  sB_true <- 2 + delta / 2
  p0A <- pfun(sA_true, MU_A)
  p0B <- pfun(sB_true, MU_B)

  max_exact_recalculation_error <- 0
  set.seed(validation_seed)

  for (i in seq_len(nrep)) {
    q_joint_A <- pfun(
      replication_matrix[i, "sJ_hat"],
      replication_matrix[i, "mu_joint_A_hat"]
    )
    q_joint_B <- pfun(
      replication_matrix[i, "sJ_hat"],
      replication_matrix[i, "mu_joint_B_hat"]
    )
    q_B_to_A <- pfun(replication_matrix[i, "sB_hat"], replication_matrix[i, "mu_B_to_A_hat"])
    q_A_to_B <- pfun(replication_matrix[i, "sA_hat"], replication_matrix[i, "mu_A_to_B_hat"])

    recalculated <- c(
      R_T(
        replication_matrix[i, "sJ_hat"],
        replication_matrix[i, "mu_joint_A_hat"], p0A
      ),
      R_T(
        replication_matrix[i, "sJ_hat"],
        replication_matrix[i, "mu_joint_B_hat"], p0B
      ),
      R_T(replication_matrix[i, "sB_hat"], replication_matrix[i, "mu_B_to_A_hat"], p0A),
      R_T(replication_matrix[i, "sA_hat"], replication_matrix[i, "mu_A_to_B_hat"], p0B)
    )
    retained <- replication_matrix[i, c(
      "R_joint_A", "R_joint_B", "R_B_to_A", "R_A_to_B"
    )]
    max_exact_recalculation_error <- max(
      max_exact_recalculation_error,
      max(abs(recalculated - retained))
    )

    cumulative_A <- integer(length(X_LEVELS))
    cumulative_B <- integer(length(X_LEVELS))
    for (v in seq_len(n_sizes)) {
      cumulative_A <- cumulative_A +
        rbinom(length(X_LEVELS), increments[v], p0A)
      cumulative_B <- cumulative_B +
        rbinom(length(X_LEVELS), increments[v], p0B)
      row <- (i - 1L) * n_sizes + v
      m <- validation_m[v]

      validation_counts_A[row, ] <- cumulative_A
      validation_counts_B[row, ] <- cumulative_B
      validation_scores[row, ] <- c(
        log_score_difference(cumulative_A, m, p0A, q_joint_A),
        log_score_difference(cumulative_B, m, p0B, q_joint_B),
        log_score_difference(cumulative_A, m, p0A, q_B_to_A),
        log_score_difference(cumulative_B, m, p0B, q_A_to_B)
      )
      conditional_variances[row, ] <- c(
        log_score_variance(m, p0A, q_joint_A),
        log_score_variance(m, p0B, q_joint_B),
        log_score_variance(m, p0A, q_B_to_A),
        log_score_variance(m, p0B, q_A_to_B)
      )
    }
  }

  fit_row <- rep(seq_len(nrep), each = n_sizes)
  retained_fit_columns <- c(
    "R_joint_A", "R_joint_B", "R_B_to_A", "R_A_to_B",
    "sA_hat", "sB_hat", "sJ_hat", "mu_joint_A_hat", "mu_joint_B_hat",
    "mu_B_to_A_hat", "mu_A_to_B_hat", "convA", "convB", "convJ",
    "reruns", "bound_replication", "bound_parameter_count"
  )
  fits <- as_tibble(replication_matrix[fit_row, retained_fit_columns, drop = FALSE]) %>%
    rename(
      R_joint_A_exact = R_joint_A,
      R_joint_B_exact = R_joint_B,
      R_B_to_A_exact = R_B_to_A,
      R_A_to_B_exact = R_A_to_B
    ) %>%
    mutate(
      across(c(convA, convB, convJ, bound_replication), ~ .x > 0),
      reruns = as.integer(reruns),
      bound_parameter_count = as.integer(bound_parameter_count)
    )

  replications <- bind_cols(
    tibble(
      delta = delta,
      replication = fit_row,
      training_trials_per_level = M_TRAIN,
      training_trials_per_task = FV_TRAIN_N,
      validation_trials_per_level = rep(validation_m, times = nrep),
      validation_trials_per_task =
        rep(validation_m, times = nrep) * length(X_LEVELS),
      train_seed = train_seed,
      validation_seed = validation_seed
    ),
    fits,
    as_tibble(validation_scores),
    as_tibble(conditional_variances),
    as_tibble(validation_counts_A),
    as_tibble(validation_counts_B)
  ) %>%
    add_decisions()

  nested_counts_ok <- TRUE
  for (i in seq_len(nrep)) {
    rows <- (i - 1L) * n_sizes + seq_len(n_sizes)
    nested_counts_ok <- nested_counts_ok &&
      all(apply(validation_counts_A[rows, , drop = FALSE], 2, function(z) {
        observed_increments <- diff(c(0, z))
        all(observed_increments >= 0 & observed_increments <= increments)
      })) &&
      all(apply(validation_counts_B[rows, , drop = FALSE], 2, function(z) {
        observed_increments <- diff(c(0, z))
        all(observed_increments >= 0 & observed_increments <= increments)
      }))
  }

  condition_checks <- tibble(
    delta = delta,
    check = c(
      "retained exact risks reproduce recomputed exact risks",
      "validation counts form valid nested binomial increments",
      "all validation scores are finite"
    ),
    value = c(
      max_exact_recalculation_error,
      as.numeric(nested_counts_ok),
      as.numeric(all(is.finite(validation_scores)))
    ),
    target = c(0, 1, 1),
    tolerance = c(1e-12, 0, 0),
    pass = c(
      max_exact_recalculation_error <= 1e-12,
      nested_counts_ok,
      all(is.finite(validation_scores))
    )
  )

  list(
    replications = replications,
    convergence = convergence,
    boundary = boundary,
    condition_checks = condition_checks
  )
}

# All gaps, each with its own training and validation seed.
simulate_validation <- function(
    nrep = FV_NREP,
    deltas = FV_DELTAS,
    validation_m = FV_VALIDATION_M) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  results <- vector("list", length(deltas))
  for (i in seq_along(deltas)) {
    key <- delta_label(deltas[i])
    results[[i]] <- simulate_condition(
      delta = deltas[i],
      nrep = nrep,
      validation_m = validation_m,
      train_seed = unname(FV_TRAIN_SEEDS[key]),
      validation_seed = unname(FV_VALIDATION_SEEDS[key])
    )
  }

  list(
    replications = bind_rows(map(results, "replications")),
    convergence = bind_rows(map(results, "convergence")),
    boundary = bind_rows(map(results, "boundary")),
    condition_checks = bind_rows(map(results, "condition_checks"))
  )
}

# Summaries and checks ----

# Add Wilson limits <prefix>_lo and <prefix>_hi for count / denominator, NA
# when the denominator is zero.
add_wilson_interval <- function(data, count, denominator, prefix) {
  n <- data[[denominator]]
  interval <- wilson_interval(data[[count]], n)
  interval[n == 0, ] <- NA
  data[[paste0(prefix, "_lo")]] <- interval[, "lower"]
  data[[paste0(prefix, "_hi")]] <- interval[, "upper"]
  data
}

# One row per experiment, validation size and criterion.
decisions_long <- function(replications) {
  bind_rows(
    replications %>% transmute(
      delta, replication, validation_trials_per_level,
      validation_trials_per_task, method = "raw_joint",
      exact_pass = exact_joint_raw,
      validation_pass = validation_joint_raw
    ),
    replications %>% transmute(
      delta, replication, validation_trials_per_level,
      validation_trials_per_task, method = "calibrated_joint",
      exact_pass = exact_joint_calibrated,
      validation_pass = validation_joint_calibrated
    ),
    replications %>% transmute(
      delta, replication, validation_trials_per_level,
      validation_trials_per_task, method = "bidirectional_substitution",
      exact_pass = exact_substitution,
      validation_pass = validation_substitution
    )
  )
}

# Agreement between validation and exact decisions, with the decisions in long
# form and `by` the column that identifies the criterion or direction.
agreement_summary <- function(long, by) {
  long %>%
    group_by(across(all_of(c(
      "delta", "validation_trials_per_level", "validation_trials_per_task", by
    )))) %>%
    summarise(
      replications = n(),
      exact_pass_count = sum(exact_pass),
      validation_pass_count = sum(validation_pass),
      exact_fail_count = sum(!exact_pass),
      validation_fail_count = sum(!validation_pass),
      validation_pass_exact_fail_count =
        sum(validation_pass & !exact_pass),
      validation_fail_exact_pass_count =
        sum(!validation_pass & exact_pass),
      disagreement_count = sum(validation_pass != exact_pass),
      .groups = "drop"
    ) %>%
    mutate(
      exact_pass_rate = exact_pass_count / replications,
      validation_pass_rate = validation_pass_count / replications,
      validation_pass_exact_fail_rate =
        validation_pass_exact_fail_count / replications,
      validation_fail_exact_pass_rate =
        validation_fail_exact_pass_count / replications,
      disagreement_rate = disagreement_count / replications,
      conditional_validation_pass_given_exact_fail = if_else(
        exact_fail_count > 0,
        validation_pass_exact_fail_count / exact_fail_count,
        NA_real_
      ),
      conditional_validation_fail_given_exact_pass = if_else(
        exact_pass_count > 0,
        validation_fail_exact_pass_count / exact_pass_count,
        NA_real_
      )
    )
}

# For each criterion.
summarize_decisions <- function(replications) {
  out <- agreement_summary(decisions_long(replications), "method")

  interval_specs <- tribble(
    ~count, ~denominator, ~prefix,
    "exact_pass_count", "replications", "exact_pass",
    "validation_pass_count", "replications", "validation_pass",
    "validation_pass_exact_fail_count", "replications",
      "validation_pass_exact_fail",
    "validation_fail_exact_pass_count", "replications",
      "validation_fail_exact_pass",
    "disagreement_count", "replications", "disagreement",
    "validation_pass_exact_fail_count", "exact_fail_count",
      "conditional_validation_pass_given_exact_fail",
    "validation_fail_exact_pass_count", "exact_pass_count",
      "conditional_validation_fail_given_exact_pass"
  )
  for (i in seq_len(nrow(interval_specs))) {
    out <- add_wilson_interval(
      out,
      interval_specs$count[i],
      interval_specs$denominator[i],
      interval_specs$prefix[i]
    )
  }
  out %>% arrange(method, delta, validation_trials_per_task)
}

# For each transfer direction.
summarize_directions <- function(replications) {
  long <- bind_rows(
    replications %>% transmute(
      delta, replication, validation_trials_per_level,
      validation_trials_per_task, direction = "B_to_A",
      exact_pass = exact_B_to_A,
      validation_pass = validation_B_to_A
    ),
    replications %>% transmute(
      delta, replication, validation_trials_per_level,
      validation_trials_per_task, direction = "A_to_B",
      exact_pass = exact_A_to_B,
      validation_pass = validation_A_to_B
    )
  )

  out <- agreement_summary(long, "direction") %>%
    select(-validation_fail_count)
  out <- add_wilson_interval(out, "disagreement_count", "replications",
                         "disagreement")
  out %>% arrange(direction, delta, validation_trials_per_task)
}

# Joint adequacy without transfer ("conflation"), directional disagreements,
# and reversals of the transfer direction.
summarize_events <- function(replications) {
  out <- replications %>%
    group_by(delta, validation_trials_per_level, validation_trials_per_task) %>%
    summarise(
      replications = n(),
      exact_conflation_count = sum(exact_conflation),
      validation_conflation_count = sum(validation_conflation),
      validation_conflation_exact_not_count =
        sum(validation_conflation & !exact_conflation),
      exact_conflation_validation_not_count =
        sum(!validation_conflation & exact_conflation),
      conflation_disagreement_count =
        sum(validation_conflation != exact_conflation),
      exact_directional_disagreement_count =
        sum(exact_B_to_A != exact_A_to_B),
      validation_directional_disagreement_count =
        sum(validation_B_to_A != validation_A_to_B),
      direction_reversal_count = sum(
        (exact_B_to_A & !exact_A_to_B &
           !validation_B_to_A & validation_A_to_B) |
          (!exact_B_to_A & exact_A_to_B &
             validation_B_to_A & !validation_A_to_B)
      ),
      exact_calibrated_joint_substitution_mismatch_count =
        sum(exact_joint_calibrated != exact_substitution),
      validation_calibrated_joint_substitution_mismatch_count =
        sum(validation_joint_calibrated != validation_substitution),
      .groups = "drop"
    )

  count_columns <- names(out)[str_detect(names(out), "_count$")]
  for (column in count_columns) {
    rate_name <- sub("_count$", "_rate", column)
    out[[rate_name]] <- out[[column]] / out$replications
  }
  out <- add_wilson_interval(
    out, "conflation_disagreement_count", "replications",
    "conflation_disagreement"
  )
  out %>% arrange(delta, validation_trials_per_task)
}

# Exact-by-validation 2 x 2 tables for each criterion.
tabulate_contingencies <- function(replications) {
  decisions <- decisions_long(replications) %>%
    transmute(
      delta, replication, validation_trials_per_level,
      validation_trials_per_task, event = method,
      exact = exact_pass, validation = validation_pass
    )
  conflation <- replications %>%
    transmute(
      delta, replication, validation_trials_per_level,
      validation_trials_per_task, event = "predictive_conflation",
      exact = exact_conflation, validation = validation_conflation
    )

  bind_rows(decisions, conflation) %>%
    count(
      delta, validation_trials_per_level, validation_trials_per_task,
      event, exact, validation, name = "count"
    ) %>%
    group_by(
      delta, validation_trials_per_level,
      validation_trials_per_task, event
    ) %>%
    tidyr::complete(
      exact = c(FALSE, TRUE),
      validation = c(FALSE, TRUE),
      fill = list(count = 0L)
    ) %>%
    mutate(denominator = sum(count)) %>%
    ungroup() %>%
    arrange(event, delta, validation_trials_per_task, exact, validation)
}

# Validation score minus exact excess loss: bias, RMSE, and the empirical
# mean squared error against the analytic variance.
summarize_scores <- function(replications) {
  score_map <- tribble(
    ~score, ~exact_column, ~validation_column, ~variance_column,
    "joint_A", "R_joint_A_exact", "R_joint_A_validation",
      "R_joint_A_conditional_variance",
    "joint_B", "R_joint_B_exact", "R_joint_B_validation",
      "R_joint_B_conditional_variance",
    "B_to_A", "R_B_to_A_exact", "R_B_to_A_validation",
      "R_B_to_A_conditional_variance",
    "A_to_B", "R_A_to_B_exact", "R_A_to_B_validation",
      "R_A_to_B_conditional_variance"
  )

  long <- map_dfr(seq_len(nrow(score_map)), function(i) {
    exact <- replications[[score_map$exact_column[i]]]
    validation <- replications[[score_map$validation_column[i]]]
    variance <- replications[[score_map$variance_column[i]]]
    tibble(
      delta = replications$delta,
      replication = replications$replication,
      validation_trials_per_level =
        replications$validation_trials_per_level,
      validation_trials_per_task =
        replications$validation_trials_per_task,
      score = score_map$score[i],
      exact_risk = exact,
      validation_score = validation,
      score_error = validation - exact,
      analytic_conditional_variance = variance
    )
  })

  long %>%
    group_by(
      delta, validation_trials_per_level, validation_trials_per_task, score
    ) %>%
    summarise(
      replications = n(),
      mean_exact_risk = mean(exact_risk),
      mean_validation_score = mean(validation_score),
      mean_score_error = mean(score_error),
      score_error_mcse = sd(score_error) / sqrt(replications),
      score_error_z = if_else(
        score_error_mcse > 0,
        mean_score_error / score_error_mcse,
        0
      ),
      score_error_rmse = sqrt(mean(score_error^2)),
      mean_analytic_conditional_variance =
        mean(analytic_conditional_variance),
      empirical_score_error_mse = mean(score_error^2),
      empirical_mse_to_analytic_variance_ratio =
        empirical_score_error_mse / mean_analytic_conditional_variance,
      negative_validation_score_count = sum(validation_score < 0),
      negative_validation_score_rate = mean(validation_score < 0),
      .groups = "drop"
    ) %>%
    arrange(score, delta, validation_trials_per_task)
}

# The new training experiments should reproduce the exact pass rates of the
# Simulation 2A main grid at the same gaps. The run stops if any rate differs
# by more than four Monte Carlo standard errors.
compare_with_main_grid <- function(decision_summary) {
  main_grid_path <- file.path(
    results_dir, "simulation2a", "finite_sample_by_delta.csv"
  )
  retained <- read_csv(main_grid_path, show_col_types = FALSE) %>%
    mutate(delta_key = round(delta, 8)) %>%
    filter(delta_key %in% round(FV_DELTAS, 8)) %>%
    transmute(
      delta_key,
      retained_replications = replications,
      raw_joint = raw_joint_pass_rate,
      calibrated_joint = calibrated_joint_pass_rate,
      bidirectional_substitution = substitution_pass_rate
    ) %>%
    pivot_longer(
      c(raw_joint, calibrated_joint, bidirectional_substitution),
      names_to = "method", values_to = "retained_exact_pass_rate"
    )

  current <- decision_summary %>%
    filter(validation_trials_per_task == min(FV_VALIDATION_N)) %>%
    transmute(
      delta = delta,
      delta_key = round(delta, 8),
      method,
      current_replications = replications,
      current_exact_pass_rate = exact_pass_rate
    )

  comparison <- inner_join(current, retained, by = c("delta_key", "method")) %>%
    mutate(
      rate_difference =
        current_exact_pass_rate - retained_exact_pass_rate,
      independent_mcse = sqrt(
        current_exact_pass_rate * (1 - current_exact_pass_rate) /
          current_replications +
          retained_exact_pass_rate * (1 - retained_exact_pass_rate) /
            retained_replications
      ),
      z_difference = case_when(
        independent_mcse > 0 ~ rate_difference / independent_mcse,
        rate_difference == 0 ~ 0,
        TRUE ~ sign(rate_difference) * Inf
      ),
      compatible_with_retained_main_grid = abs(z_difference) <= 4
    ) %>%
    select(-delta_key) %>%
    arrange(method, delta)

  expected_rows <- length(FV_DELTAS) * 3L
  if (nrow(comparison) != expected_rows ||
      any(!comparison$compatible_with_retained_main_grid)) {
    stop("The training experiments do not reproduce the Simulation 2A ",
         "pass rates.")
  }
  comparison
}

# Standard rebuild: repeat that comparison from the saved summaries.
check_training_consistency <- function(decision_summary) {
  path <- file.path(FV_RESULTS_DIR, "training_consistency.csv")
  expected <- compare_with_main_grid(decision_summary)
  stored <- read_csv(path, show_col_types = FALSE)
  if (!isTRUE(all.equal(
    as.data.frame(stored),
    as.data.frame(expected),
    check.attributes = FALSE,
    tolerance = 1e-12
  ))) {
    stop("training_consistency.csv does not match the saved summaries.")
  }
  invisible(TRUE)
}

# Optimizer-bound checks, one row per gap and parameter.
parameter_checks <- function(boundary_summary) {
  map_dfr(BOUND_PARAMETER_IDS, function(parameter) {
    boundary_summary %>%
      transmute(
        delta,
        replications,
        parameter = parameter,
        count_near_bound = .data[[paste0("bound_count_", parameter)]],
        proportion_near_bound = .data[[paste0("bound_rate_", parameter)]],
        minimum_scaled_distance =
          .data[[paste0("min_scaled_bound_distance_", parameter)]],
        minimum_optimizer_coordinate_distance =
          .data[[paste0("min_optimizer_bound_distance_", parameter)]],
        minimum_native_distance =
          .data[[paste0("min_native_bound_distance_", parameter)]],
        parameters_near_bound_total = boundary_parameter_count,
        fitting_operations_near_bound_total = boundary_fit_count,
        replications_near_bound_total = boundary_replication_count
      )
  }) %>%
    arrange(delta, parameter)
}

# Checks of the score itself on one fixed fit: aggregated counts against
# trial-level scoring, expected counts against the exact KL, and the mean and
# variance of the score over repeated validation samples.
score_checks <- function() {
  p0 <- pfun(1.7, MU_A)
  q <- pfun(2.05, -0.08)
  m <- 20L
  k <- pmin(m, pmax(0L, round(m * p0)))

  aggregate_score <- log_score_difference(k, m, p0, q)
  trial_score <- mean(unlist(map(seq_along(X_LEVELS), function(j) {
    y <- c(rep(1, k[j]), rep(0, m - k[j]))
    y * (log(p0[j]) - log(q[j])) +
      (1 - y) * (log1p(-p0[j]) - log1p(-q[j]))
  })))
  exact_risk <- mean(kl_bern(p0, q))
  fractional_score <- log_score_difference(m * p0, m, p0, q)
  identical_score <- log_score_difference(k, m, p0, p0)

  # 100,000 redraws, with their own seed.
  moment_draws <- 100000L
  moment_m <- 80L
  set.seed(FV_BASE_SEED + 20000L)
  K <- vapply(
    p0,
    function(probability) rbinom(moment_draws, moment_m, probability),
    integer(moment_draws)
  )
  success_term <- log(p0) - log(q)
  failure_term <- log1p(-p0) - log1p(-q)
  score_draws <- rowSums(
    sweep(K, 2, success_term, `*`) +
      sweep(moment_m - K, 2, failure_term, `*`)
  ) / (length(X_LEVELS) * moment_m)
  analytic_variance <- log_score_variance(moment_m, p0, q)
  mean_mcse <- sqrt(analytic_variance / moment_draws)
  variance_relative_error <-
    abs(var(score_draws) - analytic_variance) / analytic_variance

  tibble(
    check = c(
      "aggregated-binomial score equals trial-level score",
      "fractional expected counts reproduce exact KL risk",
      "q equal to p0 gives zero score for every realized count vector",
      "repeated-redraw mean agrees with exact KL risk",
      "repeated-redraw variance agrees with analytic variance"
    ),
    value = c(
      abs(aggregate_score - trial_score),
      abs(fractional_score - exact_risk),
      abs(identical_score),
      abs(mean(score_draws) - exact_risk),
      variance_relative_error
    ),
    target = 0,
    tolerance = c(1e-14, 1e-14, 1e-14, 5 * mean_mcse, 0.03),
    pass = value <= tolerance
  )
}

# Checks on the design and the score, saved as self_checks.csv. The run stops
# if any check fails.
run_self_checks <- function(
    simulation,
    nrep,
    deltas,
    validation_m) {
  x <- simulation$replications
  expected_rows <- nrep * length(deltas) * length(validation_m)
  score_summary <- summarize_scores(x)
  decision_summary <- summarize_decisions(x)

  fit_columns <- c(
    "R_joint_A_exact", "R_joint_B_exact", "R_B_to_A_exact",
    "R_A_to_B_exact", "sA_hat", "sB_hat", "sJ_hat",
    "mu_joint_A_hat", "mu_joint_B_hat",
    "mu_B_to_A_hat", "mu_A_to_B_hat",
    "convA", "convB", "convJ", "reruns", "bound_replication",
    "bound_parameter_count", "exact_joint_raw",
    "exact_joint_calibrated", "exact_substitution"
  )
  fit_invariance <- x %>%
    group_by(delta, replication) %>%
    summarise(
      across(all_of(fit_columns), n_distinct),
      .groups = "drop"
  )
  same_fits <- all(as.matrix(fit_invariance[, fit_columns]) == 1)

  expected_cells <- tidyr::crossing(
    delta = deltas,
    replication = seq_len(nrep),
    validation_trials_per_level = validation_m
  )
  observed_cells <- x %>%
    count(delta, replication, validation_trials_per_level, name = "rows")
  factorial_complete <-
    nrow(observed_cells) == expected_rows &&
    all(observed_cells$rows == 1L) &&
    nrow(anti_join(
      expected_cells, observed_cells,
      by = c("delta", "replication", "validation_trials_per_level")
    )) == 0L &&
    nrow(anti_join(
      observed_cells, expected_cells,
      by = c("delta", "replication", "validation_trials_per_level")
    )) == 0L

  decision_identity <- all(
    x$exact_joint_raw ==
      (pmax(x$R_joint_A_exact, x$R_joint_B_exact) <= EPS_PRIMARY)
  ) && all(
    x$validation_joint_raw ==
      (pmax(
        x$R_joint_A_validation, x$R_joint_B_validation
      ) <= EPS_PRIMARY)
  ) && all(
    x$exact_joint_calibrated ==
      (pmax(x$R_joint_A_exact, x$R_joint_B_exact) <=
         EPS_JOINT_CALIBRATED)
  ) && all(
    x$validation_joint_calibrated ==
      (pmax(
        x$R_joint_A_validation, x$R_joint_B_validation
      ) <= EPS_JOINT_CALIBRATED)
  ) && all(
    x$exact_substitution == (x$exact_B_to_A & x$exact_A_to_B)
  ) && all(
    x$validation_substitution ==
      (x$validation_B_to_A & x$validation_A_to_B)
  ) && all(
    x$exact_conflation == (x$exact_joint_raw & !x$exact_substitution)
  ) && all(
    x$validation_conflation ==
      (x$validation_joint_raw & !x$validation_substitution)
  )
  disagreement_identity <- all(
    decision_summary$disagreement_count ==
      decision_summary$validation_pass_exact_fail_count +
        decision_summary$validation_fail_exact_pass_count
  )

  count_columns <- names(x)[str_detect(names(x), "^k_[AB]_level_[0-9]+$")]
  retained_nesting <- x %>%
    arrange(delta, replication, validation_trials_per_level) %>%
    group_by(delta, replication) %>%
    summarise(
      validation_increment = list(diff(c(0, validation_trials_per_level))),
      across(
        all_of(count_columns),
        ~ list(diff(c(0, .x))),
        .names = "increment_{.col}"
      ),
      .groups = "drop"
    )
  retained_nested_counts_ok <- TRUE
  increment_columns <- names(retained_nesting)[str_starts(
    names(retained_nesting), "increment_k_"
  )]
  for (i in seq_len(nrow(retained_nesting))) {
    allowed <- retained_nesting$validation_increment[[i]]
    for (column in increment_columns) {
      observed <- retained_nesting[[column]][[i]]
      retained_nested_counts_ok <- retained_nested_counts_ok &&
        all(observed >= 0 & observed <= allowed)
    }
  }

  # For a fixed fit, m times the analytic variance is the same at every
  # validation size.
  variance_columns <- names(x)[str_detect(
    names(x), "_conditional_variance$"
  )]
  variance_scaling <- x %>%
    mutate(across(
      all_of(variance_columns),
      ~ .x * validation_trials_per_level
    )) %>%
    group_by(delta, replication) %>%
    summarise(
      across(all_of(variance_columns), ~ max(.x) - min(.x)),
      .groups = "drop"
    )
  max_variance_scaling_error <- max(as.matrix(
    variance_scaling[, variance_columns]
  ))

  # Quick test runs with fewer than 1,000 replications skip the two Monte
  # Carlo checks of the score.
  unbiased_check_required <- nrep >= 1000
  max_abs_score_z <- max(abs(score_summary$score_error_z), na.rm = TRUE)
  variance_ratio_deviation <- max(abs(
    score_summary$empirical_mse_to_analytic_variance_ratio - 1
  ), na.rm = TRUE)
  variance_ratio_tolerance <- if (unbiased_check_required) {
    0.12 * sqrt(FV_NREP / nrep)
  } else {
    Inf
  }

  global <- tibble(
    check = c(
      "retained validation row count matches the declared factorial design",
      "every declared factorial cell appears exactly once",
      "unique training-experiment count matches delta by replication",
      "training experiments and exact decisions are identical across validation sizes",
      "validation trials per task equal nine times trials per level",
      "retained validation counts form valid nested binomial increments",
      "joint, substitution, and conflation decisions reproduce their definitions",
      "validation-only plus exact-only counts equal total disagreement",
      "analytic score variance has exact inverse-m scaling",
      "mean score errors are compatible with conditional unbiasedness",
      "mean squared score error is compatible with mean analytic variance"
    ),
    value = c(
      nrow(x),
      as.numeric(factorial_complete),
      n_distinct(interaction(x$delta, x$replication, drop = TRUE)),
      as.numeric(same_fits),
      as.numeric(all(
        x$validation_trials_per_task ==
          length(X_LEVELS) * x$validation_trials_per_level
      )),
      as.numeric(retained_nested_counts_ok),
      as.numeric(decision_identity),
      as.numeric(disagreement_identity),
      max_variance_scaling_error,
      max_abs_score_z,
      variance_ratio_deviation
    ),
    target = c(
      expected_rows,
      1,
      nrep * length(deltas),
      1, 1, 1, 1, 1, 0, 0, 0
    ),
    tolerance = c(
      0, 0, 0, 0, 0, 0, 0, 0, 1e-14,
      if (unbiased_check_required) 6 else Inf,
      variance_ratio_tolerance
    ),
    pass = c(
      nrow(x) == expected_rows,
      factorial_complete,
      n_distinct(interaction(x$delta, x$replication, drop = TRUE)) ==
        nrep * length(deltas),
      same_fits,
      all(
        x$validation_trials_per_task ==
          length(X_LEVELS) * x$validation_trials_per_level
      ),
      retained_nested_counts_ok,
      decision_identity,
      disagreement_identity,
      max_variance_scaling_error <= 1e-14,
      !unbiased_check_required || max_abs_score_z <= 6,
      !unbiased_check_required ||
        variance_ratio_deviation <= variance_ratio_tolerance
    )
  )

  condition_checks <- simulation$condition_checks %>%
    mutate(check = paste0("delta=", sprintf("%.2f", delta), ": ", check)) %>%
    select(check, value, target, tolerance, pass)
  checks <- bind_rows(
    score_checks(),
    global,
    condition_checks
  )

  if (any(!checks$pass)) {
    failed <- checks %>% filter(!pass)
    stop(
      "Finite-validation self-checks failed: ",
      paste(failed$check, collapse = "; ")
    )
  }
  checks
}

# Save results ----

write_outputs <- function(
    simulation,
    checks,
    nrep = FV_NREP,
    deltas = FV_DELTAS,
    validation_m = FV_VALIDATION_M) {
  dir.create(FV_RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

  outputs <- list(
    specification = specification_table(
      nrep = nrep, deltas = deltas, validation_m = validation_m
    ),
    finite_validation_replications = simulation$replications,
    decision_summary = summarize_decisions(
      simulation$replications
    ),
    direction_summary = summarize_directions(
      simulation$replications
    ),
    event_summary = summarize_events(
      simulation$replications
    ),
    contingency_tables = tabulate_contingencies(
      simulation$replications
    ),
    score_audit = summarize_scores(
      simulation$replications
    ),
    training_consistency = compare_with_main_grid(
      summarize_decisions(simulation$replications)
    ),
    convergence_audit = simulation$convergence,
    parameter_boundary_checks = parameter_checks(
      simulation$boundary
    ),
    self_checks = checks
  )

  paths <- c(
    specification = file.path(FV_RESULTS_DIR, "specification.csv"),
    finite_validation_replications = FV_REPLICATION_PATH,
    decision_summary = file.path(FV_RESULTS_DIR, "decision_summary.csv"),
    direction_summary = file.path(FV_RESULTS_DIR, "direction_summary.csv"),
    event_summary = file.path(FV_RESULTS_DIR, "event_summary.csv"),
    contingency_tables = file.path(FV_RESULTS_DIR, "contingency_tables.csv"),
    score_audit = file.path(FV_RESULTS_DIR, "score_audit.csv"),
    training_consistency = file.path(
      FV_RESULTS_DIR, "training_consistency.csv"
    ),
    convergence_audit = file.path(FV_RESULTS_DIR, "convergence_audit.csv"),
    parameter_boundary_checks = file.path(
      FV_RESULTS_DIR, "parameter_boundary_checks.csv"
    ),
    self_checks = file.path(FV_RESULTS_DIR, "self_checks.csv")
  )

  walk(names(outputs), function(name) write_csv(outputs[[name]], paths[[name]]))
  invisible(outputs)
}

# Figure S2 ----

# Disagreement between validation and exact decisions, by gap and validation
# size, for the joint criterion, bidirectional transfer, and joint adequacy
# without transfer.
figure_s2_data <- function(decision_summary, event_summary) {
  decision_summary %>%
    filter(method %in% c("raw_joint", "bidirectional_substitution")) %>%
    transmute(
      delta, validation_trials_per_task,
      panel = recode(
        method,
        raw_joint = "Joint adequacy",
        bidirectional_substitution = "Bidirectional transfer\nunder substitution"
      ),
      disagreement_rate,
      disagreement_lo,
      disagreement_hi
    ) %>%
    bind_rows(
      event_summary %>%
        transmute(
          delta, validation_trials_per_task,
          panel = "Joint adequacy without transfer",
          disagreement_rate = conflation_disagreement_rate,
          disagreement_lo = conflation_disagreement_lo,
          disagreement_hi = conflation_disagreement_hi
        )
    ) %>%
    mutate(
      panel = factor(
        panel,
        levels = c(
          "Joint adequacy", "Bidirectional transfer\nunder substitution",
          "Joint adequacy without transfer"
        )
      ),
      validation_trials_per_task = factor(
        validation_trials_per_task,
        levels = FV_VALIDATION_N,
        labels = format(FV_VALIDATION_N, big.mark = ",", trim = TRUE)
      )
    )
}

make_figure_s2 <- function(decision_summary, event_summary) {
  plot_data <- figure_s2_data(decision_summary, event_summary)
  validation_colours <- c(blue, orange, green)
  names(validation_colours) <- levels(plot_data$validation_trials_per_task)
  validation_linetypes <- c("solid", "22", "dotdash")
  names(validation_linetypes) <- levels(plot_data$validation_trials_per_task)
  validation_shapes <- c(16, 17, 15)
  names(validation_shapes) <- levels(plot_data$validation_trials_per_task)

  plot <- ggplot(
    plot_data,
    aes(
      x = delta, y = disagreement_rate,
      colour = validation_trials_per_task,
      linetype = validation_trials_per_task,
      shape = validation_trials_per_task,
      group = validation_trials_per_task
    )
  ) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.7) +
    facet_wrap(~panel, nrow = 1) +
    scale_colour_manual(values = validation_colours) +
    scale_linetype_manual(values = validation_linetypes) +
    scale_shape_manual(values = validation_shapes) +
    scale_x_continuous(breaks = FV_DELTAS) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      limits = c(0, NA),
      expand = expansion(mult = c(0.025, 0.08))
    ) +
    labs(
      x = expression("True slope gap " * delta),
      y = "Decision disagreement",
      colour = "Validation responses per task",
      linetype = "Validation responses per task",
      shape = "Validation responses per task"
    ) +
    paper_theme(base_size = 9) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 9, face = "plain"),
      panel.spacing.x = unit(6, "mm"),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9)
    )

  save_figure(
    plot,
    file.path(supplement_figure_dir, "figureS2_finite_validation.pdf"),
    height_mm = 68
  )
}

run_finite_validation <- function(resimulate = FALSE) {
  if (resimulate) {
    simulation <- simulate_validation()
    checks <- run_self_checks(
      simulation,
      nrep = FV_NREP,
      deltas = FV_DELTAS,
      validation_m = FV_VALIDATION_M
    )
    outputs <- write_outputs(simulation, checks)
    decision_summary <- outputs$decision_summary
    event_summary <- outputs$event_summary
  } else {
    decision_summary <- read_csv(
      file.path(FV_RESULTS_DIR, "decision_summary.csv"),
      show_col_types = FALSE
    )
    event_summary <- read_csv(
      file.path(FV_RESULTS_DIR, "event_summary.csv"),
      show_col_types = FALSE
    )
    check_training_consistency(decision_summary)
  }
  make_figure_s2(decision_summary, event_summary)
}

if (sys.nframe() == 0) {
  run_finite_validation(resimulate = "--resimulate" %in% commandArgs(trailingOnly = TRUE))
}
