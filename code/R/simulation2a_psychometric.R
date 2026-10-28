# Simulation 2A (Figure 4; Supplement S3): psychometric slopes with
# task-specific locations.
#
# The focal slope is either shared in a joint fit or carried from one task to
# the other (substitution), while each task keeps its own location. Every
# fitted curve is scored by its population excess loss, computed from the
# known generating probabilities.
#
#   Rscript --vanilla code/R/simulation2a_psychometric.R               Figure 4 from saved results
#   Rscript --vanilla code/R/simulation2a_psychometric.R --resimulate  rerun the simulation (several minutes)
#
# A full run simulates, in this order and from one random-number stream: the
# recovery benchmark, the main grid, the sensitivity grid, and the focused
# 10,000-replication conditions.

source("code/R/common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})
source("code/R/psychometric_model.R")

SEED <- 20240115

# Pass/fail and model-selection decisions for a block of replications. BIC
# uses the total number of trials, 2 * n_per_task.
decision_flags <- function(replication_matrix, eps, n_per_task) {
  lrt <- replication_matrix[, "lrt"]
  list(
    joint_pass = replication_matrix[, "RJ"] <= eps,
    calibrated_joint_pass =
      replication_matrix[, "RJ"] <= eps / JOINT_CALIBRATION_FACTOR,
    substitution_pass = replication_matrix[, "RS"] <= eps,
    aic_shared = lrt < AIC_LRT_CUTOFF,
    bic_shared = lrt < log(2 * n_per_task) * (K_SEPARATE - K_SHARED),
    lrt_shared = lrt < LRT_05_CUTOFF
  )
}

# Main grid ----

# 2,000 replications at each of the 31 slope gaps (Figure 4).
run_finite_main <- function(nrep = 2000) {
  message("Main grid: ", nrep, " replications at ", length(DELTAS), " gaps")
  out <- vector("list", length(DELTAS))
  for (i in seq_along(DELTAS)) {
    d <- DELTAS[i]
    replication_matrix <- reps_matrix(d, nrep)
    check_fit_matrix(replication_matrix, sprintf("main grid at delta %.2f", d))
    decisions <- decision_flags(replication_matrix, EPS_PRIMARY, N_TASK)
    joint_pass <- decisions$joint_pass
    calibrated_joint_pass <- decisions$calibrated_joint_pass
    substitution_pass <- decisions$substitution_pass
    aic_shared <- decisions$aic_shared
    conflation <- joint_pass & !substitution_pass
    calibrated_joint_only <- calibrated_joint_pass & !substitution_pass
    substitution_only <- !calibrated_joint_pass & substitution_pass
    calibrated_disagreement <- calibrated_joint_only | substitution_only

    joint_count <- sum(joint_pass)
    calibrated_joint_count <- sum(calibrated_joint_pass)
    substitution_count <- sum(substitution_pass)
    aic_shared_count <- sum(aic_shared)
    conflation_count <- sum(conflation)
    aic_shared_transfer_fail <- mean(aic_shared & !substitution_pass)
    transfer_fail_given_aic_shared <- if (aic_shared_count > 0) {
      sum(aic_shared & !substitution_pass) / aic_shared_count
    } else {
      NA_real_
    }
    aic_shared_given_transfer_fail <- if (sum(!substitution_pass) > 0) {
      sum(aic_shared & !substitution_pass) / sum(!substitution_pass)
    } else {
      NA_real_
    }

    sA_true <- 2 - d / 2
    sB_true <- 2 + d / 2
    estimates <- tibble(
      delta = d, nrep = nrep,
      raw_joint_pass_rate = mean(joint_pass),
      raw_joint_wilson_half_width = wilson_interval(joint_count, nrep)[, "half_width"],
      calibrated_joint_pass_rate = mean(calibrated_joint_pass),
      calibrated_joint_wilson_half_width =
        wilson_interval(calibrated_joint_count, nrep)[, "half_width"],
      substitution_pass_rate = mean(substitution_pass),
      substitution_wilson_half_width = wilson_interval(substitution_count, nrep)[, "half_width"],
      conflation_rate = mean(conflation),
      conflation_wilson_half_width = wilson_interval(conflation_count, nrep)[, "half_width"],
      calibrated_joint_substitution_disagreement_rate =
        mean(calibrated_disagreement),
      calibrated_joint_only_rate = mean(calibrated_joint_only),
      substitution_only_rate = mean(substitution_only),
      joint_A_pass_rate =
        mean(replication_matrix[, "R_joint_A"] <= EPS_PRIMARY),
      joint_B_pass_rate =
        mean(replication_matrix[, "R_joint_B"] <= EPS_PRIMARY),
      B_to_A_pass_rate =
        mean(replication_matrix[, "R_B_to_A"] <= EPS_PRIMARY),
      A_to_B_pass_rate =
        mean(replication_matrix[, "R_A_to_B"] <= EPS_PRIMARY),
      aic_shared_rate = mean(aic_shared),
      aic_shared_wilson_half_width = wilson_interval(aic_shared_count, nrep)[, "half_width"],
      aic_shared_transfer_fail_rate = aic_shared_transfer_fail,
      transfer_fail_given_aic_shared = transfer_fail_given_aic_shared,
      aic_shared_given_transfer_fail = aic_shared_given_transfer_fail,
      sA_hat_mean = mean(replication_matrix[, "sA_hat"]),
      sA_bias = mean(replication_matrix[, "sA_hat"]) - sA_true,
      sB_hat_mean = mean(replication_matrix[, "sB_hat"]),
      sB_bias = mean(replication_matrix[, "sB_hat"]) - sB_true,
      sA_rmse = sqrt(mean((replication_matrix[, "sA_hat"] - sA_true)^2)),
      sB_rmse = sqrt(mean((replication_matrix[, "sB_hat"] - sB_true)^2)),
      sA_mcse = sd(replication_matrix[, "sA_hat"]) / sqrt(nrep),
      sB_mcse = sd(replication_matrix[, "sB_hat"]) / sqrt(nrep),
      sJ_between = mean(replication_matrix[, "sJ_between"]),
      joint_pred_rms = mean(replication_matrix[, "joint_probability_rmse"]),
      joint_pred_max = mean(
        replication_matrix[, "joint_probability_max_error"]
      ),
      substitution_pred_rms = mean(
        replication_matrix[, "substitution_probability_rmse"]
      ),
      substitution_pred_max = mean(
        replication_matrix[, "substitution_probability_max_error"]
      ),
      conv_rate = mean(replication_matrix[, c("convA", "convB", "convJ")]),
      rerun_rate = mean(replication_matrix[, "reruns"] > 0)
    )
    out[[i]] <- bind_cols(estimates, boundary_summary(replication_matrix))
    if (i %% 6 == 0) message(sprintf("  done through delta = %.2f", d))
  }
  bind_rows(out)
}

# Sensitivity grid ----

# Supplement S3. A new sample size needs new data; a new tolerance only
# reclassifies the same fits, so each block of fits is scored at all three
# tolerances.
run_sensitivity <- function(nrep = 1500,
                            trials_per_level_values = c(40, 80, 160),
                            deltas = DELTAS) {
  message("Sensitivity grid: ", nrep, " replications at ",
          length(trials_per_level_values), " sample sizes")
  rows <- list()
  k <- 0L
  for (trials_per_level in trials_per_level_values) {
    n_per_task <- trials_per_level * length(X_LEVELS)
    for (d in deltas) {
      replication_matrix <- reps_matrix(d, nrep, m = trials_per_level)
      check_fit_matrix(
        replication_matrix,
        sprintf("sensitivity grid at n=%d, delta=%.2f", n_per_task, d)
      )
      bound_checks <- boundary_summary(replication_matrix)
      for (eps in EPS_GRID) {
        decisions <- decision_flags(replication_matrix, eps, n_per_task)
        calibrated_joint_pass <- decisions$calibrated_joint_pass
        substitution_pass <- decisions$substitution_pass
        k <- k + 1L
        estimates <- tibble(
          n_per_task = n_per_task, m_per_level = trials_per_level, delta = d, eps = eps,
          nrep = nrep,
          raw_joint_pass_rate = mean(decisions$joint_pass),
          calibrated_joint_pass_rate = mean(calibrated_joint_pass),
          substitution_pass_rate = mean(substitution_pass),
          conflation_rate = mean(decisions$joint_pass & !substitution_pass),
          calibrated_joint_substitution_disagreement_rate = mean(
            calibrated_joint_pass != substitution_pass
          ),
          calibrated_joint_only_rate = mean(
            calibrated_joint_pass & !substitution_pass
          ),
          substitution_only_rate = mean(
            !calibrated_joint_pass & substitution_pass
          ),
          aic_shared_rate = mean(decisions$aic_shared),
          bic_shared_rate = mean(decisions$bic_shared),
          lrt_shared_rate = mean(decisions$lrt_shared),
          aic_shared_transfer_fail_rate = mean(
            decisions$aic_shared & !decisions$substitution_pass
          ),
          bic_shared_transfer_fail_rate = mean(
            decisions$bic_shared & !decisions$substitution_pass
          ),
          lrt_shared_transfer_fail_rate = mean(
            decisions$lrt_shared & !decisions$substitution_pass
          ),
          conv_rate = mean(replication_matrix[, c("convA", "convB", "convJ")]),
          rerun_rate = mean(replication_matrix[, "reruns"] > 0)
        )
        rows[[k]] <- bind_cols(estimates, bound_checks)
      }
    }
    message("  done n per task = ", n_per_task)
  }
  bind_rows(rows)
}

# Focused conditions ----

# 10,000 replications at five gaps around the transfer boundary. Every
# replication at delta = 0.50 is saved: the model-selection comparison in the
# paper and the open markers in Figure 4 come from it.
run_focused_conditions <- function(focus_path, nrep = 10000,
                                   deltas = seq(0.3, 0.7, by = 0.1)) {
  message("Focused conditions: ", nrep, " replications at delta = ",
          paste(deltas, collapse = ", "))
  results <- list()
  i <- 0L
  for (d in deltas) {
    replication_matrix <- reps_matrix(d, nrep)
    check_fit_matrix(
      replication_matrix,
      sprintf("focused grid at delta %.2f", d)
    )
    bound_checks <- boundary_summary(replication_matrix)
    if (abs(d - 0.50) < 1e-8) {
      write_csv(focus_diagnostics(replication_matrix, d), focus_path)
    }
    decisions <- decision_flags(replication_matrix, EPS_PRIMARY, N_TASK)
    aic_shared <- decisions$aic_shared
    calibrated_joint_pass <- decisions$calibrated_joint_pass
    substitution_pass <- decisions$substitution_pass
    n_aic_shared <- sum(aic_shared)
    n_shared_transfer_fail <- sum(aic_shared & !substitution_pass)
    interval <- wilson_interval(n_shared_transfer_fail, nrep)
    i <- i + 1L
    estimates <- tibble(
      delta = d, nrep = nrep,
      aic_shared_rate = mean(aic_shared),
      substitution_pass_rate = mean(substitution_pass),
      raw_joint_pass_rate = mean(decisions$joint_pass),
      calibrated_joint_pass_rate = mean(calibrated_joint_pass),
      conflation_rate = mean(decisions$joint_pass & !substitution_pass),
      calibrated_joint_substitution_disagreement_rate = mean(
        calibrated_joint_pass != substitution_pass
      ),
      calibrated_joint_only_rate = mean(
        calibrated_joint_pass & !substitution_pass
      ),
      substitution_only_rate = mean(
        !calibrated_joint_pass & substitution_pass
      ),
      aic_shared_transfer_fail_rate = n_shared_transfer_fail / nrep,
      aic_shared_transfer_fail_lo = interval[, "lower"],
      aic_shared_transfer_fail_hi = interval[, "upper"],
      transfer_fail_given_aic_shared = if (n_aic_shared > 0) {
        n_shared_transfer_fail / n_aic_shared
      } else {
        NA_real_
      },
      conv_rate = mean(
        replication_matrix[, c("convA", "convB", "convJ")]
      ),
      rerun_rate = mean(replication_matrix[, "reruns"] > 0)
    )
    results[[i]] <- bind_cols(estimates, bound_checks)
  }
  bind_rows(results)
}

# Replication-level record for the focused delta = 0.50 condition: losses,
# decisions, the likelihood-ratio statistic, parameter estimates and bound
# checks.
focus_diagnostics <- function(replication_matrix, delta = 0.50) {
  lrt <- replication_matrix[, "lrt"]
  aic_difference_shared_minus_separate <- lrt - AIC_LRT_CUTOFF
  decisions <- decision_flags(replication_matrix, EPS_PRIMARY, N_TASK)
  joint_pass <- decisions$joint_pass
  calibrated_joint_pass <- decisions$calibrated_joint_pass
  substitution_pass <- decisions$substitution_pass
  aic_shared <- decisions$aic_shared
  events <- tibble(
    replication = seq_len(nrow(replication_matrix)),
    delta = delta,
    R_joint = replication_matrix[, "RJ"],
    R_substitution = replication_matrix[, "RS"],
    R_joint_A = replication_matrix[, "R_joint_A"],
    R_joint_B = replication_matrix[, "R_joint_B"],
    R_B_to_A = replication_matrix[, "R_B_to_A"],
    R_A_to_B = replication_matrix[, "R_A_to_B"],
    substitution_tolerance_ratio = replication_matrix[, "RS"] / EPS_PRIMARY,
    lrt = lrt,
    aic_difference_shared_minus_separate =
      aic_difference_shared_minus_separate,
    aic_weight_shared = plogis(-aic_difference_shared_minus_separate / 2),
    joint_pass = joint_pass,
    calibrated_joint_pass = calibrated_joint_pass,
    substitution_pass = substitution_pass,
    joint_A_pass = replication_matrix[, "R_joint_A"] <= EPS_PRIMARY,
    joint_B_pass = replication_matrix[, "R_joint_B"] <= EPS_PRIMARY,
    B_to_A_pass = replication_matrix[, "R_B_to_A"] <= EPS_PRIMARY,
    A_to_B_pass = replication_matrix[, "R_A_to_B"] <= EPS_PRIMARY,
    conflation = joint_pass & !substitution_pass,
    calibrated_joint_substitution_disagreement =
      calibrated_joint_pass != substitution_pass,
    calibrated_joint_only = calibrated_joint_pass & !substitution_pass,
    substitution_only = !calibrated_joint_pass & substitution_pass,
    aic_shared = aic_shared,
    lrt_shared_not_rejected = decisions$lrt_shared,
    bic_shared = decisions$bic_shared,
    aic_joint_pass_transfer_fail =
      aic_shared & joint_pass & !substitution_pass,
    aic_joint_fail_transfer_fail =
      aic_shared & !joint_pass & !substitution_pass,
    convA = as.logical(replication_matrix[, "convA"]),
    convB = as.logical(replication_matrix[, "convB"]),
    convJ = as.logical(replication_matrix[, "convJ"]),
    bndA = as.logical(replication_matrix[, "bndA"]),
    bndB = as.logical(replication_matrix[, "bndB"]),
    bndJ = as.logical(replication_matrix[, "bndJ"]),
    bound_A_to_B_fit = as.logical(replication_matrix[, "bound_A_to_B_fit"]),
    bound_B_to_A_fit = as.logical(replication_matrix[, "bound_B_to_A_fit"]),
    reruns = replication_matrix[, "reruns"]
  )

  parameter_diagnostics <- replication_matrix[, c(
    "sA_hat", "muA_hat", "sB_hat", "muB_hat", "sJ_hat",
    "mu_joint_A_hat", "mu_joint_B_hat",
    "mu_A_to_B_hat", "mu_B_to_A_hat",
    BOUND_FLAG_COLUMNS,
    BOUND_SCALED_DISTANCE_COLUMNS, BOUND_OPTIMIZER_DISTANCE_COLUMNS,
    BOUND_NATIVE_DISTANCE_COLUMNS, "bound_parameter_count",
    "bound_fit_count", "bound_replication"
  ), drop = FALSE] %>%
    as_tibble()

  bind_cols(events, parameter_diagnostics)
}

# Event flags for the saved delta = 0.50 replications. The model-selection
# flags are recomputed from the likelihood-ratio statistic.
derive_focus_events <- function(focus_runs) {
  focus_runs %>%
    mutate(
      joint_pass = as.logical(joint_pass),
      calibrated_joint_pass = R_joint <= EPS_JOINT_CALIBRATED,
      substitution_pass = as.logical(substitution_pass),
      conflation = joint_pass & !substitution_pass,
      calibrated_joint_substitution_disagreement =
        calibrated_joint_pass != substitution_pass,
      calibrated_joint_only =
        calibrated_joint_pass & !substitution_pass,
      substitution_only =
        !calibrated_joint_pass & substitution_pass,
      aic_shared = lrt < AIC_LRT_CUTOFF,
      lrt_shared_not_rejected = lrt < LRT_05_CUTOFF,
      bic_shared = lrt < BIC_LRT_CUTOFF,
      aic_joint_pass_transfer_fail =
        aic_shared & joint_pass & !substitution_pass,
      aic_joint_fail_transfer_fail =
        aic_shared & !joint_pass & !substitution_pass
    )
}

# Proportion of replications with an event, optionally within a subset given
# by `denominator`, with its Wilson interval.
event_estimate <- function(event_id, description, definition, event,
                           denominator = rep(TRUE, length(event))) {
  keep <- as.logical(denominator)
  x <- sum(as.logical(event) & keep)
  n <- sum(keep)
  ci <- if (n > 0) wilson_interval(x, n) else cbind(lower = NA, upper = NA)
  tibble(
    event_id = event_id,
    description = description,
    definition = definition,
    count = x,
    denominator = n,
    estimate = if (n > 0) x / n else NA_real_,
    mc95_low = unname(ci[, "lower"]),
    mc95_high = unname(ci[, "upper"]),
    mc_interval = "95% Wilson"
  )
}

# J = joint criterion passes, S = bidirectional transfer passes, M_AIC = AIC
# prefers the shared slope.
focus_event_summary <- function(focus_runs) {
  d <- derive_focus_events(focus_runs)
  joint_pass <- d$joint_pass
  calibrated_joint_pass <- d$calibrated_joint_pass
  substitution_pass <- d$substitution_pass
  aic_shared <- d$aic_shared

  bind_rows(
    event_estimate("J", "Joint predictor passes", "J", joint_pass),
    event_estimate(
      "J_cal", "Joint predictor passes the calibrated threshold",
      "J_cal: R_joint <= epsilon/4", calibrated_joint_pass
    ),
    event_estimate(
      "S", "Bidirectional transfer passes", "S", substitution_pass
    ),
    event_estimate(
      "J_cal_not_equal_S",
      "Calibrated joint and transfer decisions disagree",
      "J_cal != S", calibrated_joint_pass != substitution_pass
    ),
    event_estimate(
      "J_cal_and_not_S",
      "Calibrated joint passes and transfer fails",
      "J_cal & !S", calibrated_joint_pass & !substitution_pass
    ),
    event_estimate(
      "not_J_cal_and_S",
      "Calibrated joint fails and transfer passes",
      "!J_cal & S", !calibrated_joint_pass & substitution_pass
    ),
    event_estimate(
      "C", "Joint passes and transfer fails", "J & !S",
      joint_pass & !substitution_pass
    ),
    event_estimate(
      "M_AIC", "AIC selects the shared model", "M_AIC", aic_shared
    ),
    event_estimate("M_AIC_and_not_S", "AIC selects shared and transfer fails",
                   "M_AIC & !S", aic_shared & !substitution_pass),
    event_estimate("M_AIC_and_J_and_not_S",
                   "AIC selects shared, joint passes, and transfer fails",
                   "M_AIC & J & !S",
                   aic_shared & joint_pass & !substitution_pass),
    event_estimate("M_AIC_and_not_J_and_not_S",
                   "AIC selects shared, joint fails, and transfer fails",
                   "M_AIC & !J & !S",
                   aic_shared & !joint_pass & !substitution_pass),
    event_estimate("not_M_AIC_and_S",
                   "AIC selects separate and transfer passes",
                   "!M_AIC & S", !aic_shared & substitution_pass),
    event_estimate("AIC_transfer_disagreement",
                   "AIC and transfer disagree in either direction",
                   "(M_AIC & !S) | (!M_AIC & S)",
                   (aic_shared & !substitution_pass) |
                     (!aic_shared & substitution_pass)),
    event_estimate("not_S_given_M_AIC",
                   "Transfer fails given AIC selects shared",
                   "P(!S | M_AIC)", !substitution_pass, aic_shared),
    event_estimate("J_given_M_AIC_and_not_S",
                   "Joint passes given AIC selects shared and transfer fails",
                   "P(J | M_AIC & !S)", joint_pass,
                   aic_shared & !substitution_pass)
  )
}

# 2 x 2 table of two binary decisions, with cell, row and column proportions.
two_by_two <- function(table_id, row_variable, row_values, row_labels,
                       column_variable, column_values, column_labels) {
  n <- length(row_values)
  map_dfr(c(FALSE, TRUE), function(row_level) {
    map_dfr(c(FALSE, TRUE), function(column_level) {
      cell <- row_values == row_level & column_values == column_level
      x <- sum(cell)
      row_total <- sum(row_values == row_level)
      column_total <- sum(column_values == column_level)
      ci <- wilson_interval(x, n)
      tibble(
        table_id = table_id,
        row_variable = row_variable,
        row_level = unname(row_labels[as.character(row_level)]),
        column_variable = column_variable,
        column_level = unname(column_labels[as.character(column_level)]),
        count = x,
        total = n,
        proportion_all = x / n,
        mc95_low = unname(ci[, "lower"]),
        mc95_high = unname(ci[, "upper"]),
        mc_interval = "95% Wilson",
        row_proportion = if (row_total > 0) x / row_total else NA_real_,
        column_proportion =
          if (column_total > 0) x / column_total else NA_real_
      )
    })
  })
}

# Supplementary Table S2 is the AIC-by-transfer table.
focus_contingency_tables <- function(focus_runs) {
  d <- derive_focus_events(focus_runs)
  bind_rows(
    two_by_two(
      "joint_by_transfer", "joint check", d$joint_pass,
      c(`FALSE` = "fail", `TRUE` = "pass"),
      "transfer check", d$substitution_pass,
      c(`FALSE` = "fail", `TRUE` = "pass")
    ),
    two_by_two(
      "calibrated_joint_by_transfer", "calibrated joint check",
      d$calibrated_joint_pass,
      c(`FALSE` = "fail", `TRUE` = "pass"),
      "transfer check", d$substitution_pass,
      c(`FALSE` = "fail", `TRUE` = "pass")
    ),
    two_by_two(
      "aic_by_transfer", "AIC model", d$aic_shared,
      c(`FALSE` = "separate", `TRUE` = "shared"),
      "transfer check", d$substitution_pass,
      c(`FALSE` = "fail", `TRUE` = "pass")
    )
  )
}

# The same events for AIC, the likelihood-ratio test and BIC.
focus_criterion_events <- function(focus_runs) {
  d <- derive_focus_events(focus_runs)
  criteria <- list(
    list(name = "AIC", selected = d$aic_shared,
         rule = sprintf("LRT < %.6f", AIC_LRT_CUTOFF),
         cutoff = AIC_LRT_CUTOFF),
    list(name = "LRT (alpha=.05, df=1)",
         selected = d$lrt_shared_not_rejected,
         rule = sprintf("LRT < chi-square(.95, 1) = %.6f", LRT_05_CUTOFF),
         cutoff = LRT_05_CUTOFF),
    list(name = "BIC", selected = d$bic_shared,
         rule = sprintf("LRT < ln(%d) = %.6f", N_TOTAL, BIC_LRT_CUTOFF),
         cutoff = BIC_LRT_CUTOFF)
  )

  map_dfr(criteria, function(criterion) {
    shared_selected <- criterion$selected
    joint_pass <- d$joint_pass
    substitution_pass <- d$substitution_pass
    criterion_name <- criterion$name
    criterion_rule <- criterion$rule
    criterion_cutoff <- criterion$cutoff
    bind_rows(
      event_estimate("shared", "Criterion selects or retains shared",
                     "M_C", shared_selected),
      event_estimate("shared_and_transfer_fail",
                     "Shared selected and transfer fails", "M_C & !S",
                     shared_selected & !substitution_pass),
      event_estimate("transfer_fail_given_shared",
                     "Transfer fails given shared selected", "P(!S | M_C)",
                     !substitution_pass, shared_selected),
      event_estimate("shared_and_transfer_pass",
                     "Shared selected and transfer passes", "M_C & S",
                     shared_selected & substitution_pass),
      event_estimate("separate_and_transfer_pass",
                     "Separate selected and transfer passes", "!M_C & S",
                     !shared_selected & substitution_pass),
      event_estimate("shared_and_joint_pass_and_transfer_fail",
                     "Shared selected, joint passes, and transfer fails",
                     "M_C & J & !S",
                     shared_selected & joint_pass & !substitution_pass)
    ) %>%
      mutate(
        criterion = criterion_name,
        criterion_rule = criterion_rule,
        lrt_cutoff = criterion_cutoff,
        log_base = "natural",
        n_total = N_TOTAL,
        k_shared = K_SHARED,
        k_separate = K_SEPARATE,
        .before = 1
      )
  })
}

write_focus_event_outputs <- function(focus_runs, output_dir) {
  write_csv(
    focus_event_summary(focus_runs),
    file.path(output_dir, "focus_event_summary.csv")
  )
  write_csv(
    focus_contingency_tables(focus_runs),
    file.path(output_dir, "focus_contingency_tables.csv")
  )
  write_csv(
    focus_criterion_events(focus_runs),
    file.path(output_dir, "focus_criterion_events.csv")
  )
}

# Recovery benchmark ----

# Separate-model recovery of the generating slopes and locations at
# delta = 0.60 (Supplement S3).
run_recovery_benchmark <- function(nrep = 2000, delta = 0.60) {
  replication_matrix <- reps_matrix(delta, nrep)
  check_fit_matrix(
    replication_matrix,
    sprintf("recovery check at delta %.2f", delta)
  )
  bound_checks <- boundary_summary(replication_matrix)
  sA_true <- 2 - delta / 2
  sB_true <- 2 + delta / 2
  recovery <- tibble(
    delta = delta, replications = nrep,
    parameter = c("s_A", "s_B", "mu_A", "mu_B"),
    true = c(sA_true, sB_true, MU_A, MU_B),
    estimate = c(
      mean(replication_matrix[, "sA_hat"]),
      mean(replication_matrix[, "sB_hat"]),
      mean(replication_matrix[, "muA_hat"]),
      mean(replication_matrix[, "muB_hat"])
    ),
    bias = c(
      mean(replication_matrix[, "sA_hat"]) - sA_true,
      mean(replication_matrix[, "sB_hat"]) - sB_true,
      mean(replication_matrix[, "muA_hat"]) - MU_A,
      mean(replication_matrix[, "muB_hat"]) - MU_B
    ),
    rmse = c(
      sqrt(mean((replication_matrix[, "sA_hat"] - sA_true)^2)),
      sqrt(mean((replication_matrix[, "sB_hat"] - sB_true)^2)),
      sqrt(mean((replication_matrix[, "muA_hat"] - MU_A)^2)),
      sqrt(mean((replication_matrix[, "muB_hat"] - MU_B)^2))
    ),
    conv_rate = mean(
      replication_matrix[, c("convA", "convB", "convJ")]
    ),
    boundary_rate = bound_checks$boundary_rate,
    boundary_parameter_rate = bound_checks$boundary_parameter_rate,
    boundary_replication_rate = bound_checks$boundary_replication_rate,
    min_scaled_bound_distance = bound_checks$min_scaled_bound_distance
  )
  audit_grid <- bind_cols(tibble(nrep = nrep), bound_checks)
  attr(recovery, "parameter_boundary_audit") <-
    parameter_boundary_audit_from_grid(audit_grid)
  recovery
}

# Optimizer checks ----

# The flag rule for each parameter, on the optimizer scale and on the
# parameter's own scale. On the log-slope scale a fixed tolerance corresponds
# to different slope distances at the two bounds.
parameter_boundary_specification <- function() {
  slope_parameters <- BOUND_PARAMETER_IDS %in%
    c("sep_A_s", "sep_B_s", "joint_s")
  optimizer_lower <- ifelse(slope_parameters, LOG_S_LO, MU_LO)
  optimizer_upper <- ifelse(slope_parameters, LOG_S_HI, MU_HI)
  native_lower <- ifelse(slope_parameters, S_LO, MU_LO)
  native_upper <- ifelse(slope_parameters, S_HI, MU_HI)
  optimizer_tolerance <-
    BOUND_SCALED_TOL * (optimizer_upper - optimizer_lower)

  native_lower_tolerance <- ifelse(
    slope_parameters,
    exp(optimizer_lower + optimizer_tolerance) - native_lower,
    optimizer_tolerance
  )
  native_upper_tolerance <- ifelse(
    slope_parameters,
    native_upper - exp(optimizer_upper - optimizer_tolerance),
    optimizer_tolerance
  )

  tibble(
    parameter = BOUND_PARAMETER_IDS,
    component_fit = c(
      "separate A", "separate A", "separate B", "separate B",
      "joint", "joint", "joint", "substitution A to B", "substitution B to A"
    ),
    role = c(
      "slope", "location", "slope", "location",
      "shared slope", "location A", "location B",
      "target B location", "target A location"
    ),
    optimizer_coordinate = ifelse(slope_parameters, "log(slope)", "location"),
    optimizer_lower_bound = optimizer_lower,
    optimizer_upper_bound = optimizer_upper,
    native_lower_bound = native_lower,
    native_upper_bound = native_upper,
    scaled_tolerance = BOUND_SCALED_TOL,
    optimizer_coordinate_tolerance = optimizer_tolerance,
    native_tolerance_at_lower_bound = native_lower_tolerance,
    native_tolerance_at_upper_bound = native_upper_tolerance,
    flag_rule = paste0(
      "distance to nearer optimizer bound / optimizer interval <= ",
      format(BOUND_SCALED_TOL, scientific = TRUE)
    )
  )
}

# One row per parameter: how often it was flagged and how close it came to a
# bound, from rows of grid summaries (one row per simulated condition).
parameter_boundary_audit_from_grid <- function(grid) {
  denominator <- sum(grid$nrep)
  map_dfr(BOUND_PARAMETER_IDS, function(id) {
    tibble(
      parameter = id,
      count_near_bound = sum(grid[[paste0("bound_count_", id)]]),
      denominator = denominator,
      minimum_scaled_distance =
        min(grid[[paste0("min_scaled_bound_distance_", id)]]),
      minimum_optimizer_coordinate_distance =
        min(grid[[paste0("min_optimizer_bound_distance_", id)]]),
      minimum_native_distance =
        min(grid[[paste0("min_native_bound_distance_", id)]])
    )
  }) %>%
    left_join(parameter_boundary_specification(), by = "parameter") %>%
    mutate(proportion_near_bound = count_near_bound / denominator, .after = denominator)
}

# The same, from the saved delta = 0.50 replications.
parameter_boundary_audit_from_focus <- function(focus_runs) {
  map_dfr(seq_along(BOUND_PARAMETER_IDS), function(i) {
    tibble(
      parameter = BOUND_PARAMETER_IDS[i],
      count_near_bound = sum(as.logical(focus_runs[[BOUND_FLAG_COLUMNS[i]]])),
      denominator = nrow(focus_runs),
      minimum_scaled_distance =
        min(focus_runs[[BOUND_SCALED_DISTANCE_COLUMNS[i]]]),
      minimum_optimizer_coordinate_distance =
        min(focus_runs[[BOUND_OPTIMIZER_DISTANCE_COLUMNS[i]]]),
      minimum_native_distance =
        min(focus_runs[[BOUND_NATIVE_DISTANCE_COLUMNS[i]]])
    )
  }) %>%
    left_join(parameter_boundary_specification(), by = "parameter") %>%
    mutate(proportion_near_bound = count_near_bound / denominator, .after = denominator)
}

# The three tolerance rows of the sensitivity grid share their fits, so each
# simulated condition is counted once.
sensitivity_boundary_audit <- function(sensitivity) {
  sensitivity %>%
    distinct(analysis, n_per_task, m_per_level, delta, nrep, .keep_all = TRUE) %>%
    group_split(analysis, n_per_task, m_per_level) %>%
    map_dfr(function(condition) {
      parameter_boundary_audit_from_grid(condition) %>%
        mutate(
          analysis = condition$analysis[1],
          n_per_task = condition$n_per_task[1],
          m_per_level = condition$m_per_level[1],
          .before = 1
        )
    })
}

write_parameter_boundary_checks <- function(finite_grid, focus_runs, recovery,
                                            sensitivity, output_dir) {
  add_condition <- function(audit, analysis, n_per_task, m_per_level) {
    audit %>%
      mutate(
        analysis = analysis,
        n_per_task = n_per_task,
        m_per_level = m_per_level,
        .before = 1
      )
  }

  sensitivity_audit <- sensitivity_boundary_audit(sensitivity) %>%
    mutate(
      analysis = recode(analysis, grid = "sensitivity_grid")
    )

  checks <- bind_rows(
    add_condition(
      parameter_boundary_audit_from_grid(finite_grid),
      "main_grid", N_TASK, M_TRAIN
    ),
    add_condition(
      parameter_boundary_audit_from_focus(focus_runs),
      "focus_delta_050", N_TASK, M_TRAIN
    ),
    add_condition(
      attr(recovery, "parameter_boundary_audit"),
      "recovery_delta_060", N_TASK, M_TRAIN
    ),
    sensitivity_audit
  )
  write_csv(checks, file.path(output_dir, "parameter_boundary_checks.csv"))
}

# The per-parameter bound columns go to parameter_boundary_checks.csv instead
# of the grid summaries.
compact_simulation_summary <- function(data) {
  detailed_columns <- unlist(lapply(BOUND_PARAMETER_IDS, function(parameter) {
    paste0(
      c(
        "bound_count_", "bound_rate_", "min_scaled_bound_distance_",
        "min_optimizer_bound_distance_", "min_native_bound_distance_"
      ),
      parameter
    )
  }))
  select(data, -any_of(detailed_columns))
}

# Restarts, convergence and bound checks for the main grid (Supplementary
# Table S1). Each replication has three L-BFGS-B fits with a convergence code
# and five fitting operations checked for bounds.
write_convergence_audit <- function(finite_grid, output_dir) {
  n_reps <- sum(finite_grid$nrep)
  n_convergence_fits <- 3L * n_reps
  n_boundary_fits <- 5L * n_reps
  n_parameters <- 9L * n_reps
  retried_reps <- round(sum(finite_grid$nrep * finite_grid$rerun_rate))
  converged_fits <- round(sum(
    3L * finite_grid$nrep * finite_grid$conv_rate
  ))
  boundary_fits <- sum(finite_grid$boundary_fit_count)
  boundary_parameters <- sum(finite_grid$boundary_parameter_count)
  boundary_reps <- sum(finite_grid$boundary_replication_count)
  min_scaled_distance <- min(finite_grid$min_scaled_bound_distance)

  audit <- tibble(
    metric = c(
      "replications_with_any_retry",
      "component_fits_converged_after_retry",
      "component_fits_not_converged_after_retry",
      "fitting_operations_with_parameter_near_boundary",
      "free_parameters_checked",
      "free_parameters_near_optimizer_boundary",
      "replications_with_any_parameter_near_boundary",
      "replications_discarded_or_replaced"
    ),
    count = c(
      retried_reps,
      converged_fits,
      n_convergence_fits - converged_fits,
      boundary_fits,
      n_parameters,
      boundary_parameters,
      boundary_reps,
      0L
    ),
    denominator = c(
      n_reps, n_convergence_fits, n_convergence_fits, n_boundary_fits,
      n_parameters, n_parameters, n_reps, n_reps
    ),
    description = c(
      "At least one component fit triggered the fixed restart grid",
      "Final convergence status across the two separate and one joint fit",
      "Component fits not converged after the fixed retry rule",
      paste(
        "At least one free parameter was within the declared scaled tolerance",
        "in any of five fitting operations"
      ),
      "Nine free parameters were checked in every replication",
      "Free parameters within 0.01% of their bounded optimizer interval",
      "At least one of nine free parameters met the bound flag rule",
      "All generated replications were retained"
    )
  ) %>%
    mutate(
      proportion = count / denominator,
      scaled_boundary_tolerance = BOUND_SCALED_TOL,
      minimum_scaled_boundary_distance = min_scaled_distance
    )

  write_csv(audit, file.path(output_dir, "finite_sample_convergence_audit.csv"))
}

# Full simulation ----

simulate_2a <- function(nrep_main = 2000, nrep_sens = 1500,
                                  nrep_focused = 10000,
                                  output_dir = file.path(results_dir, "simulation2a")) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(SEED)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("Population analysis")
  population <- run_population()
  write_csv(population_output_table(population$pop, population$cross),
            file.path(output_dir, "population_results.csv"))
  write_csv(reference_tolerance_anchor(),
            file.path(output_dir, "reference_tolerance_anchor.csv"))

  message("Recovery benchmark")
  recovery <- run_recovery_benchmark(nrep_main, 0.60)
  write_csv(recovery, file.path(output_dir, "recovery_benchmark.csv"))

  finite_summary <- run_finite_main(nrep_main)
  write_csv(
    compact_simulation_summary(finite_summary) %>%
      rename(replications = nrep),
    file.path(output_dir, "finite_sample_by_delta.csv")
  )
  write_convergence_audit(finite_summary, output_dir)

  sensitivity <- run_sensitivity(nrep_sens)
  focus_path <- file.path(output_dir, "focus_delta_050.csv.gz")
  focused_summary <- run_focused_conditions(focus_path, nrep_focused)
  focus_runs <- read_csv(focus_path, show_col_types = FALSE)
  write_focus_event_outputs(focus_runs, output_dir)

  # sensitivity.csv holds the sensitivity grid and the focused conditions.
  focused_block <- focused_summary %>%
    mutate(
      analysis = "focused_grid", n_per_task = N_TASK,
      m_per_level = M_TRAIN, eps = EPS_PRIMARY,
      bic_shared_rate = NA_real_, lrt_shared_rate = NA_real_,
      bic_shared_transfer_fail_rate = NA_real_, lrt_shared_transfer_fail_rate = NA_real_
    ) %>%
    select(
      analysis, n_per_task, m_per_level, delta, eps, nrep,
      raw_joint_pass_rate, calibrated_joint_pass_rate,
      substitution_pass_rate, conflation_rate,
      calibrated_joint_substitution_disagreement_rate,
      calibrated_joint_only_rate,
      substitution_only_rate,
      aic_shared_rate, bic_shared_rate, lrt_shared_rate,
      aic_shared_transfer_fail_rate, bic_shared_transfer_fail_rate, lrt_shared_transfer_fail_rate,
      conv_rate, rerun_rate, starts_with("bound"), starts_with("min_")
    )
  sensitivity_output <- bind_rows(
    sensitivity %>% mutate(analysis = "grid"), focused_block
  )
  write_csv(
    compact_simulation_summary(sensitivity_output) %>%
      rename(replications = nrep),
    file.path(output_dir, "sensitivity.csv")
  )
  write_parameter_boundary_checks(
    finite_summary, focus_runs, recovery, sensitivity_output, output_dir
  )
}

# Figure 4 ----

simulation2a_theme <- function(base = 9) {
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title = element_text(size = base + 0.5, hjust = 0, face = "plain",
                                     lineheight = 0.9, margin = margin(b = 4)),
      axis.title = element_text(size = base),
      axis.text = element_text(size = base - 1, colour = "black"),
      legend.position = "bottom",
      legend.text = element_text(size = base - 1),
      legend.key.width = unit(1.1, "cm"),
      legend.margin = margin(t = 1),
      plot.tag = element_text(size = base + 2, face = "bold"),
      plot.margin = margin(5, 7, 4, 5),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA)
    )
}

# Pass rates on the main grid, the population band in which joint adequacy
# holds but transfer fails, and open markers for the focused delta = 0.50 run.
make_figure4 <- function(input_dir = file.path(results_dir, "simulation2a")) {
  cross <- read_csv(file.path(input_dir, "population_results.csv"),
                    show_col_types = FALSE) %>%
    filter(kind == "crossing", abs(eps - EPS_PRIMARY) < 1e-8)
  stopifnot(nrow(cross) == 1L)
  delta_sub <- cross$substitution[[1]]
  finite_summary <- read_csv(file.path(input_dir, "finite_sample_by_delta.csv"),
                             show_col_types = FALSE)
  focused_summary <- read_csv(file.path(input_dir, "sensitivity.csv"),
                              show_col_types = FALSE) %>%
    filter(analysis == "focused_grid")

  method_colours <- c(
    "Joint adequacy" = orange,
    "Bidirectional transfer under substitution" = blue,
    "Joint adequacy without transfer" = purple
  )
  method_shapes <- c(
    "Joint adequacy" = 16,
    "Bidirectional transfer under substitution" = 17,
    "Joint adequacy without transfer" = 15
  )
  method_linetypes <- c(
    "Joint adequacy" = "solid",
    "Bidirectional transfer under substitution" = "22",
    "Joint adequacy without transfer" = "31"
  )

  rates <- finite_summary %>%
    filter(delta <= 1.0) %>%
    transmute(delta,
              `Joint adequacy` = 100 * raw_joint_pass_rate,
              `Bidirectional transfer under substitution` = 100 * substitution_pass_rate,
              `Joint adequacy without transfer` = 100 * conflation_rate) %>%
    pivot_longer(-delta, names_to = "method", values_to = "pass_rate") %>%
    mutate(method = factor(method, levels = names(method_colours)))
  rate_points <- rates %>%
    filter(abs(delta * 10 - round(delta * 10)) < 1e-8)

  focus_rates <- focused_summary %>%
    filter(abs(delta - 0.50) < 1e-8) %>%
    transmute(delta,
              `Joint adequacy` = 100 * raw_joint_pass_rate,
              `Bidirectional transfer under substitution` = 100 * substitution_pass_rate,
              `Joint adequacy without transfer` = 100 * conflation_rate) %>%
    pivot_longer(-delta, names_to = "method", values_to = "pass_rate") %>%
    mutate(method = factor(method, levels = names(method_colours)))
  focus_annotations <- focus_rates %>%
    mutate(
      x = delta + 0.025,
      y = case_when(
        method == "Joint adequacy" ~ pass_rate + 4,
        method == "Bidirectional transfer under substitution" ~ pass_rate + 6,
        TRUE ~ pass_rate - 5
      ),
      label = sprintf("%.2f%%", pass_rate)
    )

  figure <- ggplot(rates, aes(delta, pass_rate, colour = method)) +
    annotate("rect", xmin = cross$substitution, xmax = cross$joint,
             ymin = -Inf, ymax = Inf, fill = "grey92") +
    annotate("text",
             x = cross$joint - 0.21, y = 26,
             label = "Joint adequacy\nwithout transfer",
             colour = "grey35", size = 2.8) +
    geom_line(aes(linetype = method), linewidth = 1.1) +
    geom_point(data = rate_points, aes(shape = method), size = 1.9) +
    geom_point(data = filter(focus_rates, method == "Joint adequacy"),
               shape = 21, size = 3.4, stroke = 1.0, fill = "white",
               show.legend = FALSE) +
    geom_point(data = filter(focus_rates, method == "Bidirectional transfer under substitution"),
               shape = 24, size = 3.4, stroke = 1.0, fill = "white",
               show.legend = FALSE) +
    geom_point(data = filter(focus_rates, method == "Joint adequacy without transfer"),
               shape = 22, size = 3.4, stroke = 1.0, fill = "white",
               show.legend = FALSE) +
    geom_text(
      data = focus_annotations,
      aes(x = x, y = y, label = label, colour = method),
      inherit.aes = FALSE,
      hjust = 0, fontface = "bold", size = 2.8,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = method_colours) +
    scale_shape_manual(values = method_shapes) +
    scale_linetype_manual(values = method_linetypes) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.2),
      limits = c(0, 1),
      expand = expansion(add = c(0.01, 0.01)),
      sec.axis = sec_axis(
        ~ . / delta_sub,
        breaks = c(0, 0.5, 1, 1.5, 2),
        name = expression("Tolerance-normalized difference, " * q == delta / delta[sub])
      )
    ) +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(-3, 102),
                       expand = c(0, 0)) +
    labs(title = "Joint adequacy can persist when transfer fails",
         x = expression("True slope gap " * delta),
         y = "Replications with each outcome (%)",
         colour = NULL, shape = NULL, linetype = NULL) +
    simulation2a_theme(10.5) +
    theme(
      plot.title = element_text(size = 11.5, face = "bold", margin = margin(b = 4)),
      legend.position = "top",
      legend.justification = "left",
      legend.box.just = "left",
      legend.key.width = unit(0.62, "cm"),
      legend.margin = margin(t = 0, b = 1),
      axis.title.x.top = element_text(size = 9.5, margin = margin(b = 3)),
      axis.text.x.top = element_text(size = 8.5),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 8, 4, 6)
    ) +
    guides(colour = guide_legend(
             nrow = 1,
             override.aes = list(linewidth = 1.2, size = 2.3,
                                 shape = c(16, 17, 15),
                                 linetype = c("solid", "22", "31"))),
           shape = "none", linetype = "none")

  save_figure(figure, file.path(figure_dir, "figure4_psychometric_simulation.pdf"),
              height_mm = 82)
}

run_simulation2a <- function(resimulate = FALSE) {
  if (resimulate) simulate_2a()
  make_figure4()
}

if (sys.nframe() == 0) {
  run_simulation2a(resimulate = "--resimulate" %in% commandArgs(trailingOnly = TRUE))
}
