# Probability-scale criterion (Supplement S3 and S4).
#
# The main analyses judge a fitted curve by its mean KL divergence from the
# generating response probabilities. Here the same population fits are judged
# by the largest change in any predicted probability, in percentage points.
# Slopes and locations are still chosen by minimizing the KL loss; only the
# acceptance rule changes. Population calculations only; no random numbers are
# drawn.
#
#   Rscript --vanilla code/R/probability_scale.R

source("code/R/common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})
source("code/R/psychometric_model.R")

output_dir <- file.path(results_dir, "probability_scale")

# Settings ----

PP_CRITERIA <- c(4, 5, 6)                        # largest change, in points
CONCENTRATED_LEVELS <- seq(-0.4, 0.4, by = 0.1)  # Task B grid of Simulation 2B
FIGURE4_DELTA <- 0.70                            # gap discussed with Figure 4

# Stimulus levels of the two tasks. In the reversed design the generating
# slopes are swapped, so Task A has the steeper one.
DESIGNS <- list(
  simulation2a = list(x_A = X_LEVELS, x_B = X_LEVELS, orientation = 1),
  simulation2b = list(
    x_A = X_LEVELS, x_B = CONCENTRATED_LEVELS, orientation = 1
  ),
  simulation2b_reversed = list(
    x_A = X_LEVELS, x_B = CONCENTRATED_LEVELS, orientation = -1
  )
)

# Population fits ----

# Best location for a slope in one task, with the excess loss and the largest
# change in predicted probability that it leaves.
fit_slope <- function(slope, p0, x) {
  best <- best_mu_pop(slope, p0, x)
  p <- pfun(slope, best$mu, x)
  list(R = best$R, pp = 100 * max(abs(p - p0)), p = p)
}

# Joint fit and both transfer directions at slope gap delta, as in
# Simulations 2A and 2B.
population_fits <- function(delta, design) {
  x_A <- design$x_A
  x_B <- design$x_B
  slope_A <- 2 - design$orientation * delta / 2
  slope_B <- 2 + design$orientation * delta / 2
  p0_A <- pfun(slope_A, MU_A, x_A)
  p0_B <- pfun(slope_B, MU_B, x_B)

  joint_slope <- optimize(
    function(slope) {
      fit_slope(slope, p0_A, x_A)$R + fit_slope(slope, p0_B, x_B)$R
    },
    interval = c(S_LO, S_HI), tol = 1e-10
  )$minimum
  joint_A <- fit_slope(joint_slope, p0_A, x_A)
  joint_B <- fit_slope(joint_slope, p0_B, x_B)
  A_to_B <- fit_slope(slope_A, p0_B, x_B)  # Task A slope used in Task B
  B_to_A <- fit_slope(slope_B, p0_A, x_A)  # Task B slope used in Task A

  list(
    R_joint = max(joint_A$R, joint_B$R),
    R_A_to_B = A_to_B$R,
    R_B_to_A = B_to_A$R,
    pp_joint = max(joint_A$pp, joint_B$pp),
    pp_A_to_B = A_to_B$pp,
    pp_B_to_A = B_to_A$pp,
    p0_A = p0_A, p0_B = p0_B,
    p_A_to_B = A_to_B$p, p_B_to_A = B_to_A$p
  )
}

# Slope gap at which a measure reaches a limit.
boundary <- function(design, measure, limit) {
  uniroot(
    function(delta) population_fits(delta, design)[[measure]] - limit,
    interval = c(1e-6, 3), tol = 1e-10
  )$root
}

# Boundaries of each criterion under the KL tolerance and under each
# probability limit. Bidirectional transfer fails at the nearer of the two
# directional boundaries.
boundary_table <- function(name, design) {
  criteria <- tibble(
    criterion = c("kl", rep("pp", length(PP_CRITERIA))),
    limit = c(EPS_PRIMARY, PP_CRITERIA)
  )
  pmap_dfr(criteria, function(criterion, limit) {
    prefix <- if (criterion == "kl") "R_" else "pp_"
    A_to_B <- boundary(design, paste0(prefix, "A_to_B"), limit)
    B_to_A <- boundary(design, paste0(prefix, "B_to_A"), limit)
    joint <- boundary(design, paste0(prefix, "joint"), limit)
    tibble(
      design = name,
      criterion = criterion,
      limit = limit,
      boundary_A_to_B = A_to_B,
      boundary_B_to_A = B_to_A,
      boundary_bidirectional = min(A_to_B, B_to_A),
      boundary_joint = joint,
      joint_over_bidirectional = joint / min(A_to_B, B_to_A),
      A_to_B_over_B_to_A = A_to_B / B_to_A
    )
  })
}

# What a transfer fit changes in its target: the largest change in predicted
# probability, the range of the generating probabilities, how many of them lie
# below .10 or above .90, and the excess loss per unit of mean squared
# probability change.
direction_summary <- function(direction, target_task, R, p, p0) {
  tibble(
    direction = direction,
    target_task = target_task,
    excess_loss = R,
    max_probability_change_pp = 100 * max(abs(p - p0)),
    p0_min = min(p0),
    p0_max = max(p0),
    levels_below_10_or_above_90 = sum(p0 < 0.10 | p0 > 0.90),
    kl_per_mean_squared_change = R / mean((p - p0)^2)
  )
}

changes_at <- function(name, design, point, delta) {
  fits <- population_fits(delta, design)
  bind_rows(
    direction_summary("A_to_B", "B", fits$R_A_to_B, fits$p_A_to_B, fits$p0_B),
    direction_summary("B_to_A", "A", fits$R_B_to_A, fits$p_B_to_A, fits$p0_A)
  ) %>%
    mutate(design = name, point = point, delta = delta, .before = 1)
}

# Analysis ----

run_probability_scale <- function() {
  boundaries <- imap_dfr(DESIGNS, function(design, name) {
    boundary_table(name, design)
  })

  kl_2a <- filter(boundaries, design == "simulation2a", criterion == "kl")
  kl_2b <- filter(boundaries, design == "simulation2b", criterion == "kl")
  changes <- bind_rows(
    changes_at("simulation2a", DESIGNS$simulation2a,
               "kl_bidirectional_boundary", kl_2a$boundary_bidirectional),
    changes_at("simulation2b", DESIGNS$simulation2b,
               "kl_A_to_B_boundary", kl_2b$boundary_A_to_B),
    changes_at("simulation2b", DESIGNS$simulation2b,
               "kl_B_to_A_boundary", kl_2b$boundary_B_to_A),
    changes_at("simulation2b", DESIGNS$simulation2b,
               "figure4_gap", FIGURE4_DELTA)
  )

  # Simulation 2B at its two KL boundaries, each the gap at which one transfer
  # direction reaches the tolerance. The squared ratio of the largest changes
  # in the two targets is the ratio of their KL cost per unit of mean squared
  # change times a factor for how the changes are spread across levels.
  in_A <- filter(changes, design == "simulation2b",
                 point == "kl_B_to_A_boundary", direction == "B_to_A")
  in_B <- filter(changes, design == "simulation2b",
                 point == "kl_A_to_B_boundary", direction == "A_to_B")
  squared_ratio <- (in_B$max_probability_change_pp /
                      in_A$max_probability_change_pp)^2
  kl_ratio <- in_A$kl_per_mean_squared_change /
    in_B$kl_per_mean_squared_change
  decomposition <- tibble(
    design = "simulation2b",
    max_change_A_pp = in_A$max_probability_change_pp,
    max_change_B_pp = in_B$max_probability_change_pp,
    squared_change_ratio = squared_ratio,
    kl_ratio = kl_ratio,
    spread_factor = squared_ratio / kl_ratio
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(boundaries, file.path(output_dir, "boundaries.csv"))
  write_csv(changes, file.path(output_dir, "changes.csv"))
  write_csv(decomposition, file.path(output_dir, "scale_comparison.csv"))
  message("Saved ", output_dir)
}

if (sys.nframe() == 0) {
  run_probability_scale()
}
