# Supplementary checks:
#   - Figure S1 (Supplement S3): the one-quarter calibration applied to the
#     saved Simulation 1 and 2A results;
#   - Supplement S2: Simulation 1 repeated with three independent estimators;
#   - Supplement S5: nuisance compensation on a one-sided stimulus grid.
#
#   Rscript --vanilla code/R/supplement_checks.R

source("code/R/common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})
source("code/R/psychometric_model.R")

output_dir <- file.path(results_dir, "supplement_checks")

# Figure S1 ----

make_figure_s1 <- function() {
  gaussian <- read_csv(
    file.path(results_dir, "simulation1", "simulation1_curves.csv"),
    show_col_types = FALSE
  )
  psychometric <- read_csv(
    file.path(results_dir, "simulation2a", "finite_sample_by_delta.csv"),
    show_col_types = FALSE
  )

  method_order <- c(
    "Joint adequacy", "Calibrated joint adequacy",
    "Bidirectional transfer\nunder substitution"
  )
  finite_curves <- bind_rows(
    transmute(
      gaussian, model = "Gaussian means", delta,
      method = "Joint adequacy", pass_probability = joint_p
    ),
    transmute(
      gaussian, model = "Gaussian means", delta,
      method = "Calibrated joint adequacy", pass_probability = joint_cal_p
    ),
    transmute(
      gaussian, model = "Gaussian means", delta,
      method = "Bidirectional transfer\nunder substitution", pass_probability = substitution_p
    ),
    transmute(
      psychometric, model = "Psychometric slopes", delta,
      method = "Joint adequacy", pass_probability = raw_joint_pass_rate
    ),
    transmute(
      psychometric, model = "Psychometric slopes", delta,
      method = "Calibrated joint adequacy",
      pass_probability = calibrated_joint_pass_rate
    ),
    transmute(
      psychometric, model = "Psychometric slopes", delta,
      method = "Bidirectional transfer\nunder substitution",
      pass_probability = substitution_pass_rate
    )
  ) %>%
    filter(model != "Psychometric slopes" | delta <= 1) %>%
    mutate(
      method = factor(method, levels = method_order),
      pass_rate_percent = 100 * pass_probability
    ) %>%
    arrange(model, method, delta)

  method_colours <- c(
    "Joint adequacy" = orange,
    "Calibrated joint adequacy" = green,
    "Bidirectional transfer\nunder substitution" = blue
  )
  method_linetypes <- c(
    "Joint adequacy" = "solid",
    "Calibrated joint adequacy" = "dotdash",
    "Bidirectional transfer\nunder substitution" = "22"
  )
  calibration_theme <- paper_theme(base_size = 10.5) +
    theme(
      plot.title = element_text(
        size = 11.5, face = "bold", hjust = 0, margin = margin(b = 5)
      ),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.key.width = unit(1.25, "cm"),
      legend.margin = margin(t = 2),
      plot.margin = margin(5, 8, 4, 6)
    )

  make_panel <- function(data, title, x_title, x_breaks, x_limits,
                         show_y_title = TRUE) {
    ggplot(data, aes(delta, pass_rate_percent, colour = method,
                     linetype = method)) +
      geom_line(linewidth = 1.05, lineend = "round") +
      scale_colour_manual(values = method_colours, breaks = method_order,
                          drop = FALSE) +
      scale_linetype_manual(values = method_linetypes, breaks = method_order,
                            drop = FALSE) +
      scale_x_continuous(
        breaks = x_breaks, limits = x_limits,
        expand = expansion(mult = c(0, 0.01))
      ) +
      scale_y_continuous(
        breaks = c(0, 25, 50, 75, 100), limits = c(0, 100),
        expand = expansion(mult = c(0.01, 0.02))
      ) +
      labs(
        title = title, x = x_title,
        y = if (show_y_title) "Replications passing (%)" else NULL,
        colour = NULL, linetype = NULL
      ) +
      calibration_theme +
      guides(
        colour = guide_legend(
          nrow = 1,
          override.aes = list(
            linewidth = 1.2, linetype = unname(method_linetypes)
          )
        ),
        linetype = "none"
      )
  }

  gaussian_plot <- make_panel(
    filter(finite_curves, model == "Gaussian means"),
    "A   Gaussian means", expression("True mean gap " * delta),
    seq(0, 0.6, by = 0.1), c(0, 0.6)
  )
  psychometric_plot <- make_panel(
    filter(finite_curves, model == "Psychometric slopes"),
    "B   Psychometric slopes", expression("True slope gap " * delta),
    seq(0, 1, by = 0.2), c(0, 1), show_y_title = FALSE
  )
  figure <- (gaussian_plot + psychometric_plot +
    plot_layout(guides = "collect", widths = c(1, 1))) &
    theme(legend.position = "bottom")

  save_figure(
    figure, file.path(supplement_figure_dir, "figureS1_calibration_check.pdf"),
    height_mm = 78
  )
}

# Independent estimators (Supplement S2) ----

# Simulation 1 again, but with the joint and the two task estimates drawn
# independently, each with standard error 0.025.
run_independent_gaussian_check <- function() {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(20260824)
  tolerance <- 0.02
  radius <- sqrt(2 * tolerance)
  standard_error <- 0.025
  replications <- 20000L

  results <- map_dfr(seq(0, 0.6, by = 0.05), function(delta) {
    theta_A <- -delta / 2
    theta_B <- delta / 2
    estimate_A <- rnorm(replications, theta_A, standard_error)
    estimate_B <- rnorm(replications, theta_B, standard_error)
    estimate_joint <- rnorm(
      replications, (theta_A + theta_B) / 2, standard_error
    )
    joint_pass <- abs(estimate_joint - theta_A) <= radius &
      abs(estimate_joint - theta_B) <= radius
    substitution_pass <- abs(estimate_A - theta_B) <= radius &
      abs(estimate_B - theta_A) <= radius
    tibble(
      delta,
      replications = replications,
      raw_joint_pass_rate = mean(joint_pass),
      substitution_pass_rate = mean(substitution_pass),
      conflation_rate = mean(joint_pass & !substitution_pass)
    )
  })
  write_csv(results, file.path(output_dir, "gaussian_independent.csv"))
}

# Nuisance compensation (Supplement S5) ----

# A -> B transfer on a one-sided Task B grid, with the target location either
# held at its true value or re-optimized for each transferred slope. gamma is
# the slope-location coupling of the target loss at s = 2, from a numerical
# Hessian; the local prediction for the widening is 1 / sqrt(1 - gamma^2).
run_nuisance_policy_check <- function() {
  target_location <- 0.15
  stimulus_levels <- seq(0.2, 1.8, by = 0.2)
  loss <- function(slope, location, reference) {
    R_T(slope, location, reference, x = stimulus_levels)
  }

  step <- 1e-4
  reference <- pfun(2, target_location, x = stimulus_levels)
  local_risk <- function(slope, location) loss(slope, location, reference)
  h_ss <- (
    local_risk(2 + step, target_location) -
      2 * local_risk(2, target_location) +
      local_risk(2 - step, target_location)
  ) / step^2
  h_mm <- (
    local_risk(2, target_location + step) -
      2 * local_risk(2, target_location) +
      local_risk(2, target_location - step)
  ) / step^2
  h_sm <- (
    local_risk(2 + step, target_location + step) -
      local_risk(2 + step, target_location - step) -
      local_risk(2 - step, target_location + step) +
      local_risk(2 - step, target_location - step)
  ) / (4 * step^2)
  gamma <- h_sm / sqrt(h_ss * h_mm)

  # Excess loss of A -> B transfer at gap delta under each nuisance policy.
  profiled_risk <- function(delta) {
    target <- pfun(2 + delta / 2, target_location, x = stimulus_levels)
    best_mu_pop(2 - delta / 2, target, x = stimulus_levels)$R
  }
  fixed_risk <- function(delta) {
    target <- pfun(2 + delta / 2, target_location, x = stimulus_levels)
    loss(2 - delta / 2, target_location, target)
  }
  profiled_boundary <- uniroot(
    function(delta) profiled_risk(delta) - EPS_PRIMARY,
    c(1e-6, 4), tol = 1e-10
  )$root
  fixed_boundary <- uniroot(
    function(delta) fixed_risk(delta) - EPS_PRIMARY,
    c(1e-6, 4), tol = 1e-10
  )$root

  write_csv(
    tibble(
      design = "one-sided seq(0.2, 1.8, 0.2)",
      gamma = gamma,
      widening_predicted = 1 / sqrt(1 - gamma^2),
      boundary_fixed = fixed_boundary,
      boundary_profiled = profiled_boundary,
      widening_observed = profiled_boundary / fixed_boundary
    ),
    file.path(output_dir, "nuisance_policy_boundaries.csv")
  )
}

if (sys.nframe() == 0) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  make_figure_s1()
  run_independent_gaussian_check()
  run_nuisance_policy_check()
}
