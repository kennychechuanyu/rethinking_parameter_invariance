# Simulation 2B (Figure 5; Supplement S4): directional transfer under unequal
# task curvature.
#
# The model, lapse rate, trial count and tolerance are those of Simulation 2A;
# the tasks differ only in their stimulus levels. Task A keeps the 2A levels
# (-2 to 2) and Task B uses -0.4 to 0.4, so Task A constrains the slope more
# sharply (w_A > w_B) and the A -> B and B -> A transfer boundaries differ. Of
# four candidate designs, this pair has the profile-curvature ratio at the
# reference curve s = 2, mu = 0 closest to 4 (Supplement S4); transfer results
# play no part in the choice.
#
#   Rscript --vanilla code/R/simulation2b_directional.R               Figure 5 from saved results
#   Rscript --vanilla code/R/simulation2b_directional.R --resimulate  rerun the simulation (about a minute)

source("code/R/common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(tidyr)
})
source("code/R/psychometric_model.R")

SEED <- 20260823
S_REF <- 2 # reference slope for design selection
TARGET_CURVATURE_RATIO <- 4
POPULATION_GRID <- seq(0, 1.6, by = 0.02)
FINITE_DELTAS <- seq(0, 1.6, by = 0.1)
N_REPLICATIONS <- 2000L

CANDIDATE_DESIGNS <- list(
  moderate = X_LEVELS, # the Simulation 2A design
  informative = c(-1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6),
  near_pse = seq(-0.6, 0.6, by = 0.15),
  very_near = seq(-0.4, 0.4, by = 0.1)
)

# Design selection ----

# Effective loss of slope s: the excess loss after the location is
# re-optimized.
effective_loss <- function(s, p0, x) best_mu_pop(s, p0, x)$R

# Profile curvature: second derivative of the effective loss in the slope, by
# a central difference at the generating curve (s0, mu0).
profile_curvature <- function(x, s0 = S_REF, mu0 = 0, h = 0.01) {
  p0 <- pfun(s0, mu0, x)
  (effective_loss(s0 + h, p0, x) - 2 * effective_loss(s0, p0, x) +
     effective_loss(s0 - h, p0, x)) / h^2
}

# All ordered pairs in which Task A has the larger curvature, closest ratio
# to 4 first. The first pair is used.
rank_design_pairs <- function(design_table) {
  crossing(a = design_table$design, b = design_table$design) %>%
    filter(a != b) %>%
    left_join(design_table, by = c(a = "design")) %>%
    rename(w_a = curvature) %>%
    left_join(design_table, by = c(b = "design")) %>%
    rename(w_b = curvature) %>%
    mutate(ratio = w_a / w_b) %>%
    filter(ratio > 1) %>%
    arrange(abs(ratio - TARGET_CURVATURE_RATIO))
}

# Population analysis ----

# Population results at slope gap delta for designs x_A and x_B.
# orientation = -1 swaps the generating slopes of the two tasks.
population_2b <- function(delta, x_A, x_B, orientation = 1) {
  sA <- 2 - orientation * delta / 2
  sB <- 2 + orientation * delta / 2
  p0A <- pfun(sA, MU_A, x_A)
  p0B <- pfun(sB, MU_B, x_B)

  jobj <- function(s) effective_loss(s, p0A, x_A) + effective_loss(s, p0B, x_B)
  sJ <- optimize(jobj, interval = c(S_LO, S_HI), tol = 1e-10)$minimum

  tibble(
    delta = delta, sA = sA, sB = sB, sJ = sJ,
    R_A_to_B = effective_loss(sA, p0B, x_B), # Task A slope used in Task B
    R_B_to_A = effective_loss(sB, p0A, x_A), # Task B slope used in Task A
    R_joint_A = effective_loss(sJ, p0A, x_A),
    R_joint_B = effective_loss(sJ, p0B, x_B)
  ) %>%
    mutate(
      across(
        c(R_A_to_B, R_B_to_A, R_joint_A, R_joint_B),
        ~ pmax(.x, 0)
      ),
      R_J = pmax(R_joint_A, R_joint_B),
      R_S = pmax(R_A_to_B, R_B_to_A)
    )
}

# Gap at which a criterion's population excess loss reaches the tolerance.
boundary_2b <- function(column, x_A, x_B, orientation = 1) {
  uniroot(
    function(d) {
      population_2b(d, x_A, x_B, orientation = orientation)[[column]] -
        EPS_PRIMARY
    },
    interval = c(1e-6, 4),
    tol = 1e-10
  )$root
}

# Finite-sample simulation ----

# One simulated experiment, as run_rep() in psychometric_model.R but with a
# separate design for each task. The target-location fits use a tighter
# tolerance than Simulation 2A, which keeps optimize()'s default.
run_directional_replication <- function(delta, x_A, x_B, m = M_TRAIN) {
  sA <- 2 - delta / 2
  sB <- 2 + delta / 2
  p0A <- pfun(sA, MU_A, x_A)
  p0B <- pfun(sB, MU_B, x_B)
  kA <- rbinom(length(x_A), m, p0A)
  kB <- rbinom(length(x_B), m, p0B)

  fA <- fit_sep(kA, m, x_A)
  fB <- fit_sep(kB, m, x_B)
  fJ <- fit_joint(kA, kB, m, x_A, x_B)
  target_B_fit <- fit_target_location(fA$s, kB, m, x_B, tolerance = 1e-10)
  target_A_fit <- fit_target_location(fB$s, kA, m, x_A, tolerance = 1e-10)

  c(
    R_A_to_B = R_T(fA$s, target_B_fit$mu, p0B, x_B),
    R_B_to_A = R_T(fB$s, target_A_fit$mu, p0A, x_A),
    R_joint_A = R_T(fJ$s, fJ$muA, p0A, x_A),
    R_joint_B = R_T(fJ$s, fJ$muB, p0B, x_B),
    sA_hat = fA$s, sB_hat = fB$s, sJ_hat = fJ$s,
    muA_hat = fA$mu, muB_hat = fB$mu,
    mu_joint_A_hat = fJ$muA, mu_joint_B_hat = fJ$muB,
    mu_A_to_B_hat = target_B_fit$mu,
    mu_B_to_A_hat = target_A_fit$mu,
    conv = as.numeric(fA$conv && fB$conv && fJ$conv),
    reruns = fA$reruns + fB$reruns + fJ$reruns,
    near_bound = as.numeric(
      fA$bnd || fB$bnd || fJ$bnd ||
        target_A_fit$bnd || target_B_fit$bnd
    )
  )
}

simulate_2b <- function(x_A, x_B) {
  message(sprintf(
    "Running %s replications at each of %d slope gaps.",
    format(N_REPLICATIONS, big.mark = ","), length(FINITE_DELTAS)
  ))
  replications <- map_dfr(FINITE_DELTAS, function(d) {
    replication_matrix <- replicate_rows(
      N_REPLICATIONS, function() run_directional_replication(d, x_A, x_B)
    )
    message(sprintf("  delta %.2f done", d))
    as_tibble(replication_matrix) %>%
      mutate(delta = d, replication = seq_len(N_REPLICATIONS), .before = 1)
  })
  # As in Simulation 2A, no replication is discarded.
  stopifnot(all(replications$conv == 1), all(replications$near_bound == 0))
  replications
}

# Wilson half-width for the pass rate of a logical vector.
wilson_hw <- function(pass) wilson_interval(sum(pass), length(pass))[, "half_width"]

summarize_finite <- function(replications) {
  replications %>%
    group_by(delta) %>%
    summarise(
      replications = n(),
      A_to_B_pass_rate = mean(R_A_to_B <= EPS_PRIMARY),
      A_to_B_wilson_half_width = wilson_hw(R_A_to_B <= EPS_PRIMARY),
      B_to_A_pass_rate = mean(R_B_to_A <= EPS_PRIMARY),
      B_to_A_wilson_half_width = wilson_hw(R_B_to_A <= EPS_PRIMARY),
      bidirectional_substitution_pass_rate = mean(
        R_A_to_B <= EPS_PRIMARY & R_B_to_A <= EPS_PRIMARY
      ),
      bidirectional_substitution_wilson_half_width = wilson_hw(
        R_A_to_B <= EPS_PRIMARY & R_B_to_A <= EPS_PRIMARY
      ),
      raw_joint_pass_rate = mean(
        pmax(R_joint_A, R_joint_B) <= EPS_PRIMARY
      ),
      raw_joint_wilson_half_width = wilson_hw(
        pmax(R_joint_A, R_joint_B) <= EPS_PRIMARY
      ),
      conflation_rate = mean(
        pmax(R_joint_A, R_joint_B) <= EPS_PRIMARY &
          !(R_A_to_B <= EPS_PRIMARY & R_B_to_A <= EPS_PRIMARY)
      ),
      conflation_wilson_half_width = wilson_hw(
        pmax(R_joint_A, R_joint_B) <= EPS_PRIMARY &
          !(R_A_to_B <= EPS_PRIMARY & R_B_to_A <= EPS_PRIMARY)
      ),
      sA_bias = mean(sA_hat) - (2 - first(delta) / 2),
      sB_bias = mean(sB_hat) - (2 + first(delta) / 2),
      sA_rmse = sqrt(mean((sA_hat - (2 - first(delta) / 2))^2)),
      sB_rmse = sqrt(mean((sB_hat - (2 + first(delta) / 2))^2)),
      convergence_rate = mean(conv),
      retry_rate = mean(reruns > 0),
      near_bound_rate = mean(near_bound),
      .groups = "drop"
    )
}

# Figure 5 ----

# Top: population boundary of each criterion. Bottom: finite-sample pass rates
# for the two transfer directions. The dashed line marks delta = 0.70.
make_figure5 <- function(finite, boundary_joint, boundary_A_to_B,
                         boundary_B_to_A, example_delta = 0.70,
                         plot_max = 1.30) {
  criterion_colours <- c(
    "Joint adequacy" = orange,
    "A→B transfer" = blue,
    "B→A transfer" = "#4D4D4D"
  )

  population_limits <- tibble(
    criterion = factor(
      c("Joint adequacy", "A→B transfer", "B→A transfer"),
      levels = rev(c("Joint adequacy", "A→B transfer", "B→A transfer"))
    ),
    boundary = c(boundary_joint, boundary_A_to_B, boundary_B_to_A)
  )

  direction_rates <- finite %>%
    filter(delta <= plot_max + 1e-8) %>%
    transmute(
      delta,
      `A→B transfer` = 100 * A_to_B_pass_rate,
      `B→A transfer` = 100 * B_to_A_pass_rate
    ) %>%
    pivot_longer(-delta, names_to = "direction", values_to = "pass_rate") %>%
    mutate(
      direction = factor(
        direction,
        levels = c("A→B transfer", "B→A transfer")
      )
    )

  focus_rates <- direction_rates %>%
    filter(abs(delta - example_delta) < 1e-8)

  focus_annotations <- focus_rates %>%
    mutate(
      x = delta + if_else(direction == "A→B transfer", -0.035, 0.035),
      y = if_else(direction == "A→B transfer", 94.5, 25),
      hjust = if_else(direction == "A→B transfer", 1, 0),
      label = sprintf(
        "%s  %.0f%%",
        if_else(direction == "A→B transfer", "A→B", "B→A"),
        pass_rate
      )
    )

  baseline_annotation <- direction_rates %>%
    filter(abs(delta) < 1e-8, direction == "B→A transfer") %>%
    transmute(
      x = 0.035,
      y = 90,
      direction,
      label = sprintf("%.0f%% at δ = 0", pass_rate)
    )

  # The labels assume that delta = 0.70 lies between the B -> A and A -> B
  # boundaries.
  stopifnot(
    nrow(focus_rates) == 2L,
    nrow(baseline_annotation) == 1L,
    boundary_B_to_A < example_delta,
    example_delta < boundary_A_to_B,
    boundary_A_to_B < boundary_joint
  )

  p_population <- ggplot(population_limits, aes(y = criterion)) +
    geom_vline(
      xintercept = example_delta,
      linetype = "22",
      colour = "grey55",
      linewidth = 0.45
    ) +
    geom_segment(
      aes(x = 0, xend = boundary, yend = criterion, colour = criterion),
      linewidth = 2.1,
      lineend = "round"
    ) +
    geom_point(aes(x = boundary, colour = criterion), size = 2.6) +
    geom_text(
      aes(
        x = boundary + 0.035, label = sprintf("%.3f", boundary),
        colour = criterion
      ),
      hjust = 0,
      size = 3.0
    ) +
    annotate(
      "text",
      x = example_delta + 0.025,
      y = 3.48,
      label = "δ = 0.70",
      hjust = 0,
      colour = "grey35",
      size = 2.9
    ) +
    scale_colour_manual(values = criterion_colours) +
    scale_x_continuous(
      limits = c(0, plot_max),
      expand = expansion(add = c(0.015, 0.015))
    ) +
    labs(
      title = "Population discrepancy boundaries",
      x = NULL,
      y = NULL
    ) +
    simulation_theme() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10.5, face = "bold", margin = margin(b = 2)),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = margin(10, 9, 1, 5)
    ) +
    coord_cartesian(ylim = c(0.65, 3.58), clip = "off")

  p_finite <- ggplot(
    direction_rates,
    aes(delta, pass_rate, colour = direction, linetype = direction)
  ) +
    geom_vline(
      xintercept = example_delta,
      linetype = "22",
      colour = "grey55",
      linewidth = 0.45
    ) +
    geom_line(linewidth = 1.05) +
    geom_point(size = 1.1) +
    geom_point(
      data = focus_rates,
      shape = 21,
      fill = "white",
      size = 2.8,
      stroke = 0.9,
      show.legend = FALSE
    ) +
    geom_text(
      data = focus_annotations,
      aes(x = x, y = y, label = label, colour = direction, hjust = hjust),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 3.0,
      show.legend = FALSE
    ) +
    geom_text(
      data = baseline_annotation,
      aes(x = x, y = y, label = label, colour = direction),
      inherit.aes = FALSE,
      hjust = 0,
      size = 2.8,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = criterion_colours[c("A→B transfer", "B→A transfer")]) +
    scale_linetype_manual(values = c("A→B transfer" = "solid", "B→A transfer" = "22")) +
    scale_x_continuous(
      breaks = seq(0, 1.2, by = 0.2),
      limits = c(0, plot_max),
      expand = expansion(add = c(0.015, 0.015))
    ) +
    scale_y_continuous(
      breaks = c(0, 25, 50, 75, 100),
      limits = c(0, 100),
      expand = expansion(mult = c(0.015, 0.025))
    ) +
    labs(
      title = "Finite-sample transfer",
      x = "True slope difference, δ",
      y = "Transfer pass rate (%)",
      colour = NULL,
      linetype = NULL
    ) +
    simulation_theme() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10.5, face = "bold", margin = margin(b = 1)),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35),
      plot.margin = margin(1, 9, 3, 5)
    )

  figure <- p_population / p_finite +
    patchwork::plot_layout(heights = c(0.56, 1.12))
  save_figure(figure, file.path(figure_dir, "figure5_directional_simulation.pdf"),
              width_mm = 155, height_mm = 105)
}

# Run ----

run_simulation2b <- function(resimulate = FALSE) {
  output_dir <- file.path(results_dir, "simulation2b")
  finite_path <- file.path(output_dir, "finite_sample_by_delta.csv")
  replication_path <- file.path(output_dir, "finite_sample_replications.csv.gz")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  design_table <- tibble(
    design = names(CANDIDATE_DESIGNS),
    curvature = vapply(CANDIDATE_DESIGNS, profile_curvature, numeric(1))
  )
  pairs <- rank_design_pairs(design_table)
  x_A <- CANDIDATE_DESIGNS[[pairs$a[1]]]
  x_B <- CANDIDATE_DESIGNS[[pairs$b[1]]]
  curvature_ratio <- pairs$ratio[1]
  message(sprintf(
    "Selected %s for Task A and %s for Task B (curvature ratio %.3f).",
    pairs$a[1], pairs$b[1], curvature_ratio
  ))

  boundary_B_to_A <- boundary_2b("R_B_to_A", x_A, x_B)
  boundary_A_to_B <- boundary_2b("R_A_to_B", x_A, x_B)
  boundary_joint <- boundary_2b("R_J", x_A, x_B)
  boundary_substitution <- boundary_2b("R_S", x_A, x_B)
  message(sprintf(
    "Population boundaries: B->A %.5f, A->B %.5f, bidirectional %.5f, joint %.5f.",
    boundary_B_to_A, boundary_A_to_B, boundary_substitution, boundary_joint
  ))

  write_csv(
    design_table %>%
      mutate(
        stimulus_levels = vapply(
          CANDIDATE_DESIGNS[design],
          function(x) paste(format(x, trim = TRUE), collapse = ", "),
          character(1)
        ),
        selected_as = case_when(
          design == pairs$a[1] ~ "Task A",
          design == pairs$b[1] ~ "Task B",
          TRUE ~ "not selected"
        ),
        reference_slope = S_REF,
        reference_location = 0
      ),
    file.path(output_dir, "design_curvatures.csv")
  )
  write_csv(
    pairs %>%
      transmute(
        design_A = a,
        design_B = b,
        curvature_A = w_a,
        curvature_B = w_b,
        curvature_ratio = ratio,
        distance_from_target = abs(ratio - TARGET_CURVATURE_RATIO),
        selected = row_number() == 1L
      ),
    file.path(output_dir, "candidate_design_pairs.csv")
  )
  write_csv(
    map_dfr(POPULATION_GRID, population_2b, x_A = x_A, x_B = x_B),
    file.path(output_dir, "population_curves.csv")
  )

  # Curvature at the generating curve of each gap, to show how far the ratio
  # moves away from its reference value.
  curvature_by_delta <- map_dfr(POPULATION_GRID, function(d) {
    sA <- 2 - d / 2
    sB <- 2 + d / 2
    wA <- profile_curvature(x_A, s0 = sA, mu0 = MU_A)
    wB <- profile_curvature(x_B, s0 = sB, mu0 = MU_B)
    tibble(
      delta = d, slope_A = sA, slope_B = sB,
      curvature_A = wA, curvature_B = wB,
      curvature_ratio = wA / wB
    )
  })
  write_csv(
    curvature_by_delta,
    file.path(output_dir, "profile_curvature_by_delta.csv")
  )

  # Supplement S4: the ordering survives swapping the generating slopes.
  write_csv(
    tibble(
      slope_orientation = "s_A = 2 + delta/2; s_B = 2 - delta/2",
      boundary_B_to_A = boundary_2b("R_B_to_A", x_A, x_B, orientation = -1),
      boundary_A_to_B = boundary_2b("R_A_to_B", x_A, x_B, orientation = -1),
      boundary_joint = boundary_2b("R_J", x_A, x_B, orientation = -1)
    ),
    file.path(output_dir, "reversed_orientation_boundaries.csv")
  )

  # Exact boundaries next to the local quadratic predictions: sqrt(w_A / w_B)
  # for the ratio of the directional boundaries and 1 + w_B / w_A for the
  # joint boundary relative to A -> B.
  write_csv(
    tibble(
      boundary_B_to_A,
      boundary_A_to_B,
      boundary_substitution,
      boundary_joint,
      realized_curvature_ratio = curvature_ratio,
      design_A = pairs$a[1], design_B = pairs$b[1],
      predicted_directional_boundary_ratio = sqrt(curvature_ratio),
      observed_directional_boundary_ratio = boundary_A_to_B / boundary_B_to_A,
      predicted_joint_to_A_to_B_boundary_ratio = 1 + 1 / curvature_ratio,
      observed_joint_to_A_to_B_boundary_ratio = boundary_joint / boundary_A_to_B
    ),
    file.path(output_dir, "population_boundaries.csv")
  )

  if (resimulate) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(SEED)
    replications <- simulate_2b(x_A, x_B)
    finite <- summarize_finite(replications)
    write_csv(replications, replication_path)
    write_csv(finite, finite_path)
  } else {
    # The saved summary must agree with the saved replications.
    finite <- read_csv(finite_path, show_col_types = FALSE)
    recomputed <- summarize_finite(
      read_csv(replication_path, show_col_types = FALSE)
    )
    if (!isTRUE(all.equal(as.data.frame(finite), as.data.frame(recomputed),
                          check.attributes = FALSE, tolerance = 1e-12))) {
      stop("finite_sample_by_delta.csv does not match the saved replications.")
    }
  }

  make_figure5(finite, boundary_joint, boundary_A_to_B, boundary_B_to_A)
}

if (sys.nframe() == 0) {
  run_simulation2b(resimulate = "--resimulate" %in% commandArgs(trailingOnly = TRUE))
}
