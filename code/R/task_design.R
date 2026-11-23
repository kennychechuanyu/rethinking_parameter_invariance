# Task-design example (Figure 5; Supplement S9): tasks designed from accepted
# joint slopes, the pairwise tests among them, how the slope used to design a
# task sets the width of its acceptable region, and which slopes pass a test
# designed from them. Population calculations only; no random numbers are
# drawn.
#
#   Rscript --vanilla code/R/task_design.R

source("code/R/common.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})
source("code/R/psychometric_model.R")

output_dir <- file.path(results_dir, "task_design")

# Settings ----

CONCENTRATED_LEVELS <- seq(-0.4, 0.4, by = 0.1)  # Task B grid of Simulation 2B
TARGET_PROBABILITIES <- seq(0.25, 0.75, length.out = 9)
PP_TOLERANCE <- 5  # largest change in a predicted probability, in points
DESIGN_SLOPES <- seq(1.2, 2.9, by = 0.01)

TASK_SLOPES <- c(A = 1.70, B = 2.30, C = 2.45, D = 1.75)
ALTERNATIVE_DESIGN_SLOPE <- 1.70  # Task A's slope, the smallest of A to C

KAPPA_STEP <- 1e-3  # step for the derivative of the log curvature

# With the location fixed, the Fisher information about the slope at a
# stimulus is slope_information(z) / slope^2 with z = slope * stimulus, so
# the most informative placement corresponds to one predicted probability
# for every slope, which the locally optimal rule aims at.
slope_information <- function(z) {
  p <- LAMBDA / 2 + (1 - LAMBDA) * plogis(z)
  (1 - LAMBDA)^2 * z^2 * dlogis(z)^2 / (p * (1 - p))
}
OPTIMAL_Z <- optimize(slope_information, c(0.01, 10), maximum = TRUE,
                      tol = 1e-12)$maximum
OPTIMAL_PROBABILITY <- LAMBDA / 2 + (1 - LAMBDA) * plogis(OPTIMAL_Z)

# Design rules for the self-confirming analysis: the target probabilities
# from which Task D's stimuli are placed. The first rule is the one used in
# the example.
SELF_CONFIRMING_RULES <- list(
  "probabilities .25 to .75" = TARGET_PROBABILITIES,
  "probabilities .10 to .90" = seq(0.10, 0.90, length.out = 9),
  "probabilities .05 to .95" = seq(0.05, 0.95, length.out = 9),
  "two points .15 and .85" = c(0.15, 0.85),
  "two points, locally optimal" =
    c(1 - OPTIMAL_PROBABILITY, OPTIMAL_PROBABILITY),
  "two points .08 and .92" = c(0.08, 0.92),
  "two points .04 and .96" = c(0.04, 0.96)
)
OPTIMAL_RULE <- "two points, locally optimal"
KAPPA_CHECK_SLOPES <- c(1.20, 2.45)  # other generating slopes for kappa
SPREAD_EXTREME <- seq(0.04, 0.96, length.out = 9)  # nine points, .04 to .96

# Tasks and tests ----

# Stimulus levels at which a curve with the given slope, location zero and the
# fixed lapse rate reaches the target probabilities.
design_levels <- function(slope, probabilities = TARGET_PROBABILITIES) {
  qlogis((probabilities - LAMBDA / 2) / (1 - LAMBDA)) / slope
}

# A task's generating curve and stimulus levels. Its location is re-estimated
# when free_location is TRUE and fixed at zero otherwise.
make_task <- function(name, slope, location, x, free_location) {
  list(name = name, slope = slope, location = location, x = x,
       p0 = pfun(slope, location, x), free_location = free_location)
}

# Excess loss of a slope in a task (mean KL per trial; the generating curve is
# in the model family) and the largest change in a predicted probability, in
# percentage points. A free location is fitted by minimizing the loss.
evaluate <- function(task, slope) {
  location <- if (task$free_location) {
    best_mu_pop(slope, task$p0, task$x)$mu
  } else {
    0
  }
  p <- pfun(slope, location, task$x)
  c(R = mean(kl_bern(task$p0, p)), pp = 100 * max(abs(p - task$p0)))
}

# Slope that minimizes the summed excess loss of the tasks.
joint_slope <- function(tasks) {
  optimize(
    function(slope) sum(map_dbl(tasks, ~ evaluate(.x, slope)[["R"]])),
    interval = c(S_LO, S_HI), tol = 1e-10
  )$minimum
}

# Curvature of a task's excess loss at its own slope, by central differences.
curvature <- function(task, h = 1e-3) {
  loss <- function(slope) evaluate(task, slope)[["R"]]
  (loss(task$slope + h) - 2 * loss(task$slope) + loss(task$slope - h)) / h^2
}

# Joint fit and both transfer directions for a pair of tasks. Joint adequacy
# uses the calibrated tolerance, transfer the common one. The pair's own
# calibration factor, w1 * w2 / (w1 + w2)^2, checks the equal-curvature
# factor of one quarter.
pair_test <- function(first, second) {
  slope <- joint_slope(list(first, second))
  joint <- pmax(evaluate(first, slope), evaluate(second, slope))
  forward <- evaluate(second, first$slope)
  backward <- evaluate(first, second$slope)
  w_first <- curvature(first)
  w_second <- curvature(second)
  tibble(
    pair = paste(first$name, second$name, sep = "-"),
    joint_slope = slope,
    R_joint = joint[["R"]],
    R_first_to_second = forward[["R"]],
    R_second_to_first = backward[["R"]],
    pp_joint = joint[["pp"]],
    pp_first_to_second = forward[["pp"]],
    pp_second_to_first = backward[["pp"]],
    calibrated_joint = R_joint <= EPS_JOINT_CALIBRATED,
    bidirectional_transfer =
      max(R_first_to_second, R_second_to_first) <= EPS_PRIMARY,
    pair_calibration_factor = w_first * w_second / (w_first + w_second)^2,
    pair_calibrated_joint = R_joint <= EPS_PRIMARY * pair_calibration_factor
  )
}

# Acceptable region of a task with fixed location: the slopes whose excess
# loss (criterion "kl") or largest probability change (criterion "pp") stays
# within tolerance.
acceptable_region <- function(task, criterion = c("kl", "pp")) {
  criterion <- match.arg(criterion)
  excess <- if (criterion == "kl") {
    function(slope) evaluate(task, slope)[["R"]] - EPS_PRIMARY
  } else {
    function(slope) evaluate(task, slope)[["pp"]] - PP_TOLERANCE
  }
  c(
    lower = uniroot(excess, c(S_LO, task$slope), tol = 1e-10)$root,
    upper = uniroot(excess, c(task$slope, S_HI), tol = 1e-10)$root
  )
}

# Task D designed with a given slope: its regions under both criteria and the
# range of response probabilities it realizes.
task_D_summary <- function(design_slope) {
  task <- make_task("D", TASK_SLOPES[["D"]], 0, design_levels(design_slope),
                    FALSE)
  kl <- acceptable_region(task, "kl")
  pp <- acceptable_region(task, "pp")
  tibble(
    design_slope = design_slope,
    kl_lower = kl[["lower"]], kl_upper = kl[["upper"]],
    kl_width = kl[["upper"]] - kl[["lower"]],
    pp_lower = pp[["lower"]], pp_upper = pp[["upper"]],
    pp_width = pp[["upper"]] - pp[["lower"]],
    p_min = min(task$p0), p_max = max(task$p0)
  )
}

# Self-confirming slopes ----

# Task D with the given generating slope, its stimuli placed from the design
# slope under a design rule.
designed_task <- function(slope, design,
                          probabilities = TARGET_PROBABILITIES) {
  make_task("D", slope, 0, design_levels(design, probabilities), FALSE)
}

# Excess loss of a slope in a task designed from that same slope. With the
# location fixed, the designed curve reaches the target probabilities, so
# this is the mean KL between the generating and target probabilities.
self_test_loss <- function(value, slope,
                           probabilities = TARGET_PROBABILITIES) {
  evaluate(designed_task(slope, value, probabilities), value)[["R"]]
}

# Slopes that pass the test designed from them, under the tolerance
# (criterion "kl") or the five-point limit on any probability change ("pp").
self_confirming_range <- function(slope, probabilities,
                                  criterion = c("kl", "pp")) {
  criterion <- match.arg(criterion)
  excess <- function(value) {
    task <- designed_task(slope, value, probabilities)
    if (criterion == "kl") {
      evaluate(task, value)[["R"]] - EPS_PRIMARY
    } else {
      evaluate(task, value)[["pp"]] - PP_TOLERANCE
    }
  }
  c(lower = uniroot(excess, c(S_LO, slope), tol = 1e-10)$root,
    upper = uniroot(excess, c(slope, S_HI), tol = 1e-10)$root)
}

# Derivative of the log curvature at the generating slope with respect to
# the slope used for the design.
design_kappa <- function(slope, probabilities) {
  log_curvature <- function(design) {
    log(curvature(designed_task(slope, design, probabilities)))
  }
  (log_curvature(slope + KAPPA_STEP) - log_curvature(slope - KAPPA_STEP)) /
    (2 * KAPPA_STEP)
}

# For each design rule: kappa, the acceptable region of a design from the
# generating slope, the self-confirming range, and the shift of that range
# relative to the region, exact and to first order (-kappa r^2 / 2).
self_confirming_table <- function(slope = TASK_SLOPES[["D"]]) {
  bind_rows(imap(SELF_CONFIRMING_RULES, function(probabilities, rule) {
    own_task <- designed_task(slope, slope, probabilities)
    own <- acceptable_region(own_task, "kl")
    self <- self_confirming_range(slope, probabilities)
    kappa <- design_kappa(slope, probabilities)
    radius <- sqrt(2 * EPS_PRIMARY / curvature(own_task))
    tibble(
      rule = rule, levels = length(probabilities), kappa = kappa,
      own_lower = own[["lower"]], own_upper = own[["upper"]],
      self_lower = self[["lower"]], self_upper = self[["upper"]],
      shift = mean(self) - mean(own),
      first_order_shift = -kappa * radius^2 / 2
    )
  }))
}

# Checks quoted in the text: the identity behind the self-test, the range of
# Task D slopes over which the accepted value still confirms itself, the
# losses, in units of the tolerance, if Task D's slope were 1.70, kappa under
# the locally optimal rule at other generating slopes, the ranges under the
# five-point criterion, and Task D designed from the B-C joint slope with the
# locally optimal rule.
self_confirming_checks <- function(design_BC) {
  s <- TASK_SLOPES
  alt <- ALTERNATIVE_DESIGN_SLOPE
  task_D <- designed_task(s[["D"]], design_BC)
  confirms <- function(slope) {
    self_test_loss(design_BC, slope) - EPS_PRIMARY
  }
  in_low_D <- function(design, value) {
    evaluate(designed_task(alt, design), value)[["R"]] / EPS_PRIMARY
  }
  own_task <- designed_task(s[["D"]], s[["D"]])
  optimal_D <- designed_task(s[["D"]], design_BC,
                             SELF_CONFIRMING_RULES[[OPTIMAL_RULE]])
  optimal_region <- unname(acceptable_region(optimal_D, "kl"))
  pp_own <- unname(acceptable_region(own_task, "pp"))
  pp_self <- unname(self_confirming_range(s[["D"]], TARGET_PROBABILITIES,
                                         "pp"))
  tibble(
    quantity = c(
      "self_test_loss_BC_joint_in_D",
      "kl_D_to_target_probabilities",
      "lowest_D_slope_BC_joint_confirms_itself",
      "highest_D_slope_BC_joint_confirms_itself",
      "BC_joint_in_low_D_designed_from_low_slope",
      "C_in_low_D_designed_from_BC_joint",
      "C_in_low_D_designed_from_low_slope",
      "optimal_predicted_probability",
      paste0("kappa_optimal_rule_slope_", KAPPA_CHECK_SLOPES),
      "kappa_nine_points_04_to_96",
      "pp_own_lower", "pp_own_upper", "pp_self_lower", "pp_self_upper",
      "BC_joint_in_D_optimal_rule", "C_in_D_optimal_rule",
      "D_optimal_rule_region_lower", "D_optimal_rule_region_upper"
    ),
    value = c(
      evaluate(task_D, design_BC)[["R"]],
      mean(kl_bern(task_D$p0, TARGET_PROBABILITIES)),
      uniroot(confirms, c(0.8, design_BC), tol = 1e-10)$root,
      uniroot(confirms, c(design_BC, S_HI), tol = 1e-10)$root,
      in_low_D(alt, design_BC),
      in_low_D(design_BC, s[["C"]]),
      in_low_D(alt, s[["C"]]),
      OPTIMAL_PROBABILITY,
      map_dbl(KAPPA_CHECK_SLOPES,
              ~ design_kappa(.x, SELF_CONFIRMING_RULES[[OPTIMAL_RULE]])),
      design_kappa(s[["D"]], SPREAD_EXTREME),
      pp_own, pp_self,
      evaluate(optimal_D, design_BC)[["R"]] / EPS_PRIMARY,
      evaluate(optimal_D, s[["C"]])[["R"]] / EPS_PRIMARY,
      optimal_region
    ),
    unit = c("nats per trial", "nats per trial", "slope", "slope",
             "tolerance", "tolerance", "tolerance", "probability",
             rep("per unit slope", length(KAPPA_CHECK_SLOPES) + 1),
             rep("slope", 4), "tolerance", "tolerance", "slope", "slope")
  )
}

# Analysis ----

run_task_design <- function() {
  s <- TASK_SLOPES
  task_A <- make_task("A", s[["A"]], MU_A, CONCENTRATED_LEVELS, TRUE)
  task_B <- make_task("B", s[["B"]], MU_B, CONCENTRATED_LEVELS, TRUE)
  slope_AB <- joint_slope(list(task_A, task_B))
  task_C <- make_task("C", s[["C"]], 0, design_levels(slope_AB), FALSE)
  slope_BC <- joint_slope(list(task_B, task_C))
  task_D <- make_task("D", s[["D"]], 0, design_levels(slope_BC), FALSE)

  # Alternatives: Task D designed from Task A's slope, the smallest optimum so
  # far, and Tasks C and D on the broad grid of Simulation 2A.
  task_D_alternative <- make_task(
    "D_alternative", s[["D"]], 0, design_levels(ALTERNATIVE_DESIGN_SLOPE), FALSE
  )
  broad_C <- make_task("C_broad", s[["C"]], 0, X_LEVELS, FALSE)
  broad_D <- make_task("D_broad", s[["D"]], 0, X_LEVELS, FALSE)

  tasks <- list(task_A, task_B, task_C, task_D)
  tests <- bind_rows(
    map(combn(4, 2, simplify = FALSE),
        ~ pair_test(tasks[[.x[1]]], tasks[[.x[2]]])),
    pair_test(task_C, task_D_alternative),
    pair_test(broad_C, broad_D)
  )

  slope_ABCD <- joint_slope(tasks)
  four_task_joint <- tibble(
    task = map_chr(tasks, "name"),
    joint_slope = slope_ABCD,
    R_joint = map_dbl(tasks, ~ evaluate(.x, slope_ABCD)[["R"]]),
    pp_joint = map_dbl(tasks, ~ evaluate(.x, slope_ABCD)[["pp"]])
  )

  designs <- bind_rows(map(
    list(task_A, task_B, task_C, task_D, task_D_alternative),
    ~ tibble(task = .x$name, level = seq_along(.x$x), x = .x$x, p0 = .x$p0)
  ))
  design_slopes <- tibble(
    task = c("C", "D", "D_alternative"),
    design_slope = c(slope_AB, slope_BC, ALTERNATIVE_DESIGN_SLOPE)
  )
  D_regions <- bind_rows(map(c(DESIGN_SLOPES, slope_BC), task_D_summary)) %>%
    arrange(design_slope)

  # Task C under the design it received (the A-B joint slope, below its own
  # slope) and under a design from its own slope.
  task_C_regions <- bind_rows(map(c(slope_AB, s[["C"]]), function(design) {
    task <- make_task("C", s[["C"]], 0, design_levels(design), FALSE)
    kl <- acceptable_region(task, "kl")
    tibble(design_slope = design, kl_lower = kl[["lower"]],
           kl_upper = kl[["upper"]], kl_width = kl[["upper"]] - kl[["lower"]],
           p_min = min(task$p0), p_max = max(task$p0))
  }))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv(tests, file.path(output_dir, "pairwise_tests.csv"))
  write_csv(four_task_joint, file.path(output_dir, "four_task_joint.csv"))
  write_csv(left_join(designs, design_slopes, by = "task"),
            file.path(output_dir, "designs.csv"))
  write_csv(D_regions, file.path(output_dir, "task_D_regions.csv"))
  write_csv(task_C_regions, file.path(output_dir, "task_C_regions.csv"))

  self_confirming <- self_confirming_table()
  write_csv(self_confirming, file.path(output_dir, "self_confirming.csv"))
  write_csv(self_confirming_checks(slope_BC),
            file.path(output_dir, "self_confirming_checks.csv"))

  make_figure5(D_regions, tests, unlist(self_confirming[1, c("self_lower",
                                                           "self_upper")]))
}

# Figure 5 ----

# (A) The inverse stimulus-placement rule. (B) Population acceptable slopes
# as the design slope varies, with the two selected designs and the
# diagonal on which the design slope is also the slope evaluated.
# Rounds an interval inward to two decimals, so that every printed value lies
# inside it.
inward_round <- function(interval) {
  c(ceiling(interval[1] * 100), floor(interval[2] * 100)) / 100
}

make_figure5 <- function(regions, pairs, self_confirming) {
  design_BC <- pairs$joint_slope[pairs$pair == "B-C"]
  slope_C <- TASK_SLOPES[["C"]]
  selected <- tibble(
    task = c("D_alternative", "D"),
    design_slope = c(ALTERNATIVE_DESIGN_SLOPE, design_BC),
    pair = c("C-D_alternative", "C-D")
  ) %>%
    left_join(regions, by = "design_slope") %>%
    mutate(loss = pairs$R_first_to_second[match(pair, pairs$pair)] /
             EPS_PRIMARY,
           outcome = if_else(loss <= 1, "passes", "fails"),
           label = sprintf("%.2f \u03b5: %s", loss, outcome),
           label_y = slope_C + c(0.24, -0.24))

  colours <- c(D = orange, D_alternative = blue)
  line_types <- c(D = "solid", D_alternative = "44")
  design_labels <- c(D = sprintf("Design Slope %.2f", design_BC),
                     D_alternative = sprintf("Design Slope %.2f",
                                             ALTERNATIVE_DESIGN_SLOPE))

  # Invert the assumed response curve: the same targets yield different inputs.
  targets <- seq(min(TARGET_PROBABILITIES), max(TARGET_PROBABILITIES),
                 length.out = 201)
  placements <- bind_rows(lapply(seq_len(nrow(selected)), function(i) {
    tibble(
      task = selected$task[i], target = targets,
      stimulus = qlogis((targets - LAMBDA / 2) / (1 - LAMBDA)) /
        selected$design_slope[i]
    )
  }))
  endpoints <- placements[placements$target %in% range(targets), ]

  placement <- ggplot(placements, aes(target, stimulus, colour = task,
                                      linetype = task)) +
    geom_hline(yintercept = 0, colour = "grey88", linewidth = 0.35) +
    geom_line(linewidth = 0.7) +
    geom_point(data = endpoints, size = 1.6, show.legend = FALSE) +
    scale_colour_manual(values = colours, labels = design_labels,
                        breaks = c("D", "D_alternative")) +
    scale_linetype_manual(values = line_types, labels = design_labels,
                          breaks = c("D", "D_alternative")) +
    scale_x_continuous(breaks = c(0.25, 0.5, 0.75),
                       labels = c(".25", ".50", ".75"),
                       expand = expansion(mult = 0)) +
    scale_y_continuous(breaks = c(-0.5, 0, 0.5),
                       expand = expansion(mult = 0)) +
    coord_cartesian(xlim = c(0.23, 0.77), ylim = c(-0.76, 0.88)) +
    labs(title = "A  Stimulus placement", x = "Intended probability",
         y = "Selected stimulus value") +
    paper_theme(9) +
    theme(
      plot.title = element_text(size = 9.5, face = "bold"),
      legend.position = "inside",
      legend.position.inside = c(0.03, 0.98),
      legend.justification = c(0, 1),
      legend.text = element_text(size = 8.5),
      legend.key.width = unit(0.75, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.background = element_blank(),
      plot.margin = margin(6, 8, 5, 4)
    )

  # The generating D curve stays fixed as its stimulus placement changes.
  adequacy <- ggplot(regions, aes(design_slope)) +
    geom_ribbon(aes(ymin = kl_lower, ymax = kl_upper), fill = "grey92") +
    geom_line(aes(y = kl_lower), colour = "grey60", linewidth = 0.5) +
    geom_line(aes(y = kl_upper), colour = "grey60", linewidth = 0.5) +
    geom_hline(yintercept = slope_C, colour = "grey40", linetype = "22",
               linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, colour = "grey35",
                linetype = "13", linewidth = 0.45) +
    annotate("point", x = self_confirming, y = self_confirming, shape = 21,
             size = 1.8, stroke = 0.6, fill = "white", colour = "grey25") +
    annotate("text", x = self_confirming + c(0.05, -0.05),
             y = self_confirming + c(-0.06, 0.06),
             label = sprintf("%.2f", inward_round(self_confirming)),
             hjust = c(0, 1),
             size = 2.9, colour = "grey25") +
    geom_segment(data = selected,
                 aes(xend = design_slope, y = kl_lower, yend = kl_upper,
                     colour = task), linewidth = 1.1) +
    geom_text(data = selected,
              aes(x = design_slope + 0.035, y = label_y, label = label,
                  colour = task), hjust = 0, size = 3.0) +
    annotate("text", x = 1.2, y = slope_C + 0.05,
             label = sprintf("C's slope: %.2f", slope_C),
             hjust = 0, vjust = 0, size = 3.0, colour = "grey30") +
    annotate("text", x = 2.05, y = 1.55, label = "Acceptable in D",
             size = 2.9, colour = "grey40") +
    scale_colour_manual(values = colours, guide = "none") +
    scale_x_continuous(
      breaks = c(1.2, ALTERNATIVE_DESIGN_SLOPE, 2, design_BC, 2.8),
      labels = c("1.2", sprintf("%.2f\nEarlier A", ALTERNATIVE_DESIGN_SLOPE),
                 "2.0", sprintf("%.2f\nB\u2013C joint", design_BC), "2.8"),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(breaks = c(1, 1.5, 2, 2.5),
                       expand = expansion(mult = 0)) +
    coord_cartesian(xlim = c(1.15, 2.95), ylim = c(0.65, 2.95)) +
    labs(title = "B  Acceptable slopes",
         subtitle = sprintf("D's generating slope stays at %.2f",
                            TASK_SLOPES[["D"]]),
         x = "Slope used to design D", y = "Slope evaluated in D") +
    paper_theme(9) +
    theme(
      plot.title = element_text(size = 9.5, face = "bold"),
      plot.subtitle = element_text(size = 8.5, colour = "grey35",
                                   margin = margin(b = 4)),
      plot.margin = margin(6, 8, 5, 5)
    )

  figure <- (placement | adequacy) + plot_layout(widths = c(1, 1.65))

  save_figure(figure, file.path(figure_dir, "figure5_task_design.pdf"),
              width_mm = 170, height_mm = 86)
  invisible(figure)
}

if (sys.nframe() == 0) {
  run_task_design()
}
