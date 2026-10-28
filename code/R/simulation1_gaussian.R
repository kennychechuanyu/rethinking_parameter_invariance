# Simulation 1 (Figure 3; Supplement S2): Gaussian means with matched
# precision.
#
# Each task mean is estimated from 1,600 observations. The joint estimate
# averages a fixed subset of 800 observations from each task, so all three
# estimators have standard error 0.025. Because the subsets are shared, the
# joint estimate has correlation .5 with each task estimate.
#
# The paper's acceptable radius r and standard error tau are a and s here.

source("code/R/common.R")

run_simulation1 <- function() {
  # Settings ----
  seed <- 20240115
  sigma <- 1               # observation SD, known
  n_total <- 1600          # observations per task
  n_joint_each <- 800      # observations per task in the joint estimate
  eps <- 0.02              # predictive tolerance, nats per observation
  a <- sqrt(2 * eps)       # acceptable radius, 0.20
  eps_joint_cal <- eps / 4 # calibrated joint tolerance
  a_joint_cal <- sqrt(2 * eps_joint_cal)
  delta_grid <- seq(0, 0.60, by = 0.01)
  n_rep <- 20000           # replications per delta

  output_results_dir <- file.path(results_dir, "simulation1")
  dir.create(output_results_dir, recursive = TRUE, showWarnings = FALSE)
  data_path <- function(filename) file.path(output_results_dir, filename)

  # Both estimators have variance 1/1600.
  var_S <- sigma^2 / n_total
  var_J <- (1 / 4) * (sigma^2 / n_joint_each + sigma^2 / n_joint_each)
  s2 <- var_S
  s <- sqrt(s2)

  # Closed-form results ----

  # Excess loss (mean Gaussian NLL, unit variance) of predicting with c when the
  # task optimum is theta.
  risk <- function(c, theta) 0.5 * (c - theta)^2

  # Expected excess loss of the two procedures, including sampling error.
  E_R_substitution <- function(delta) 0.5 * (delta^2 + s2)
  E_R_joint <- function(delta) 0.5 * (delta^2 / 4 + s2)

  # Pass probabilities (Supplement S2). The joint estimate is centred at zero
  # and passes when it lies within a - delta/2 of zero.
  P_joint_success <- function(delta) {
    ifelse(delta <= 2 * a, 2 * pnorm((a - delta / 2) / s) - 1, 0)
  }
  # Under eps/4 the joint radius is a/2, which moves the joint boundary to the
  # transfer boundary delta = a.
  P_joint_cal_success <- function(delta) {
    ifelse(
      delta <= 2 * a_joint_cal,
      2 * pnorm((a_joint_cal - delta / 2) / s) - 1,
      0
    )
  }
  # One transfer direction, equal to Phi((r - delta)/tau) - Phi((-r - delta)/tau)
  # in Supplement S2. The two directions use independent task data.
  p_S <- function(delta) pnorm((delta + a) / s) - pnorm((delta - a) / s)
  P_substitution_success <- function(delta) p_S(delta)^2

  # Monte Carlo simulation ----

  # The four half-sample means are independent. Each task mean averages the two
  # halves of its task; the joint mean averages the first half of each task,
  # which is also the minimizer of the summed joint NLL.
  run_mc <- function(delta_grid, n_rep, seed) {
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(seed)
    out <- vector("list", length(delta_grid))
    for (i in seq_along(delta_grid)) {
      delta <- delta_grid[i]
      theta_A <- -delta / 2
      theta_B <- delta / 2

      subset_sd <- sigma / sqrt(n_joint_each)
      mA_used <- rnorm(n_rep, mean = theta_A, sd = subset_sd)
      mA_rest <- rnorm(n_rep, mean = theta_A, sd = subset_sd)
      mB_used <- rnorm(n_rep, mean = theta_B, sd = subset_sd)
      mB_rest <- rnorm(n_rep, mean = theta_B, sd = subset_sd)

      mA_S <- 0.5 * (mA_used + mA_rest) # Task A mean, 1,600 observations
      mB_S <- 0.5 * (mB_used + mB_rest) # Task B mean
      mJ <- 0.5 * (mA_used + mB_used)   # joint mean, 800 + 800

      rB_of_A <- risk(mA_S, theta_B) # A -> B
      rA_of_B <- risk(mB_S, theta_A) # B -> A
      rA_of_J <- risk(mJ, theta_A)
      rB_of_J <- risk(mJ, theta_B)

      substitution_success <- (rB_of_A <= eps) & (rA_of_B <= eps)
      joint_success <- (rA_of_J <= eps) & (rB_of_J <= eps)
      joint_cal_success <-
        (rA_of_J <= eps_joint_cal) & (rB_of_J <= eps_joint_cal)
      conflation <- joint_success & !substitution_success

      substitution_risk_rep <- 0.5 * (rB_of_A + rA_of_B)
      joint_risk_rep <- 0.5 * (rA_of_J + rB_of_J)

      out[[i]] <- data.frame(
        delta = delta,
        q = delta / a,
        bias_A_S = mean(mA_S) - theta_A,
        bias_B_S = mean(mB_S) - theta_B,
        bias_J = mean(mJ) - (theta_A + theta_B) / 2,
        sd_A_S = sd(mA_S),
        sd_B_S = sd(mB_S),
        sd_J = sd(mJ),
        cor_J_A = cor(mJ, mA_S),
        cor_J_B = cor(mJ, mB_S),
        substitution_risk_mc = mean(substitution_risk_rep),
        joint_risk_mc = mean(joint_risk_rep),
        substitution_risk_se = sd(substitution_risk_rep) / sqrt(n_rep),
        joint_risk_se = sd(joint_risk_rep) / sqrt(n_rep),
        substitution_succ_k = sum(substitution_success),
        joint_succ_k = sum(joint_success),
        joint_cal_succ_k = sum(joint_cal_success),
        conflation_k = sum(conflation),
        both_pass_k = sum(joint_success & substitution_success),
        joint_only_k = sum(joint_success & !substitution_success),
        substitution_only_k = sum(!joint_success & substitution_success),
        neither_pass_k = sum(!joint_success & !substitution_success),
        cal_joint_and_substitution_k = sum(joint_cal_success & substitution_success),
        cal_joint_only_k = sum(joint_cal_success & !substitution_success),
        substitution_only_vs_cal_joint_k = sum(!joint_cal_success & substitution_success),
        neither_cal_joint_nor_substitution_k = sum(!joint_cal_success & !substitution_success),
        cal_joint_outside_raw_joint_k = sum(joint_cal_success & !joint_success),
        replications = n_rep
      )
    }
    do.call(rbind, out)
  }

  mc <- run_mc(delta_grid, n_rep, seed)

  mc$substitution_risk_exact <- E_R_substitution(mc$delta)
  mc$joint_risk_exact <- E_R_joint(mc$delta)
  mc$P_substitution_exact <- P_substitution_success(mc$delta)
  mc$P_joint_exact <- P_joint_success(mc$delta)
  mc$P_joint_cal_exact <- P_joint_cal_success(mc$delta)

  # Simulation 1 uses the exact normal quantile for its Wilson intervals.
  ci_substitution <- wilson_interval(mc$substitution_succ_k, mc$replications,
                                     z = qnorm(0.975))
  ci_joint <- wilson_interval(mc$joint_succ_k, mc$replications,
                              z = qnorm(0.975))
  ci_cal <- wilson_interval(mc$joint_cal_succ_k, mc$replications,
                            z = qnorm(0.975))
  mc$substitution_p <- mc$substitution_succ_k / mc$replications
  mc$substitution_p_lo <- ci_substitution[, "lower"]
  mc$substitution_p_hi <- ci_substitution[, "upper"]
  mc$joint_p <- mc$joint_succ_k / mc$replications
  mc$joint_p_lo <- ci_joint[, "lower"]
  mc$joint_p_hi <- ci_joint[, "upper"]
  mc$joint_cal_p <- mc$joint_cal_succ_k / mc$replications
  mc$joint_cal_p_lo <- ci_cal[, "lower"]
  mc$joint_cal_p_hi <- ci_cal[, "upper"]

  # Checks on the closed-form results ----
  `%approx%` <- function(x, y) abs(x - y) < 5e-5
  scenario_report <- function(delta) {
    q <- delta / a
    data.frame(
      q = q,
      delta = delta,
      substitution_risk_exact = E_R_substitution(delta),
      joint_risk_exact = E_R_joint(delta),
      P_joint = P_joint_success(delta),
      P_joint_cal = P_joint_cal_success(delta),
      P_substitution = P_substitution_success(delta)
    )
  }
  sc <- do.call(rbind, lapply(c(0.10, 0.30, 0.50), scenario_report))

  checks <- c(
    substitution_q0_5 = sc$substitution_risk_exact[1] %approx% 0.00531,
    joint_q0_5 = sc$joint_risk_exact[1] %approx% 0.00156,
    substitution_q1_5 = sc$substitution_risk_exact[2] %approx% 0.04531,
    joint_q1_5 = sc$joint_risk_exact[2] %approx% 0.01156,
    substitution_q2_5 = sc$substitution_risk_exact[3] %approx% 0.12531,
    joint_q2_5 = sc$joint_risk_exact[3] %approx% 0.03156,
    Pjoint_q1_5 = abs(sc$P_joint[2] - 0.9545) < 2e-3,
    Pcal_q0_5 = abs(sc$P_joint_cal[1] - 0.9545) < 2e-3,
    cal_subset = all(mc$cal_joint_outside_raw_joint_k == 0)
  )
  stopifnot(all(checks))

  # Diagnostics ----
  emp_sd <- mean(c(mc$sd_A_S, mc$sd_B_S, mc$sd_J))
  emp_bias <- mean(abs(c(mc$bias_A_S, mc$bias_B_S, mc$bias_J)))
  # Fourfold loss ratio, after removing the common sampling floor s^2/2.
  row10 <- mc[which(abs(mc$delta - 0.10) < 1e-9), ]
  row30 <- mc[which(abs(mc$delta - 0.30) < 1e-9), ]
  disc_substitution_emp <- row30$substitution_risk_mc - s2 / 2
  disc_joint_emp <- row30$joint_risk_mc - s2 / 2
  ratio_emp <- disc_substitution_emp / disc_joint_emp

  audit <- data.frame(
    quantity = c(
      "SE_substitution (theory)", "SE_joint (theory)", "common SE s (theory)",
      "empirical mean |bias| (all 3 estimators)",
      "empirical mean SD (all 3 estimators)",
      "target SD = s",
      "mean corr(joint, Task A estimate)",
      "mean corr(joint, Task B estimate)",
      "common delta=0 risk floor (theory s^2/2)",
      "delta=0 substitution risk (MC)", "delta=0 joint risk (MC)",
      "discrepancy ratio substitution/joint (theory)",
      "discrepancy ratio substitution/joint (empirical, delta=.30)",
      "substitution population boundary q",
      "joint population boundary q",
      "calibrated joint population boundary q",
      "eps (original rule, nats/obs)",
      "eps (calibrated joint rule, nats/obs)",
      "a = sqrt(2 eps)",
      "a_cal = sqrt(2 eps_cal)",
      "scenario q=0.5 (delta=.10) substitution risk",
      "scenario q=0.5 (delta=.10) joint risk",
      "scenario q=0.5 (delta=.10) P(calibrated joint)",
      "scenario q=0.5 (delta=.10) calibrated-joint passes (MC)",
      "scenario q=1.5 (delta=.30) substitution risk",
      "scenario q=1.5 (delta=.30) joint risk",
      "scenario q=1.5 (delta=.30) P(joint success)",
      "scenario q=1.5 (delta=.30) P(substitution bidirectional)",
      "scenario q=1.5 (delta=.30) P(calibrated joint)",
      "scenario q=1.5 (delta=.30) calibrated-joint passes (MC)",
      "scenario q=2.5 (delta=.50) substitution risk",
      "scenario q=2.5 (delta=.50) joint risk",
      "calibrated passes outside original joint (all conditions)",
      "MC reps per condition", "n delta conditions", "seed",
      "discarded or replaced replications"
    ),
    value = c(
      sqrt(var_S), sqrt(var_J), s,
      emp_bias, emp_sd, s,
      mean(mc$cor_J_A), mean(mc$cor_J_B),
      s2 / 2,
      mc$substitution_risk_mc[1], mc$joint_risk_mc[1],
      4, ratio_emp,
      1, 2, 1,
      eps, eps_joint_cal, a, a_joint_cal,
      sc$substitution_risk_exact[1], sc$joint_risk_exact[1],
      sc$P_joint_cal[1], row10$joint_cal_succ_k,
      sc$substitution_risk_exact[2], sc$joint_risk_exact[2], sc$P_joint[2], sc$P_substitution[2],
      sc$P_joint_cal[2], row30$joint_cal_succ_k,
      sc$substitution_risk_exact[3], sc$joint_risk_exact[3],
      sum(mc$cal_joint_outside_raw_joint_k),
      n_rep, length(delta_grid), seed, 0
    ),
    note = c(
      "sigma^2/n_total", "(1/4)(sigma^2/800 + sigma^2/800)", "matched",
      "expect ~0 (all unbiased)",
      "expect ~0.025", "0.025",
      "theory 0.5 from shared subset", "theory 0.5 from shared subset",
      "0.0003125", "MC", "MC",
      "delta^2/2 vs delta^2/8", "MC risk minus common floor",
      "delta = a", "delta = 2a", "delta = a",
      "joint and substitution", "one quarter of the original tolerance",
      "original acceptable absolute mean discrepancy",
      "calibrated acceptable absolute mean discrepancy",
      "both succeed", "both succeed",
      "exact", sprintf("out of %s", format(n_rep, big.mark = ",")),
      "joint succeeds, substitution fails", "joint succeeds", "~0.9545", "~1e-9",
      "exactly zero beyond delta=a",
      sprintf("out of %s", format(n_rep, big.mark = ",")),
      "both fail", "both fail",
      "must be zero because eps/4 is stricter",
      "", "", "", "none"
    ),
    stringsAsFactors = FALSE
  )
  # Six significant digits, except for the seed.
  is_seed <- audit$quantity == "seed"
  audit$value[!is_seed] <- signif(audit$value[!is_seed], 6)
  write.csv(audit, data_path("simulation1_diagnostics.csv"), row.names = FALSE)
  write.csv(mc, data_path("simulation1_curves.csv"), row.names = FALSE)

  # Figure 3 ----

  # Panel A shows the estimators at delta = 0.30.
  delta_A <- 0.30
  thA <- -delta_A / 2
  thB <- delta_A / 2

  col_shade <- "grey92"
  col_accept <- "#E7F1E2"
  method_colours <- c(
    "Joint value" = orange,
    "Transferred values" = blue
  )
  method_shapes <- c(
    "Joint value" = 16,
    "Transferred values" = 17
  )
  method_linetypes <- c(
    "Joint value" = "solid",
    "Transferred values" = "22"
  )
  method_labels <- c(
    "Joint value" = "Joint estimate",
    "Transferred values" = "Transferred estimate"
  )

  panelA_data <- data.frame(
    task = c(
      "Predict Task A", "Predict Task A",
      "Predict Task B", "Predict Task B"
    ),
    y = c(2, 2, 1, 1),
    method = factor(
      rep(c("Joint value", "Transferred values"), 2),
      levels = c("Joint value", "Transferred values")
    ),
    mean = c(0, thB, 0, thA),
    pass_lower = c(thA - a, thA - a, thB - a, thB - a),
    pass_upper = c(thA + a, thA + a, thB + a, thB + a)
  )
  panelA_data$lower <- panelA_data$mean - 1.96 * s
  panelA_data$upper <- panelA_data$mean + 1.96 * s
  pA_bands <- panelA_data[!duplicated(panelA_data$task), ]
  pA <- ggplot(panelA_data, aes(y = y, colour = method, shape = method)) +
    geom_rect(
      data = pA_bands,
      aes(
        xmin = pass_lower, xmax = pass_upper,
        ymin = y - 0.24, ymax = y + 0.24
      ),
      inherit.aes = FALSE,
      fill = col_accept, colour = "grey70", linewidth = 0.25
    ) +
    geom_segment(aes(x = lower, xend = upper, yend = y), linewidth = 1.15) +
    geom_point(aes(x = mean), size = 3.0) +
    scale_colour_manual(values = method_colours, labels = method_labels) +
    scale_shape_manual(values = method_shapes, labels = method_labels) +
    scale_x_continuous(
      breaks = c(-0.30, -0.15, 0, 0.15, 0.30),
      limits = c(-0.40, 0.40), expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = c(1, 2),
      labels = c("Task B target", "Task A target"),
      limits = c(0.55, 2.45), expand = c(0, 0)
    ) +
    labs(
      title = "Equal precision, different estimator centers",
      x = "Mean used for prediction", y = NULL,
      colour = NULL, shape = NULL
    ) +
    simulation_theme() +
    theme(
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.35),
      legend.position = "top",
      legend.justification = "left",
      legend.box.just = "left",
      legend.margin = margin(t = 0, b = 1)
    ) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2.5)))

  # Panel B: exact pass probabilities, with Monte Carlo estimates at three gaps.
  exact_curves <- rbind(
    data.frame(
      delta = mc$delta, p = 100 * mc$P_substitution_exact,
      method = "Transferred values"
    ),
    data.frame(
      delta = mc$delta, p = 100 * mc$P_joint_exact,
      method = "Joint value"
    )
  )
  exact_curves$method <- factor(
    exact_curves$method,
    levels = c("Joint value", "Transferred values")
  )
  mc_points <- rbind(
    data.frame(
      delta = mc$delta, p = 100 * mc$substitution_p,
      method = "Transferred values"
    ),
    data.frame(
      delta = mc$delta, p = 100 * mc$joint_p,
      method = "Joint value"
    )
  )
  mc_points$method <- factor(
    mc_points$method,
    levels = c("Joint value", "Transferred values")
  )
  mc_points <- subset(
    mc_points,
    abs(delta - 0.10) < 1e-8 |
      abs(delta - 0.30) < 1e-8 |
      abs(delta - 0.50) < 1e-8
  )
  focus_points <- subset(mc_points, abs(delta - 0.30) < 1e-8)
  focus_joint <- subset(focus_points, method == "Joint value")
  focus_joint$label <- sprintf("%.1f%%", focus_joint$p)
  focus_substitution <- subset(focus_points, method == "Transferred values")
  focus_substitution$label <- sprintf(
    "%s / %s",
    format(round(focus_substitution$p / 100 * n_rep), big.mark = ",", trim = TRUE),
    format(n_rep, big.mark = ",", trim = TRUE)
  )

  pB <- ggplot(exact_curves, aes(delta, p, colour = method)) +
    annotate("rect",
      xmin = a, xmax = 2 * a,
      ymin = -Inf, ymax = Inf, fill = col_shade
    ) +
    annotate(
      "text", x = 1.5 * a, y = 34,
      label = "Joint adequacy\nwithout transfer",
      colour = "grey40", size = 2.8, lineheight = 0.95
    ) +
    geom_line(aes(linetype = method), linewidth = 1.05) +
    geom_point(data = mc_points, aes(shape = method), size = 1.7, alpha = 0.8) +
    geom_point(
      data = focus_joint,
      size = 4.0, shape = 21, stroke = 1.1, fill = "white",
      show.legend = FALSE
    ) +
    geom_point(
      data = focus_substitution,
      size = 4.0, shape = 24, stroke = 1.1, fill = "white",
      show.legend = FALSE
    ) +
    geom_text(
      data = focus_joint,
      aes(label = label),
      nudge_x = -0.012, nudge_y = -10,
      hjust = 1, fontface = "bold", size = 3.0,
      show.legend = FALSE
    ) +
    geom_text(
      data = focus_substitution,
      aes(label = label),
      nudge_x = 0.018, nudge_y = 13,
      hjust = 0, fontface = "bold", size = 3.0,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = method_colours) +
    scale_shape_manual(values = method_shapes) +
    scale_linetype_manual(values = method_linetypes) +
    scale_x_continuous(
      breaks = seq(0, 0.6, by = 0.1), limits = c(0, 0.6),
      expand = expansion(mult = c(0, 0.015)),
      sec.axis = sec_axis(
        ~ . / a,
        breaks = 0:3,
        name = expression(q == delta / delta[sub])
      )
    ) +
    scale_y_continuous(
      breaks = c(0, 50, 100), limits = c(-8, 102),
      expand = c(0, 0)
    ) +
    labs(
      title = "Joint adequacy persists after transfer fails",
      x = "True difference between task means",
      y = "Pass rate (%)",
      colour = NULL, shape = NULL, linetype = NULL
    ) +
    simulation_theme() +
    theme(
      legend.position = "none",
      axis.title.x.top = element_text(size = 9.5, margin = margin(b = 5)),
      axis.text.x.top = element_text(size = 8.5),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.35)
    )

  fig <- (pA / pB) +
    plot_layout(heights = c(0.74, 1.26)) +
    plot_annotation(tag_levels = "A")

  save_figure(fig, file.path(figure_dir, "figure3_gaussian_simulation.pdf"),
              height_mm = 118)

  invisible(list(curves = mc, audit = audit, figure = fig))
}

if (sys.nframe() == 0) run_simulation1()
