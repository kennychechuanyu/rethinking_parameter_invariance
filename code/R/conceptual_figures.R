# Figures 1 and 2: acceptable regions and the geometry of the joint compromise.
# Both figures are schematic; nothing here is simulated.

source("code/R/common.R")

# Figure 1. Each task has an acceptable region of radius 0.20 around its
# optimum. The three panels move the Task B optimum away from Task A: each
# optimum inside the other region, regions overlapping without containing the
# other optimum, and no overlap.
make_figure1 <- function() {
  scenarios <- data.frame(
    panel = factor(
      c("Bidirectional transfer", "Accommodation only", "Neither"),
      levels = c("Bidirectional transfer", "Accommodation only", "Neither")
    ),
    theta_A = c(0.30, 0.30, 0.30),
    theta_B = c(0.40, 0.60, 0.80),
    radius_A = 0.20,
    radius_B = 0.20
  )

  scenarios$A_left <- scenarios$theta_A - scenarios$radius_A
  scenarios$A_right <- scenarios$theta_A + scenarios$radius_A
  scenarios$B_left <- scenarios$theta_B - scenarios$radius_B
  scenarios$B_right <- scenarios$theta_B + scenarios$radius_B
  scenarios$overlap_left <- pmax(scenarios$A_left, scenarios$B_left)
  scenarios$overlap_right <- pmin(scenarios$A_right, scenarios$B_right)
  overlap <- subset(scenarios, overlap_left < overlap_right)

  p <- ggplot() +
    geom_rect(
      data = overlap,
      aes(xmin = overlap_left, xmax = overlap_right,
          ymin = 0.22, ymax = 0.78),
      fill = "#D9EEE8", colour = "#72AA9A", linewidth = 0.4
    ) +
    geom_segment(
      data = scenarios,
      aes(x = A_left, xend = A_right, y = 0.32, yend = 0.32),
      colour = blue, linewidth = 4.2, lineend = "butt"
    ) +
    geom_segment(
      data = scenarios,
      aes(x = B_left, xend = B_right, y = 0.68, yend = 0.68),
      colour = orange, linewidth = 4.2, lineend = "butt"
    ) +
    geom_point(
      data = scenarios,
      aes(theta_A, 0.32), shape = 21, size = 4.1,
      stroke = 1.0, colour = blue, fill = "white"
    ) +
    geom_point(
      data = scenarios,
      aes(theta_B, 0.68), shape = 21, size = 4.1,
      stroke = 1.0, colour = orange, fill = "white"
    ) +
    geom_text(
      data = scenarios,
      aes(theta_A, 0.12), label = expression(theta[A]^"*"),
      colour = blue, size = 4.0
    ) +
    geom_text(
      data = scenarios,
      aes(theta_B, 0.88), label = expression(theta[B]^"*"),
      colour = orange, size = 4.0
    ) +
    geom_text(
      data = subset(scenarios, panel == "Bidirectional transfer"),
      aes(x = 0.04, y = 0.32), label = "Task A",
      colour = blue, hjust = 1, size = 4.0
    ) +
    geom_text(
      data = subset(scenarios, panel == "Bidirectional transfer"),
      aes(x = 0.04, y = 0.68), label = "Task B",
      colour = orange, hjust = 1, size = 4.0
    ) +
    facet_wrap(~panel, nrow = 1) +
    scale_x_continuous(limits = c(-0.30, 1.12), breaks = seq(0, 1, 0.2)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = expression("Focal value  " * theta), y = NULL) +
    paper_theme(11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(size = 10.5),
      axis.title.x = element_text(size = 11.5, margin = margin(t = 5)),
      strip.background = element_blank(),
      strip.text = element_text(size = 11.5, face = "bold",
                                margin = margin(b = 6)),
      panel.spacing.x = unit(0.7, "lines"),
      legend.position = "none"
    )

  save_figure(p, file.path(figure_dir, "figure1_transfer_regions.pdf"),
              height_mm = 64)
}

# Figure 2. (A) With equal quadratic losses, the joint value lies halfway
# between the task optima. (B) Half the distance to the Task B optimum means
# one quarter of the Task B excess loss. (C) A common tolerance cuts the more
# curved loss into the narrower acceptable interval.
make_figure2 <- function() {
  theta_A <- 0
  theta_B <- 2
  theta_joint <- 1
  x <- seq(-0.35, 2.35, length.out = 500)

  curves <- data.frame(
    theta = rep(x, 3),
    loss = c(
      0.5 * (x - theta_A)^2,
      0.5 * (x - theta_B)^2,
      0.5 * (x - theta_A)^2 + 0.5 * (x - theta_B)^2
    ),
    curve = factor(
      rep(c("Task A", "Task B", "Combined"), each = length(x)),
      levels = c("Task A", "Task B", "Combined")
    )
  )
  curve_colours <- c(
    "Task A" = blue, "Task B" = orange, "Combined" = charcoal
  )
  curve_types <- c("Task A" = "solid", "Task B" = "22", "Combined" = "31")

  minima <- data.frame(
    theta = c(theta_A, theta_B, theta_joint),
    loss = c(0, 0, 1),
    curve = factor(
      c("Task A", "Task B", "Combined"),
      levels = levels(curves$curve)
    )
  )

  panel_a <- ggplot(curves, aes(theta, loss, colour = curve, linetype = curve)) +
    geom_line(linewidth = 0.95) +
    annotate(
      "segment", x = theta_joint, xend = theta_joint, y = 0, yend = 1,
      colour = charcoal, linewidth = 0.4, linetype = "dotted"
    ) +
    geom_point(data = minima, size = 2.7, show.legend = FALSE) +
    annotate("text", theta_A, -0.13, label = expression(theta[A]^"*"),
             colour = blue, size = 3.2) +
    annotate("text", theta_joint, -0.13,
             label = expression(theta[joint]^"*"),
             colour = charcoal, size = 3.2) +
    annotate("text", theta_B, -0.13, label = expression(theta[B]^"*"),
             colour = orange, size = 3.2) +
    annotate("text", 1.66, 1.70, label = expression(L[A]),
             colour = blue, size = 3.1) +
    annotate("text", 0.34, 1.70, label = expression(L[B]),
             colour = orange, size = 3.1) +
    annotate("text", 1.00, 1.27, label = expression(L[A] + L[B]),
             colour = charcoal, size = 3.1) +
    scale_colour_manual(values = curve_colours) +
    scale_linetype_manual(values = curve_types) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = 0:2) +
    coord_cartesian(xlim = c(-0.4, 2.4), ylim = c(-0.2, 2.25)) +
    labs(title = "The joint value is a compromise",
         x = expression(theta), y = "Loss") +
    paper_theme(9) +
    theme(
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )

  task_b_loss <- data.frame(theta = x, loss = 0.5 * (x - theta_B)^2)
  evaluations <- data.frame(
    theta = c(theta_A, theta_joint, theta_B),
    loss = c(2, 0.5, 0),
    value = factor(
      c("Substitution", "Joint", "Task B optimum"),
      levels = c("Substitution", "Joint", "Task B optimum")
    )
  )
  evaluation_colours <- c(
    "Substitution" = blue,
    "Joint" = charcoal,
    "Task B optimum" = orange
  )

  panel_b <- ggplot(task_b_loss, aes(theta, loss)) +
    geom_line(colour = charcoal, linewidth = 0.95) +
    geom_segment(
      data = evaluations,
      aes(x = theta, xend = theta, y = 0, yend = loss, colour = value),
      linewidth = 0.45, linetype = "dotted"
    ) +
    geom_point(data = evaluations, aes(colour = value), size = 2.8) +
    annotate("text", theta_A + 0.08, 2.10,
             label = expression(Delta * L[B]^sub), colour = blue,
             hjust = 0, size = 3.2) +
    annotate("text", theta_joint + 0.10, 0.66,
             label = expression(Delta * L[B]^joint), colour = charcoal,
             hjust = 0, size = 3.2) +
    annotate(
      "text", 1.48, 1.43,
      label = expression(Delta * L[B]^joint == frac(1, 4) * Delta * L[B]^sub),
      colour = charcoal, size = 3.15
    ) +
    annotate("text", theta_A, -0.13, label = expression(theta[A]^"*"),
             colour = blue, size = 3.2) +
    annotate("text", theta_joint, -0.13,
             label = expression(theta[joint]^"*"),
             colour = charcoal, size = 3.2) +
    annotate("text", theta_B, -0.13, label = expression(theta[B]^"*"),
             colour = orange, size = 3.2) +
    scale_colour_manual(values = evaluation_colours) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = 0:2) +
    coord_cartesian(xlim = c(-0.4, 2.4), ylim = c(-0.2, 2.25)) +
    labs(
      title = "Half the distance, one-quarter the loss",
      x = expression(theta),
      y = "Excess loss on Task B"
    ) +
    paper_theme(9) +
    theme(
      axis.ticks.x = element_blank(),
      legend.position = "none"
    )

  # Panel C: curvatures 6 and 1.2 with tolerance 0.5.
  w_a <- 6
  w_b <- 1.2
  eps <- 0.5
  t_a <- -0.35; t_b <- 0.35
  r_a <- sqrt(2 * eps / w_a)
  r_b <- sqrt(2 * eps / w_b)
  xs <- seq(-1.55, 1.85, length.out = 481)
  curves <- rbind(
    data.frame(task = "A", x = xs, y = w_a / 2 * (xs - t_a)^2),
    data.frame(task = "B", x = xs, y = w_b / 2 * (xs - t_b)^2)
  )
  curves <- subset(curves, y <= 1.5)
  bars <- data.frame(
    task = c("A", "B"),
    xmin = c(t_a - r_a, t_b - r_b),
    xmax = c(t_a + r_a, t_b + r_b),
    y = c(-0.13, -0.27),
    optimum = c(t_a, t_b)
  )

  panel_c <- ggplot() +
    geom_hline(yintercept = eps, linetype = "22", colour = "grey45",
               linewidth = 0.4) +
    geom_line(data = curves, aes(x, y, colour = task), linewidth = 1.05) +
    geom_segment(
      data = bars,
      aes(x = xmin, xend = xmax, y = y, yend = y, colour = task),
      linewidth = 4.2, lineend = "butt"
    ) +
    geom_point(
      data = bars, aes(optimum, y), shape = 21, size = 3.6,
      stroke = 1.0, colour = c(blue, orange), fill = "white"
    ) +
    scale_colour_manual(values = c(A = blue, B = orange), guide = "none") +
    annotate("text", x = -0.35, y = 0.78, label = "Task A",
             colour = blue, size = 3.1, fontface = "bold") +
    annotate("text", x = 0.35, y = 0.30, label = "Task B",
             colour = orange, size = 3.1, fontface = "bold") +
    annotate("text", x = 1.78, y = eps + 0.09, label = "tolerance~epsilon",
             parse = TRUE, colour = "grey35", size = 3.0, hjust = 1) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL, limits = c(-0.38, 1.5),
                       expand = c(0, 0)) +
    labs(x = expression("Focal value " * theta), y = "Excess loss") +
    paper_theme(9)

  figure <- ((panel_a | panel_b) / panel_c) +
    plot_annotation(tag_levels = "A") +
    plot_layout(heights = c(1, 0.72))
  save_figure(figure, file.path(figure_dir, "figure2_compromise_geometry.pdf"),
              height_mm = 112)
}

if (sys.nframe() == 0) {
  make_figure1()
  make_figure2()
}
