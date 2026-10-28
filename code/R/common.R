# Settings and helpers shared by the analysis scripts.
#
# All scripts are run from the repository root, either through code/run_all.R
# or one at a time, e.g. Rscript --vanilla code/R/simulation1_gaussian.R.

project_root <- normalizePath(".")

# --vanilla skips .Rprofile, so the library created by renv::restore() is
# loaded here. Processes started by run_all.R inherit it from the parent.
# Disabling the renv sandbox keeps the system library visible.
Sys.setenv(RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
if (requireNamespace("renv", quietly = TRUE) &&
    dir.exists(renv::paths$library(project = project_root)) &&
    !identical(Sys.getenv("RENV_PROJECT"), project_root)) {
  renv::load(project = project_root, quiet = TRUE)
}

required_packages <- c(
  "digest", "dplyr", "ggplot2", "patchwork", "purrr", "readr",
  "scales", "stringr", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       ". See code/README.md.", call. = FALSE)
}
if (!capabilities("cairo")) {
  stop("The figures need an R build with Cairo support.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

figure_dir <- "figures"
supplement_figure_dir <- "supplement_figures"
results_dir <- file.path("code", "results")
dir.create(figure_dir, showWarnings = FALSE)
dir.create(supplement_figure_dir, showWarnings = FALSE)

# Figure colours. Blue, orange and green are from the Okabe-Ito palette.
blue <- "#0072B2"
orange <- "#D55E00"
green <- "#009E73"
purple <- "#332288"
charcoal <- "#3A3A3A"

# Wilson score interval for a binomial proportion.
wilson_interval <- function(successes, trials, z = 1.96) {
  p <- successes / trials
  denominator <- 1 + z^2 / trials
  centre <- (p + z^2 / (2 * trials)) / denominator
  half_width <- z * sqrt(
    p * (1 - p) / trials + z^2 / (4 * trials^2)
  ) / denominator
  cbind(
    lower = pmax(0, centre - half_width),
    upper = pmin(1, centre + half_width),
    half_width = half_width
  )
}

# Call simulate_one() n times and stack the named vectors it returns into a
# matrix with one row per replication.
replicate_rows <- function(n, simulate_one) {
  do.call(rbind, replicate(n, simulate_one(), simplify = FALSE))
}

# Theme for the conceptual figures and the supplementary figures.
paper_theme <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      text = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      axis.title = element_text(size = base_size),
      plot.title = element_text(size = base_size, face = "plain", hjust = 0,
                                margin = margin(b = 5)),
      plot.tag = element_text(size = base_size + 2, face = "bold"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.width = unit(1.3, "cm"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

# Theme for the simulation figures (Figures 3 and 5).
simulation_theme <- function() {
  theme_classic(base_size = 10.5, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(size = 9, colour = "black"),
      axis.title = element_text(size = 10, colour = "black"),
      plot.title = element_text(
        size = 10.5, hjust = 0, face = "bold",
        margin = margin(b = 4)
      ),
      plot.tag = element_text(size = 11.5, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 9),
      legend.key.width = unit(1.1, "cm"),
      legend.margin = margin(t = 1),
      plot.margin = margin(5, 8, 4, 6),
      plot.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA)
    )
}

save_figure <- function(plot, path, width_mm = 170, height_mm) {
  ggsave(
    filename = path,
    plot = plot,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    device = cairo_pdf,
    bg = "white"
  )
  message("Saved ", path)
  invisible(path)
}
