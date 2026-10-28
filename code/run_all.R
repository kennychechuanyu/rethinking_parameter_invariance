#!/usr/bin/env Rscript

# Rebuild all figures and check the saved results. Run from the repository root:
#
#   Rscript --vanilla code/run_all.R                          standard rebuild
#   Rscript --vanilla code/run_all.R --verify                 check saved results only
#   Rscript --vanilla code/run_all.R --resimulate             also rerun Simulations 2A and 2B
#   Rscript --vanilla code/run_all.R --resimulate-validation  also rerun the finite-validation analysis
#   Rscript --vanilla code/run_all.R --resimulate-all         rerun every simulation
#
# The standard rebuild reruns the fast analyses and draws the remaining figures
# from the saved Monte Carlo results. Each analysis runs in its own R process.
# The saved results are compared with the SHA-256 checksums in
# code/results/manifest.csv before and after the rebuild.

if (!file.exists("code/run_all.R")) {
  stop("Run from the repository root: Rscript --vanilla code/run_all.R",
       call. = FALSE)
}

usage <- paste(
  "Usage: Rscript --vanilla code/run_all.R [--verify | --resimulate |",
  "--resimulate-validation | --resimulate-all]\n"
)
args <- commandArgs(trailingOnly = TRUE)
unknown <- setdiff(args, c(
  "--verify", "--resimulate", "--resimulate-validation", "--resimulate-all",
  "--update-manifest", "--help"
))
if (length(unknown)) {
  stop("Unknown option: ", toString(unknown), "\n", usage, call. = FALSE)
}
if ("--help" %in% args) {
  cat(usage)
  quit(save = "no")
}

source("code/R/common.R")
suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(readr)
  library(tibble)
})

manifest_path <- file.path(results_dir, "manifest.csv")
software_path <- file.path(results_dir, "software_versions.csv")

# Manifest ----

analysis_scripts <- c(
  simulation1 = "simulation1_gaussian.R",
  simulation2a = "simulation2a_psychometric.R",
  simulation2b = "simulation2b_directional.R",
  finite_validation = "finite_validation.R",
  supplement_checks = "supplement_checks.R"
)
analysis_seeds <- c(
  simulation1 = "20240115",
  simulation2a = "20240115",
  simulation2b = "20260823",
  finite_validation = "20260821",
  supplement_checks = "20260824"
)

# Monte Carlo output that is kept because it is slow to regenerate. Every
# other file is a summary of these, a fast simulation, or a calculation.
retained_files <- c(
  "simulation2a/finite_sample_by_delta.csv",
  "simulation2a/sensitivity.csv",
  "simulation2a/recovery_benchmark.csv",
  "simulation2a/parameter_boundary_checks.csv",
  "simulation2a/focus_delta_050.csv.gz",
  "simulation2b/finite_sample_replications.csv.gz",
  "finite_validation/finite_validation_replications.csv.gz",
  "finite_validation/parameter_boundary_checks.csv"
)

# Files that use no random numbers; the manifest records no seed for them.
deterministic_files <- c(
  "simulation2a/population_results.csv",
  "simulation2a/reference_tolerance_anchor.csv",
  "simulation2b/candidate_design_pairs.csv",
  "simulation2b/design_curvatures.csv",
  "simulation2b/population_boundaries.csv",
  "simulation2b/population_curves.csv",
  "simulation2b/profile_curvature_by_delta.csv",
  "simulation2b/reversed_orientation_boundaries.csv",
  "finite_validation/specification.csv",
  "supplement_checks/nuisance_policy_boundaries.csv"
)

result_files <- function() {
  files <- list.files(results_dir, recursive = TRUE, all.files = TRUE)
  setdiff(files, c("manifest.csv", "software_versions.csv"))
}

open_text <- function(path) {
  if (grepl("[.]gz$", path)) gzfile(path, "rt") else file(path, "rt")
}

count_lines <- function(path) {
  connection <- open_text(path)
  on.exit(close(connection))
  total <- 0L
  repeat {
    lines <- readLines(connection, n = 10000L, warn = FALSE)
    if (!length(lines)) break
    total <- total + length(lines)
  }
  total
}

header_line <- function(path) {
  connection <- open_text(path)
  on.exit(close(connection))
  readLines(connection, n = 1L, warn = FALSE)
}

sha256 <- function(paths) {
  vapply(paths, digest, character(1), file = TRUE, algo = "sha256")
}

build_manifest <- function() {
  files <- result_files()
  paths <- file.path(results_dir, files)
  analysis <- sub("/.*$", "", files)
  seed <- unname(analysis_seeds[analysis])
  seed[files %in% deterministic_files] <- NA_character_

  tibble(
    file = files,
    analysis = analysis,
    role = ifelse(
      files %in% retained_files,
      "retained Monte Carlo result", "derived summary"
    ),
    rows = pmax(0L, vapply(paths, count_lines, integer(1)) - 1L),
    columns = lengths(strsplit(
      vapply(paths, header_line, character(1)), ",", fixed = TRUE
    )),
    bytes = as.numeric(file.info(paths)$size),
    sha256 = sha256(paths),
    analysis_seed = seed,
    produced_by = unname(analysis_scripts[analysis])
  ) %>%
    arrange(analysis, file)
}

verify_results <- function(exclude = character()) {
  expected <- read_csv(manifest_path, show_col_types = FALSE) %>%
    filter(!analysis %in% exclude)
  present <- result_files()
  present <- present[!sub("/.*$", "", present) %in% exclude]

  missing <- setdiff(expected$file, present)
  unexpected <- setdiff(present, expected$file)
  if (length(missing) || length(unexpected)) {
    stop(
      "The results folder does not match the manifest.",
      if (length(missing)) paste0(" Missing: ", toString(missing), "."),
      if (length(unexpected)) paste0(" Not in manifest: ", toString(unexpected), "."),
      call. = FALSE
    )
  }

  paths <- file.path(results_dir, expected$file)
  changed <- file.info(paths)$size != expected$bytes |
    sha256(paths) != expected$sha256
  if (any(changed)) {
    stop("Saved results differ from the manifest: ",
         toString(expected$file[changed]), call. = FALSE)
  }
  message("Saved results match the manifest (", nrow(expected), " files).")
}

# Software ----

software_versions <- function() {
  packages <- required_packages
  tibble(
    component = c("R", "Platform", "Cairo PDF", packages, "RNG",
                  "renv.lock SHA-256"),
    version = c(
      as.character(getRversion()),
      R.version$platform,
      as.character(capabilities("cairo")),
      vapply(packages, function(p) as.character(packageVersion(p)), ""),
      paste(RNGkind(), collapse = " / "),
      digest("renv.lock", file = TRUE, algo = "sha256")
    )
  )
}

# Package versions and renv.lock must match the record exactly. A different R
# version or platform only produces a warning, and the results are still
# checked against their checksums.
check_software <- function() {
  recorded <- read_csv(software_path, show_col_types = FALSE)
  changed <- full_join(
    recorded, software_versions(),
    by = "component", suffix = c("_recorded", "_current")
  ) %>%
    filter(is.na(version_recorded) | is.na(version_current) |
             version_recorded != version_current)

  informational <- c("R", "Platform")
  strict <- filter(changed, !component %in% informational)
  if (nrow(strict)) {
    stop(
      "Installed software differs from code/results/software_versions.csv: ",
      toString(strict$component), ". Run renv::restore() first.",
      call. = FALSE
    )
  }
  if (nrow(changed)) {
    warning("Differs from the recorded environment: ",
            toString(changed$component), call. = FALSE)
  } else {
    message("Software versions match the record.")
  }
}

# Run ----

if ("--update-manifest" %in% args) {
  # For the authors, after an intended change to the saved results.
  write_csv(build_manifest(), manifest_path)
  write_csv(software_versions(), software_path)
  message("Updated ", manifest_path, " and ", software_path)
  quit(save = "no")
}

check_software()
if ("--verify" %in% args) {
  verify_results()
  quit(save = "no")
}

resimulate_main <- any(c("--resimulate", "--resimulate-all") %in% args)
resimulate_validation <- any(
  c("--resimulate-validation", "--resimulate-all") %in% args
)

# Simulation 1 and the supplementary checks take a few seconds and are always
# rerun. Saved results that will not be regenerated are checked first.
rerun <- c(
  "simulation1", "supplement_checks",
  if (resimulate_main) c("simulation2a", "simulation2b"),
  if (resimulate_validation) "finite_validation"
)
if (length(setdiff(names(analysis_scripts), rerun))) {
  verify_results(exclude = rerun)
}

run_script <- function(label, script, resimulate = FALSE) {
  cat("\n", label, "\n", sep = "")
  started <- proc.time()[["elapsed"]]
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", file.path("code", "R", script),
      if (resimulate) "--resimulate")
  )
  if (status != 0) stop(label, " failed.", call. = FALSE)
  cat(sprintf("Done in %.1f s.\n", proc.time()[["elapsed"]] - started))
}

run_script("Figures 1 and 2", "conceptual_figures.R")
run_script("Simulation 1 (Figure 3)", "simulation1_gaussian.R")
run_script("Simulation 2A (Figure 4)", "simulation2a_psychometric.R",
           resimulate_main)
run_script("Simulation 2B (Figure 5)", "simulation2b_directional.R",
           resimulate_main)
run_script("Supplementary checks (Figure S1)", "supplement_checks.R")
run_script("Finite-validation analysis (Figure S2)", "finite_validation.R",
           resimulate_validation)

cat("\n")
verify_results()
cat("All figures rebuilt.\n")
