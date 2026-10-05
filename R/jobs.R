# The scenario x replicate job list for each fixed-K block.
#
# Shared by run_cell.R, which executes a slice of it, and combine.R, which needs
# it to know how many cells SHOULD be present.
#
# That second use is the reason this is a separate file. A completeness check
# that derives the expected count from the rows it can see is circular: it will
# report 100 percent on a single task out of an array of a hundred, which is
# precisely the failure it exists to catch. Expected comes from the config.

# Vectorised: it is called with a whole column of levels, and an if() there
# silently uses only the first element in older R and errors in current R.
d_of <- function(cfg, level) {
  ifelse(level == "hard", gnum(cfg, "d_hard"), gnum(cfg, "d_easy"))
}

build_jobs <- function(cfg, block, n_reps = NA_integer_) {
  if (block == "pilot") {
    s <- cfg$pilot
    reps <- seq_len(if (is.na(n_reps)) gint(cfg, "n_reps_pilot") else n_reps)
    j <- expand.grid(row = seq_len(nrow(s)), rep = reps)
    out <- data.frame(scenario_id = s$scenario_id[j$row], arm = s$arm[j$row],
                      d = s$d[j$row], rep = j$rep, contrast = "parallel",
                      M_true = 2L, lambda = NA, error_sd = NA, n = NA,
                      pi = NA, window = "common", count_mode = "fixed",
                      obs_min_override = NA, include = "include_pilot",
                      stringsAsFactors = FALSE)
  } else if (block == "core") {
    s <- cfg$core
    reps <- seq_len(if (is.na(n_reps)) gint(cfg, "n_reps_core") else n_reps)
    j <- expand.grid(row = seq_len(nrow(s)), rep = reps)
    out <- data.frame(scenario_id = s$scenario_id[j$row], arm = s$arm[j$row],
                      d = d_of(cfg, s$d_level[j$row]), rep = j$rep,
                      contrast = "parallel", M_true = 2L, lambda = NA,
                      error_sd = NA, n = NA, pi = NA, window = "common",
                      count_mode = "fixed", obs_min_override = NA,
                      include = "include_core", stringsAsFactors = FALSE)
  } else if (block == "variations") {
    s <- cfg$variations
    arms <- c(2L, 4L)          # decided 2026-08-30: BOTH, so the factors are
    reps <- seq_len(if (is.na(n_reps)) gint(cfg, "n_reps") else n_reps)
    j <- expand.grid(row = seq_len(nrow(s)), arm = arms, rep = reps)
    out <- data.frame(
      scenario_id = paste0(s$config_id[j$row], "_a", j$arm), arm = j$arm,
      # the mean-identical contrast IS d = 0; it is not a change of basis
      d = ifelse(s$contrast[j$row] == "identical", 0, d_of(cfg, s$d_level[j$row])),
      rep = j$rep, contrast = s$contrast[j$row], M_true = s$M_true[j$row],
      lambda = s$lambda[j$row], error_sd = s$error_sd[j$row], n = s$n[j$row],
      pi = s$pi[j$row], window = s$window[j$row],
      count_mode = s$count_mode[j$row],
      obs_min_override = s$obs_min_override[j$row],
      include = "include_var", stringsAsFactors = FALSE)
  } else {
    stop("unknown block: ", block, ". Use pilot, core or variations.", call. = FALSE)
  }
  out$block <- block
  out
}

# How many rows a complete block should produce: cells x the methods it runs.
expected_rows <- function(cfg, block, n_reps = NA_integer_) {
  j <- build_jobs(cfg, block, n_reps)
  nrow(j) * sum(as.logical(cfg$methods[[unique(j$include)]]))
}
