# Data generator for every block of the univariate RKCFC revision.
#
# Implements manuscript Equation (17):
#
#   Y_ij,(k) = mu_(k)(t_ij) + sum_{m=1}^{M} xi_im,(k) phi_m,(k)(t_ij) + eps_ij
#
# Everything is read from study/config. Nothing here is hard-coded, so changing
# a parameter never means editing this file.
#
# THREE THINGS THAT ARE EASY TO GET WRONG, AND ARE GUARDED
#
# 1. The design's arm numbering is NOT the package's gen_scores() numbering.
#    Design arm 4 is `t nu = 1.5`, which is gen_scores arm 6. Arms 5 and 6 of
#    the design use GAUSSIAN scores (gen_scores arm 1) and put their departure
#    in the errors and the trajectories instead. The mapping lives in
#    config/arms.csv and is asserted by config/validate_config.R.
#
# 2. Observation counts are REDRAWN, not clamped, and the floor is 4.
#    The ORIGINAL script drew Poisson(7) and DELETED subjects outside 4-12,
#    drawing replacements until 500 survived. So the
#    analysed distribution was Poisson(7) conditioned on [4,12], the stated
#    n = 500 is correct, and the manuscript's "yielding 4-12" describes a
#    filter rather than a property of the Poisson.
#
#    The revision KEEPS the floor and DROPS the ceiling. A floor is justified:
#    predicted_yhat() cannot fit a singleton and very sparse subjects are
#    unstable. A ceiling is not - discarding well-observed subjects only throws
#    information away. Clamping would leave a spike at the floor; redrawing
#    gives a genuine truncated Poisson.
#
#    LIMITATION TO STATE IN THE MANUSCRIPT: with a floor of 4 the design never
#    tests subjects with 2 or 3 observations, which is arguably the sparse
#    regime the title refers to.
#
# 3. The contaminant is EXCLUDED from the true latent curve. See below.

# THE LIBRARY PATH IS NOT INCIDENTAL. Two different packages are currently named
# RKCFC: the standalone univariate build that this study needs, and the combined
# univariate+multivariate build that preceded the split into RKCFC and MRKCFC. They are NOT
# interchangeable -- `gen_scores` is exported by the standalone and internal in
# the combined one, so a plain library(RKCFC) picks up whichever happens to be
# installed and fails, or worse, silently runs the wrong code.
#
# A study-local library is used when present: a folder `lib` here, or the
# folder named by the environment variable RKCFC_LIB. Otherwise the package is
# loaded from the default library, where RKCFC 0.2.0 should be installed.
local({
  lib <- Sys.getenv("RKCFC_LIB", unset = "lib")
  if (dir.exists(lib)) .libPaths(c(normalizePath(lib), .libPaths()))
})
suppressMessages(library(RKCFC))
stopifnot(
  "the STANDALONE RKCFC is required; the combined build keeps gen_scores internal" =
    "gen_scores" %in% getNamespaceExports("RKCFC"),
  "this is the combined build, not the standalone one" =
    !("MKCFC" %in% getNamespaceExports("RKCFC"))
)
source(file.path("config", "validate_config.R"))

# ---------------------------------------------------------------------------
# Per-cell seed. Derived from the master seed plus the cell identity, so any
# scenario x replicate reproduces in isolation regardless of execution order or
# which tasks are re-run after a failure. This is what makes an HPC array
# restartable without replaying the whole grid.
cell_seed <- function(seed_master, scenario_id, rep) {
  key <- paste(scenario_id, rep, sep = "|")
  code <- utf8ToInt(key)
  h <- sum(code * seq_along(code) * 7919)
  as.integer((as.numeric(seed_master) + h) %% 2147483647)
}

# ---------------------------------------------------------------------------
# Observation counts, in two modes.
#
#   "fixed"  n_i ~ Poisson(lambda), truncated below at obs_min by REDRAWING.
#            The count is independent of the follow-up window, so a window level
#            using this mode varies WHERE observations fall but not HOW MANY.
#            That isolates coverage from count, which is a cleaner experiment
#            but is not the real mechanism.
#
#   "rate"   n_i ~ Poisson(lambda * T_i). Visits arrive at a rate, so shorter
#            follow-up means FEWER observations as well as narrower coverage --
#            which is what shorter follow-up actually is. This is the only mode
#            that reproduces the over-dispersion real cohorts show: a mixed
#            Poisson has var/mean > 1, whereas plain truncation pushes it BELOW
#            1, in the wrong direction.
#
# Rows using "rate" relax the floor to 2 via obs_min_override. That is not a
# convenience. Demanding four observations from a subject followed for four
# months forces short-follow-up subjects back up to the floor, which is exactly
# the selection that affects the cohort analysis: the submitted manuscript
# excluded patients with fewer than four CEA measurements, and short follow-up
# is why they had few. A floor of 4 here would make the simulation REPRODUCE
# that bias rather than test it.
draw_counts <- function(n, lambda, obs_min, Tmax = NULL, mode = "fixed") {
  lam <- if (mode == "rate") lambda * Tmax else rep(lambda, n)
  stopifnot(length(lam) == n, all(lam > 0))
  ni <- stats::rpois(n, lam)
  short <- which(ni < obs_min)
  while (length(short)) {
    ni[short] <- stats::rpois(length(short), lam[short])
    short <- short[ni[short] < obs_min]
  }
  as.integer(ni)
}

# ---------------------------------------------------------------------------
# Follow-up windows. "common" gives every subject the whole domain; "class"
# gives class 2 a systematically shorter window, so follow-up is INFORMATIVE
# about class -- mirroring the CEA cohort, where patients who die earlier accrue
# fewer measurements over a shorter period. That is the harder question
# about subject-specific follow-up windows, and the one that connects to the
# guarantee-time problem of the cohort.
draw_windows <- function(cfg, labels, window) {
  n <- length(labels)
  switch(window,
    common = rep(gnum(cfg, "t_max"), n),
    # Windows vary but carry NO information about class. This is the control:
    # it isolates the effect of uneven follow-up from the effect of follow-up
    # being informative, which the class-dependent level below confounds.
    noninf = stats::runif(n, gnum(cfg, "window_ni_lo"), gnum(cfg, "window_ni_hi")),
    class  = {
      lo <- c(gnum(cfg, "window_c1_lo"), gnum(cfg, "window_c2_lo"))
      hi <- c(gnum(cfg, "window_c1_hi"), gnum(cfg, "window_c2_hi"))
      stats::runif(n, lo[labels], hi[labels])
    },
    stop("unknown window: ", window, call. = FALSE))
}

# ---------------------------------------------------------------------------
# Measurement errors, calibrated on IQR so that the t3 arm carries the same
# interquartile spread as N(0, sd^2) and differs only in its tails. The factor
# below is 0.1 * IQR_SCALE[5] from the package's sim_arms.R, i.e. the same
# calibration the score arms use, so errors and scores are matched the same way.
draw_error <- function(n, sd, dist) {
  switch(dist,
    gaussian = stats::rnorm(n, 0, sd),
    t3       = sd * (stats::qnorm(0.75) - stats::qnorm(0.25)) /
                    (stats::qt(0.75, 3) - stats::qt(0.25, 3)) * stats::rt(n, 3),
    stop("unknown error_dist: ", dist, call. = FALSE)
  )
}

# ---------------------------------------------------------------------------
# The contaminating function for arm 6: a smooth squared-exponential Gaussian
# process draw, orthogonalised against EVERY class's eigenfunctions so that the
# contaminant cannot be absorbed by either class's own basis, then normalised to
# unit norm on the grid.
contam_shape <- function(tgrid, lengthscale, Phi_all) {
  D <- outer(tgrid, tgrid, "-")
  S <- exp(-(D^2) / (2 * lengthscale^2)) + diag(1e-8, length(tgrid))
  g <- as.vector(crossprod(chol(S), stats::rnorm(length(tgrid))))
  for (j in seq_len(ncol(Phi_all))) {                 # Gram-Schmidt
    p <- Phi_all[, j]
    g <- g - sum(g * p) / sum(p * p) * p
  }
  g / sqrt(mean(g^2))
}

# ---------------------------------------------------------------------------
#' Generate one replicate.
#'
#' @return A list with `Ly`, `Lt` (the observed data), `labels` (true classes),
#'   `xi` (true scores), `X_true` (true latent curves on the dense grid),
#'   `tgrid`, `contaminated` (logical), and `n_obs`.
#'
#' ON THE TRUE LATENT CURVE, WHICH IS A DELIBERATE CHOICE.
#'
#' For contaminated subjects, `X_true` holds the CLEAN class signal and excludes
#' the contaminant, while `Ly` includes it. The alternative -- treating the
#' contaminated curve as the truth -- would reward a method for chasing the
#' outlier and penalise one for resisting it, which inverts the property the arm
#' exists to measure. Contamination is a nuisance drawn from a different
#' process, not part of the estimand, which is the standard framing in robust
#' statistics. Consequence to report: the 5 percent of contaminated subjects
#' carry a large true-curve error for every method, so report that metric for
#' the contaminated and uncontaminated subsets separately rather than pooling.
generate_replicate <- function(cfg, arm_row, d, rep, scenario_id,
                               contrast = "parallel", M_true = 2L,
                               n = NULL, pi = NULL, lambda = NULL,
                               error_sd = NULL, window = "common",
                               count_mode = "fixed", obs_min_override = NA) {
  set.seed(cell_seed(gint(cfg, "seed_master"), scenario_id, rep))

  n        <- if (is.null(n)) gint(cfg, "n") else as.integer(n)
  pi       <- if (is.null(pi)) gvec(cfg, "pi") else pi
  lambda   <- if (is.null(lambda)) gnum(cfg, "obs_poisson_lambda") else lambda
  error_sd <- if (is.null(error_sd)) gnum(cfg, "error_sd") else error_sd
  # K is overridable because the K-recovery block needs true K in {1, 2, 3} and
  # every other block needs the configured 2. It is taken from `pi` when that is
  # supplied, so the two can never disagree -- passing a length-3 pi with a
  # config K of 2 previously reached sample.int() and died there with "incorrect
  # number of probabilities", which names the symptom and not the cause.
  K        <- if (!is.null(pi)) length(pi) else gint(cfg, "K_true")
  stopifnot("the design defines classes 1..3 only; see config/basis.csv" =
              K >= 1L && K <= 3L)
  obs_min  <- if (is.na(obs_min_override)) gint(cfg, "obs_min") else as.integer(obs_min_override)
  M_true   <- as.integer(M_true)

  tgrid <- seq(gnum(cfg, "t_min"), gnum(cfg, "t_max"),
               length.out = gint(cfg, "grid_n"))

  # class membership, true basis, true scores -------------------------------
  labels <- sample.int(K, n, replace = TRUE, prob = pi)
  mu  <- lapply(seq_len(K), function(k) eval_mean(cfg, k, tgrid, d, contrast))
  Phi <- lapply(seq_len(K), function(k) eval_eigen(cfg, k, M_true, tgrid, "parallel"))
  xi  <- lapply(seq_len(K), function(k)
                gen_scores(n, arm_row$score_arm, Z_matrix(cfg, k, M_true)))

  # true latent curves, CLEAN by construction --------------------------------
  X_true <- matrix(0, n, length(tgrid))
  for (i in seq_len(n)) {
    k <- labels[i]
    X_true[i, ] <- mu[[k]] + as.vector(Phi[[k]] %*% xi[[k]][i, ])
  }

  # contamination ------------------------------------------------------------
  contaminated <- rep(FALSE, n)
  X_obs <- X_true
  if (isTRUE(arm_row$contaminate)) {
    Phi_all <- do.call(cbind, Phi)
    rms <- vapply(seq_len(K), function(k)
      sqrt(mean(rowSums(sweep(Phi[[k]]^2, 2L, diag(Z_matrix(cfg, k, M_true)), "*")))),
      numeric(1))
    contaminated <- stats::runif(n) < gnum(cfg, "contam_frac")
    a <- gnum(cfg, "contam_mult")
    for (i in which(contaminated)) {
      g <- contam_shape(tgrid, gnum(cfg, "contam_lengthscale"), Phi_all)
      X_obs[i, ] <- X_true[i, ] + a * rms[labels[i]] * g
    }
  }

  # sparse irregular sampling within each subject's window --------------------
  Tmax  <- draw_windows(cfg, labels, window)
  n_obs <- draw_counts(n, lambda, obs_min, Tmax, count_mode)
  Ly <- Lt <- vector("list", n)
  for (i in seq_len(n)) {
    avail <- which(tgrid <= Tmax[i])
    if (length(avail) < n_obs[i]) n_obs[i] <- length(avail)
    idx <- sort(sample(avail, n_obs[i]))
    Lt[[i]] <- tgrid[idx]
    Ly[[i]] <- X_obs[i, idx] + draw_error(n_obs[i], error_sd, arm_row$error_dist)
  }

  list(Ly = Ly, Lt = Lt, labels = labels, xi = xi, X_true = X_true,
       tgrid = tgrid, contaminated = contaminated, n_obs = n_obs,
       scenario_id = scenario_id, rep = rep)
}
