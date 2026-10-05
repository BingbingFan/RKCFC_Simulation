# Estimation accuracy on the TRUE partition -- block E.
#
#   Rscript R/run_estim.R <task_id> [n_tasks] [outdir] [n_reps] [n_override]
#
# WHY THIS BLOCK HAS ITS OWN RUNNER, AND IS NOT A METHOD ROW.
#
# Every other block scores CLUSTERING, and clustering on the heavy-tailed arms is
# dominated by initialisation failure, which the 0.1.x start made
# worse. The estimator is then being judged through a failure mode that has
# nothing to do with it. Fitting under the TRUE partition removes the confound
# completely: no initialiser, no K-centers loop, no convergence, no stop_reason.
# What is left is the question the method is actually about -- does Kendall's tau
# estimate the class-specific mean and eigenfunctions better than PACE does?
#
# It is also cheap, because there is no iteration: one FPCA per class per
# replicate. The output shape is per (replicate, estimator, class, component),
# not per (cell, method), which is why it cannot share run_cell.R.
#
# FOUR ENDPOINTS, and each answers something specific:
#
#   mean_err     relative L2 error of mu_hat against the TRUE mu.
#                The mean is local-linear least squares and is NOT
#                robust, and this measures what that costs, arm by arm.
#   angle_m      angle between phi_hat_m and the true phi_m, in degrees. The
#                method's central claim. Meaningful HERE because the eigen-gap is
#                55.6 percent; with a small gap the angles would be unstable
#                and canonical correlation would be needed instead.
#   angle_sub    largest principal angle between the estimated and true spans.
#                Sign- and order-invariant companion to the above.
#   cancor       mean canonical correlation between estimated and true scores.
#                Invariant to any nonsingular linear transform of either set, so
#                it needs no alignment.
#
# NOTHING IS TRUNCATED: M_true = nK retained for this block, so the comparison is
# the full estimated eigenspace against the full true one. The truncation is a
# property of the MEASUREMENT and is deliberately NOT tied to nK_fit -- tying
# them would make the target move with the setting under test, which is the
# M_true / nK_fit conflation.

source(file.path("R", "methods.R"))

args    <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 1L
outdir  <- if (length(args) >= 3) args[3] else file.path("results", "estim")
n_reps  <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_
n_over  <- if (length(args) >= 5) as.integer(args[5]) else NA_integer_

cfg <- load_config()
K   <- gint(cfg, "K_true")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

M_EVAL <- 2L        # see "NOTHING IS TRUNCATED" above; NOT gint(cfg, "nK_fit")
reps   <- seq_len(if (is.na(n_reps)) gint(cfg, "n_reps") else n_reps)

# 6 arms x the HARD separation only. The easy separation adds nothing here: the
# estimation endpoints do not involve the class labels at all, so `d` shifts the
# means apart without changing how well either is estimated.
jobs <- expand.grid(arm = 1:6, rep = reps, KEEP.OUT.ATTRS = FALSE)
jobs$d <- gnum(cfg, "d_hard")
jobs$scenario_id <- sprintf("E%02d", jobs$arm)
jobs <- jobs[order(jobs$rep, jobs$arm), ]
mine <- jobs[(seq_len(nrow(jobs)) - 1L) %% n_tasks + 1L == task_id, , drop = FALSE]

est <- cfg$estimators
cat(sprintf("block estim | task %d of %d | %d of %d cells | %d estimators: %s\n",
            task_id, n_tasks, nrow(mine), nrow(jobs), nrow(est),
            paste(est$estimator_id, collapse = " ")))
if (!nrow(mine))
  stop(sprintf(paste("task %d of %d was assigned NO cells (the block has %d).",
                     "Pass the ORIGINAL array size as n_tasks when re-running one",
                     "task; with n_tasks = 1 every cell is reassigned."),
               task_id, n_tasks, nrow(jobs)), call. = FALSE)

# --- the row shape, defined ONCE ----------------------------------------------
#
# Same discipline as empty_scores() in methods.R, and for the same reason: a
# failure frame that differs from the success frame by one column takes down the
# whole task at the final rbind, after all of its compute has been spent.
empty_estim <- function(estimator = NA_character_, class = NA_integer_,
                        error = NA_character_) {
  data.frame(estimator = estimator, class = class, error = error,
             seconds = NA_real_, mean_err = NA_real_,
             angle_1 = NA_real_, angle_2 = NA_real_, angle_sub = NA_real_,
             cancor = NA_real_, n_k = NA_integer_, stringsAsFactors = FALSE)
}

# Angle in degrees between two vectors, sign-invariant. |cos| because an
# eigenfunction and its negation are the same direction, and the sign returned by
# eigen() is arbitrary.
vec_angle <- function(a, b) {
  ca <- abs(sum(a * b)) / sqrt(sum(a^2) * sum(b^2))
  acos(min(1, max(-1, ca))) * 180 / pi
}

# --- run ----------------------------------------------------------------------

rows <- vector("list", nrow(mine) * nrow(est) * K)
kk <- 0L
t_start <- Sys.time()

for (i in seq_len(nrow(mine))) {
  jb   <- mine[i, ]
  arm  <- cfg$arms[cfg$arms$arm == jb$arm, ]
  seed <- cell_seed(gint(cfg, "seed_master"), jb$scenario_id, jb$rep)

  dat <- tryCatch(
    generate_replicate(cfg, arm, d = jb$d, rep = jb$rep,
                       scenario_id = jb$scenario_id, M_true = M_EVAL,
                       n = if (is.na(n_over)) NULL else n_over),
    error = function(e) structure(conditionMessage(e), class = "genfail"))

  if (inherits(dat, "genfail")) {
    for (m in seq_len(nrow(est))) for (k in seq_len(K)) {
      kk <- kk + 1L
      rows[[kk]] <- cbind(jb, empty_estim(est$estimator_id[m], k,
                                          paste("GENERATION:", as.character(dat))),
                          row.names = NULL)
    }
    next
  }

  # The TRUE mean and eigenfunctions, on the same grid the fits are evaluated on.
  tg    <- dat$tgrid
  mu_t  <- lapply(seq_len(K), function(k) eval_mean(cfg, k, tg, jb$d))
  Phi_t <- lapply(seq_len(K), function(k) eval_eigen(cfg, k, M_EVAL, tg))

  for (m in seq_len(nrow(est))) {
    eid <- est$estimator_id[m]
    rob <- as.logical(est$robust[m])

    # E2 and E3 were built on 2026-09-07 and no longer
    # record a stated NOT IMPLEMENTED.
    set.seed(cell_seed(gint(cfg, "seed_master"),
                       paste0(jb$scenario_id, "|fit"), jb$rep))
    # Set ONLY the option the engine in use recognises. build_optns() already
    # picks the right name per branch; setting both put an `nK` into the fdapace
    # list, and fdapace VALIDATES option names, so every E0 fit failed with
    # "Invalid option names: nK" while the robust rows sailed through.
    optns <- build_optns(cfg, rob)
    if (rob) optns$nK <- M_EVAL else optns$methodSelectK <- M_EVAL

    # PACKAGE ITEM 1. E1M is E1 with the mean estimated robustly and NOTHING else
    # changed -- same Kendall scatter, same bandwidths, same grid -- so the E1M
    # against E1 contrast isolates the mean estimator and only that. The option
    # reaches KFPCA() through univar_FPCA_k(), which forwards
    # isTRUE(optns$meanrobust); an optns list without the field selects least
    # squares, which is why E0 and E1 keep reproducing their banked numbers.
    if (eid == "E1M") optns$meanrobust <- TRUE

    # E2: spherical / spatial-sign scatter. Signs each subject's deviation from
    # the mean rather than pairwise differences, at the SAME mean E1 uses, so
    # E2 against E1 isolates the scatter estimator.
    if (eid == "E2") optns$scatter <- "sign"

    # E3: normal-scores marginal transform, then ordinary PACE. The transform is
    # applied to the MEASUREMENTS, not to an option, which is why it appears here
    # rather than in `optns`.
    #
    # READ E3'S ROWS DIFFERENTLY FROM THE OTHERS. A monotone marginal transform
    # changes the estimand: the eigenfunctions of the transformed process are not
    # those of X, so E3's `mean_err` and angles measure how far the rank-transformed
    # eigenspace sits from the UNTRANSFORMED truth. That is a real and useful
    # quantity -- it prices what the transform distorts -- but it is NOT estimation
    # accuracy and must not be tabulated beside E0/E1/E1M as though it were.
    # A PARTITION is invariant to the transform, so the clustering comparison in
    # the settings block is where E3 is judged on equal terms.
    Ly_fit <- if (eid == "E3") normal_scores(dat$Ly) else dat$Ly

    t0 <- Sys.time()
    fits <- tryCatch(
      univar_FPCA_k(Ly_fit, dat$Lt, K, dat$labels, rob, optns),   # TRUE labels
      error = function(e) structure(conditionMessage(e), class = "efail"))
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (inherits(fits, "efail")) {
      for (k in seq_len(K)) {
        kk <- kk + 1L
        r <- empty_estim(eid, k, as.character(fits)); r$seconds <- secs
        rows[[kk]] <- cbind(jb, r, row.names = NULL)
      }
      next
    }

    for (k in seq_len(K)) {
      r <- empty_estim(eid, k); r$seconds <- secs
      r$n_k <- sum(dat$labels == k)
      out <- tryCatch({
        f <- fits[[k]]
        g   <- if (rob) f$RegGrid else f$workGrid
        mu  <- if (rob) f$mean    else f$mu
        Phi <- if (rob) f$FPC_dis else f$phi
        # Interpolate the FIT onto the truth's grid, never the other way round:
        # the truth is exact everywhere and the fit is not.
        mu_g  <- stats::spline(g, mu, xout = tg)$y
        Phi_g <- vapply(seq_len(M_EVAL), function(j)
          stats::spline(g, Phi[, j], xout = tg)$y, numeric(length(tg)))

        r$mean_err  <- sqrt(sum((mu_g - mu_t[[k]])^2)) / sqrt(sum(mu_t[[k]]^2))
        r$angle_1   <- vec_angle(Phi_g[, 1], Phi_t[[k]][, 1])
        r$angle_2   <- vec_angle(Phi_g[, 2], Phi_t[[k]][, 2])
        # principal_angles() ALREADY RETURNS DEGREES -- it applies 180/pi itself.
        # An earlier version converted a second time and reported subspace angles
        # of 550 to 1450 degrees, which is impossible and is why the guard below
        # exists. Every angle here is in [0, 90] by construction, so anything
        # outside it means a units or ordering error rather than a bad fit.
        r$angle_sub <- max(RKCFC::principal_angles(Phi_g, Phi_t[[k]]))
        stopifnot("an angle outside [0, 90] degrees means a units error" =
                    all(c(r$angle_1, r$angle_2, r$angle_sub) >= 0 &
                        c(r$angle_1, r$angle_2, r$angle_sub) <= 90 + 1e-8))

        xi_hat <- if (rob) f$score else f$xiEst
        xi_tru <- dat$xi[[k]][dat$labels == k, , drop = FALSE]
        if (!is.null(xi_hat) && nrow(xi_hat) == nrow(xi_tru))
          r$cancor <- mean(stats::cancor(xi_hat, xi_tru)$cor)
        r
      }, error = function(e) { r$error <- conditionMessage(e); r })
      kk <- kk + 1L
      rows[[kk]] <- cbind(jb, out, row.names = NULL)
    }
  }
  if (i %% 5L == 0L || i == nrow(mine))
    cat(sprintf("  %d/%d cells, %.1f min elapsed\n", i, nrow(mine),
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

out <- do.call(rbind, rows[seq_len(kk)])

# PROVENANCE. Two task files in one directory can differ in their slicing or in
# the code that produced them, at exactly 100 percent completeness, and nothing
# else here can see it. The combiner refuses to merge files that disagree.
out$n_tasks   <- n_tasks
out$code_hash <- source_hash(file.path("R", "run_estim.R"))
f <- file.path(outdir, sprintf("estim_task%04d.rds", task_id))
saveRDS(out, f)
cat(sprintf("wrote %s  (%d rows, %.1f min)\n", f, nrow(out),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
