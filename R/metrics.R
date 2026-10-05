# Metrics, in the two tiers the design specifies.
#
#   ALL SEVEN METHODS   classification accuracy: CCR, aRI, NMI. These need only
#                       labels, which every method returns.
#
#   M1, M2, M3 ONLY     reconstruction error, three ways. Only these three are
#                       FPCA-based and directly comparable. M4, M5 and M7 return
#                       labels and nothing else, and M6 fits a different model
#                       family (parametric splines with a random intercept), so
#                       scoring it alongside would confound partition quality
#                       with model family.
#
# WHY NO COMMON RECONSTRUCTION IS IMPOSED. Scoring every partition with one
# shared estimator sounds neutral and is not. RKCFC selects its partition by
# minimising reconstruction error under Kendall FPCA, so scoring it with Kendall
# FPCA scores it on its own training objective -- and PACE does exactly the same
# favour to KCFC. The two estimators are biased toward opposite halves of the
# method list, and PACE additionally has no population target at all on the
# infinite-variance arm. Reporting each method's own fit dissolves the question
# instead of adjudicating it.

# ---------------------------------------------------------------------------
# Correct classification rate, with labels aligned by the permutation that
# maximises agreement.
#
# R1 minor comment 10 is right that this is required: aRI is permutation
# invariant but CCR is not, so without alignment CCR would depend on which
# arbitrary integer the algorithm happened to assign to each class. Exhaustive
# over permutations rather than Hungarian, because K <= 6 here (720 cases) and
# an exact enumeration needs no dependency and cannot be subtly wrong.
ccr <- function(truth, est) {
  # REJECT A MALFORMED PARTITION BY NAME. With an NA in `est`, K below is NA and
  # `if (K > 8L)` stops with "missing value where TRUE/FALSE needed" -- an error
  # that names neither the method nor the real problem. hlme produces exactly
  # this on a degenerate fit; it cost 72 core tasks on 2026-09-07.
  if (length(truth) != length(est))
    stop(sprintf("ccr(): %d true labels against %d estimated", length(truth),
                 length(est)), call. = FALSE)
  if (anyNA(est) || anyNA(truth))
    stop(sprintf("ccr(): partition has %d missing labels of %d", sum(is.na(est)),
                 length(est)), call. = FALSE)
  truth <- as.integer(as.factor(truth)); est <- as.integer(as.factor(est))
  K <- max(c(truth, est))
  if (K > 8L) stop("exhaustive label matching is only sane up to K = 8.", call. = FALSE)
  best <- 0
  grid <- .permutations(K)
  for (r in seq_len(nrow(grid))) {
    best <- max(best, mean(grid[r, ][est] == truth))
  }
  best
}

.permutations <- function(n) {
  if (n == 1L) return(matrix(1L, 1, 1))
  sub <- .permutations(n - 1L)
  do.call(rbind, lapply(seq_len(n), function(i)
    cbind(i, matrix(seq_len(n)[-i][sub], nrow(sub)))))
}

# ---------------------------------------------------------------------------
# Normalised mutual information, symmetric ("sum") normalisation:
#   NMI = 2 I(A;B) / (H(A) + H(B)),  which is 1 for identical partitions and
#   0 when they are independent. State the normalisation in Methods: the
#   "sqrt" and "max" variants give different numbers on the same partitions.
nmi <- function(truth, est) {
  tab <- table(truth, est); n <- sum(tab)
  p   <- tab / n; pa <- rowSums(p); pb <- colSums(p)
  H   <- function(x) -sum(x[x > 0] * log(x[x > 0]))
  idx <- which(p > 0, arr.ind = TRUE)
  I   <- sum(p[idx] * log(p[idx] / (pa[idx[, 1]] * pb[idx[, 2]])))
  den <- H(pa) + H(pb)
  if (den <= 0) return(NA_real_)
  2 * I / den
}

# ---------------------------------------------------------------------------
# Relative reconstruction error: residual sum of squares over total sum of
# squares about the mean. 0 is perfect, 1 is no better than the sample mean.
# Unit-free, so the simulation and the application are directly comparable --
# which is what repairs the five-order-of-magnitude scale gap between the
# submitted Table 3 and Section 3.2.
rre_vec <- function(y, yhat) sum((y - yhat)^2) / sum((y - mean(y))^2)

# ---------------------------------------------------------------------------
# Split each subject's observations into a fitting part and a held-out part.
#
# Subjects with fewer than `min_ni` observations are excluded from the held-out
# metric ONLY: they cannot spare 20 percent and still be fitted. Withholding a
# fixed single observation instead would make the withheld fraction vary from
# 50 percent to 20 percent across subjects and turn the metric into a mixture of
# different quantities. Those subjects still contribute to fitting and to every
# other metric, and the number excluded is reported.
holdout_split <- function(dat, frac, min_ni, seed) {
  set.seed(seed)
  n <- length(dat$Ly)
  fit_y <- fit_t <- out_y <- out_t <- vector("list", n)
  scored <- logical(n)
  for (i in seq_len(n)) {
    ni <- length(dat$Ly[[i]])
    if (ni < min_ni) {
      fit_y[[i]] <- dat$Ly[[i]]; fit_t[[i]] <- dat$Lt[[i]]
      next
    }
    k <- max(1L, round(frac * ni))
    if (ni - k < 3L) k <- ni - 3L                   # a 2-component fit needs 3
    if (k < 1L) { fit_y[[i]] <- dat$Ly[[i]]; fit_t[[i]] <- dat$Lt[[i]]; next }
    hold <- sample.int(ni, k)
    out_y[[i]] <- dat$Ly[[i]][hold];  out_t[[i]] <- dat$Lt[[i]][hold]
    fit_y[[i]] <- dat$Ly[[i]][-hold]; fit_t[[i]] <- dat$Lt[[i]][-hold]
    scored[i] <- TRUE
  }
  list(fit = list(Ly = fit_y, Lt = fit_t),
       out = list(Ly = out_y, Lt = out_t),
       scored = scored, n_excluded = sum(!scored))
}

# ---------------------------------------------------------------------------
# Evaluate a fitted class model on an arbitrary grid.
#
# univar_predicted_yhat() evaluates only at the subject's own observation times,
# because it passes newLt and newLy together. For the true-curve error the fit
# is needed on the dense grid instead, so the scores are taken from that call
# and the mean and eigenfunctions are interpolated to the grid separately.
curve_on_grid <- function(fit, y_i, t_i, tgrid, robust, optns) {
  pr <- univar_predicted_yhat(fit, y_i, t_i, robust, optns)
  xi <- as.vector(pr$xi_Est_i)
  if (robust) { g <- fit$RegGrid; mu <- fit$mean; Phi <- fit$FPC_dis }
  else        { g <- fit$workGrid; mu <- fit$mu;  Phi <- fit$phi }
  M   <- length(xi)
  mu_g  <- stats::spline(g, mu, xout = tgrid)$y
  Phi_g <- vapply(seq_len(M), function(m)
    stats::spline(g, Phi[, m], xout = tgrid)$y, numeric(length(tgrid)))
  as.vector(mu_g + Phi_g %*% xi)
}

# ---------------------------------------------------------------------------
# The three reconstruction errors, for a method that produces a fitted model.
#
# A TRANSFORMED METHOD NEEDS BOTH DIRECTIONS, AND ONLY ONE OF THEM ON EACH SIDE.
#
# M3 fits on log(Y + c). Scoring it therefore has three separate obligations,
# and getting any one of them wrong is silent:
#
#   `fwd`   the model was estimated on the log scale, so the data handed to
#           univar_predicted_yhat() to compute a subject's scores must be on
#           that scale too. Passing raw Y projects raw data onto a log-scale
#           mean and eigenfunctions, which is meaningless but throws no error.
#   `back`  the FITTED values come back on the log scale and must be returned
#           to the original one.
#   neither the TARGET is the observed data as it stands. Applying `back` to it
#           exponentiates data that was never transformed.
#
# An earlier version applied `back` to both sides and passed raw Y to the model.
# Measured on the wiring-test replicate, that put M3's target on a scale running
# to 2.5e7 and reported rre_true = 3.9e13 against M2's 0.034. Corrected, M3
# reads 0.069 against M2's 0.015 -- the same order, as it must be once all
# comparative error metrics are on exactly the same scale.
reconstruction_errors <- function(dat, labels, fits, robust, optns,
                                  back = identity, fwd = identity) {
  K <- length(fits)
  ids <- split(seq_along(labels), labels)

  yhat_in <- numeric(0); y_in <- numeric(0)
  for (k in seq_len(K)) for (i in ids[[as.character(k)]]) {
    p <- univar_predicted_yhat(fits[[k]], fwd(dat$Ly[[i]]), dat$Lt[[i]],
                               robust, optns)
    yhat_in <- c(yhat_in, back(p$fittedX)); y_in <- c(y_in, dat$Ly[[i]])
  }
  in_sample <- rre_vec(y_in, yhat_in)

  # true-curve: no estimated quantity on the target side, so no choice of
  # estimator can tilt it. X_true is already on the original scale, so it takes
  # neither transform. Simulation only.
  true_curve <- NA_real_
  if (!is.null(dat$X_true)) {
    xt <- numeric(0); xh <- numeric(0)
    for (k in seq_len(K)) for (i in ids[[as.character(k)]]) {
      xh <- c(xh, back(curve_on_grid(fits[[k]], fwd(dat$Ly[[i]]), dat$Lt[[i]],
                                     dat$tgrid, robust, optns)))
      xt <- c(xt, dat$X_true[i, ])
    }
    true_curve <- rre_vec(xt, xh)
  }

  list(rre_in = in_sample, rre_true = true_curve)
}

# Held-out error needs a refit on the retained data, so it is computed by the
# caller, which owns the fitting engine. This helper only scores.
holdout_error <- function(split, labels, fits, robust, optns,
                          back = identity, fwd = identity) {
  K <- length(fits); y <- numeric(0); yh <- numeric(0)
  for (i in which(split$scored)) {
    k <- labels[i]
    p <- univar_predicted_yhat(fits[[k]], fwd(split$fit$Ly[[i]]),
                               split$fit$Lt[[i]], robust, optns)
    # scores from the retained points, evaluated at the withheld times
    xi <- as.vector(p$xi_Est_i)
    if (robust) { g <- fits[[k]]$RegGrid; mu <- fits[[k]]$mean; Phi <- fits[[k]]$FPC_dis }
    else        { g <- fits[[k]]$workGrid; mu <- fits[[k]]$mu;  Phi <- fits[[k]]$phi }
    tt <- split$out$Lt[[i]]
    mu_t  <- stats::spline(g, mu, xout = tt)$y
    Phi_t <- vapply(seq_along(xi), function(m)
      stats::spline(g, Phi[, m], xout = tt)$y, numeric(length(tt)))
    yh <- c(yh, back(as.vector(mu_t + Phi_t %*% xi)))
    y  <- c(y,  split$out$Ly[[i]])          # the target is already raw
  }
  if (!length(y)) return(NA_real_)
  rre_vec(y, yh)
}

# ---------------------------------------------------------------------------
# Monte Carlo standard error of a median, by NONPARAMETRIC BOOTSTRAP.
#
# Every median carries its Monte Carlo uncertainty; reporting a median without one invites
# the question of whether a difference is real.
#
# WHY NOT THE TEXTBOOK FORMULA. `sd(x) * sqrt(pi/2) / sqrt(B)` is the
# large-sample SE of a median UNDER NORMALITY: the general form is
# 1 / (2 * f(m) * sqrt(B)), and `sd * sqrt(pi/2)` is what 1 / (2 f(m)) becomes
# for a normal f. Both of this study's interesting cells break that assumption,
# in OPPOSITE directions.
#
# Measured 2026-09-05 by brute force -- 4,000 repetitions of a B = 100 study, so
# the SD of those 4,000 medians IS the true MCSE by definition:
#
#   replicate distribution              TRUE     normal   bootstrap   normal err
#   unimodal, interior (arms 1-3)     0.00609   0.00619    0.00720        +2%
#   ATOM at the ceiling (CCR = 1)     0.00128   0.00562    0.00139      +340%
#   BIMODAL (arm 4, t nu = 1.5)       0.27849   0.04435    0.27610       -84%
#
# The ceiling case errs LARGE, which merely wastes power. The bimodal case errs
# SMALL, which is the dangerous direction: applying the "a difference must exceed
# twice the MCSE" rule, the true threshold on that cell is 0.556 and the normal
# form reports 0.088 -- six times too permissive. A gap of 0.2 in aRI, which is
# pure noise there, would be declared a real difference. Arm 4 is the
# infinite-variance arm the method exists for.
#
# The mechanism: on a bimodal cell the median sits in the empty valley BETWEEN
# the modes, so f(m) is near zero and the true wobble is large -- a few
# replicates changing basin swing the median across the gap. `sd` meanwhile is
# modest because the values are bunched at two poles, so the formula reads the
# bunching as stability. The bootstrap resamples the actual distribution and
# tracks the truth in all three cases.
#
# SEEDED, because a summary table must not change between two runs of
# combine.R on identical input. The RNG state is saved and restored so that
# calling this never perturbs a caller's stream.
mcse_median <- function(x, n_boot = 2000L, seed = 20260905L) {
  x <- x[is.finite(x)]; B <- length(x)
  if (B < 5L) return(NA_real_)
  if (diff(range(x)) == 0) return(0)          # a constant column has no wobble

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  }
  set.seed(seed)
  stats::sd(vapply(seq_len(n_boot),
                   function(i) stats::median(x[sample.int(B, B, replace = TRUE)]),
                   numeric(1)))
}

# --- PROVENANCE: what the existing guards cannot see --------------------------
#
# Two hazards have already occurred, three days apart, and NEITHER is visible to
# any check that was in place. Both leave a results directory whose task files
# were computed under DIFFERENT ASSUMPTIONS, at exactly 100 percent completeness:
#
#   2026-09-07  `estim` was submitted at NTASKS = 25 and again at 100. Task files
#               are named by task id, so the second run OVERWROTE the first by
#               name -- but the slicing differs, so the two runs assign different
#               cells to the same file. Harmless only because timestamps proved
#               the second overwrote every file of the first.
#   2026-09-08  `run_settings.R` was uploaded while job 3370777 was reading it.
#               `Rscript file.R` parses as it reads, so task 164 died on a
#               truncated source. That failure was LOUD. The silent version is a
#               task that starts AFTER the upload finishes and quietly runs
#               different code from its siblings in the same array.
#
# The completeness check passes in both cases (the row count is right), the
# stale-file check passes (no surplus files), and the failure table passes (no
# errors). Every signal reads clean. Hence a stamp.

# A short digest of the code that can change a result. DELIBERATELY NARROW:
# combiners, table scripts, tests and diagnostics are excluded because they
# cannot affect a run, and including them would trip this guard every time a
# table script is edited while a block is in flight. A guard that cries wolf gets
# ignored, which is worse than no guard.
source_hash <- function(extra = character(0)) {
  fs <- c(file.path("R", c("methods.R", "generate.R", "metrics.R", "jobs.R",
                           "separation.R")),
          list.files("config", pattern = "\\.csv$", full.names = TRUE),
          extra)
  fs <- sort(unique(fs[file.exists(fs)]))
  if (!length(fs)) return(NA_character_)
  h  <- tools::md5sum(fs)              # base R; no extra dependency on the cluster
  tf <- tempfile(); on.exit(unlink(tf), add = TRUE)
  writeLines(paste(basename(names(h)), h), tf)
  substr(unname(tools::md5sum(tf)), 1L, 12L)
}

# Read a block's task files and REFUSE to combine ones that do not belong
# together. Replaces a bare do.call(rbind, lapply(files, readRDS)), which either
# fails with "numbers of columns of arguments do not match" -- true but
# unhelpful -- or, worse, succeeds on files that should never have been mixed.
read_task_files <- function(files) {
  if (!length(files)) stop("no task files to read.", call. = FALSE)
  xs <- lapply(files, readRDS)

  cols <- lapply(xs, function(z) sort(names(z)))
  if (length(unique(cols)) > 1L) {
    k <- which(!vapply(cols, identical, logical(1), cols[[1]]))
    stop(sprintf(paste("task files disagree on their COLUMNS, so they were written by",
                       "different versions of the code.\n  %d of %d differ, first: %s\n",
                       " Clear this directory and re-run the block."),
                 length(k), length(files), basename(files[k[1]])), call. = FALSE)
  }

  raw <- do.call(rbind, xs)

  # Present only on files written after 2026-09-08; older banked results predate
  # the stamp and are checked on columns alone, which is all that can be asked of
  # them.
  for (fld in c("n_tasks", "code_hash")) {
    if (!fld %in% names(raw)) next
    u <- unique(raw[[fld]])
    if (length(u) > 1L) {
      per <- vapply(xs, function(z) as.character(z[[fld]][1]), character(1))
      tb  <- table(per)
      stop(sprintf(paste("task files disagree on `%s`: %s.\n",
                         " This directory mixes %d runs made under different",
                         "assumptions, and the row count cannot see it.\n",
                         " Counts: %s\n  Clear the directory and re-run the block."),
                   fld, paste(u, collapse = " vs "), length(tb),
                   paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = ", ")),
           call. = FALSE)
    }
  }
  raw
}
