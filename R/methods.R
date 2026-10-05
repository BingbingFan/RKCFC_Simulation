# The seven Study A methods, plus the machinery they share.
#
# Every method returns the SAME shape: a partition, a convergence flag, timing,
# and an error string if it failed. Failures are RECORDED, never fatal -- one
# method blowing up on one replicate must not take down an array task, and the
# failure rate is itself a result.
#
# Only M1, M2, M3 and M6 produce a fitted trajectory. M4, M5 and M7 return
# labels and nothing else, so they carry no reconstruction error. That is by
# construction rather than omission; see metrics.R.

source(file.path("R", "generate.R"))
source(file.path("R", "metrics.R"))

# mclust must be ATTACHED, not merely available. Mclust() calls mclustBIC()
# unqualified from its own body, so it resolves it from the search path rather
# than from its namespace; a purely namespaced mclust::Mclust() call fails with
# "could not find function mclustBIC" even though mclust exports it. Passing a
# pre-computed BIC object does not help, because the internal call happens
# anyway. This is a property of mclust, not of this code.
suppressMessages(library(mclust))

# --- shared machinery --------------------------------------------------------

build_optns <- function(cfg, robust) {
  iv <- c(gnum(cfg, "t_min"), gnum(cfg, "t_max"))
  if (robust) {
    list(interval = iv,
         nK       = gint(cfg, "nK_fit"),
         bw       = gnum(cfg, "bw"),
         bwK      = cfg$global$bwK,
         bwmean   = cfg$global$bwmean,
         nRegGrid = gint(cfg, "nRegGrid"),
         basis    = fda::create.bspline.basis(
           iv, nbasis = gint(cfg, "basis_nbasis"),
           norder = gint(cfg, "basis_norder"),
           breaks = seq(iv[1], iv[2], length.out = gint(cfg, "basis_nbreaks"))))
  } else {
    list(dataType = "Sparse", methodMuCovEst = "smooth",
         nRegGrid = gint(cfg, "nRegGrid"),
         methodSelectK = gint(cfg, "nK_fit"),
         useBinnedData = "OFF")
  }
}

# Global FPC scores from ONE fit over all subjects, which the score-then-cluster
# comparators consume.
#
# STATE THIS AS A HANDICAP IN METHODS, because it is one. The design gives the
# two classes DIFFERENT eigenfunctions -- sines against cosines -- so a single
# global basis is wrong for both. M1 and M2 re-estimate class-specific bases
# inside their loop and M4, M5 and M7 cannot. That is exactly the difference the
# ablation is meant to expose, but a reader should be told it is structural
# rather than a deficiency of k-means.
global_scores <- function(dat, cfg, robust, optns) {
  n <- length(dat$Ly)
  fit <- univar_FPCA_k(dat$Ly, dat$Lt, 1L, rep(1L, n), robust, optns)[[1]]
  if (robust) fit$score else fit$xiEst
}

.timed <- function(expr) {
  t0 <- Sys.time()
  val <- tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "mfail"))
  list(value = val, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

.result <- function(labels = NULL, converged = NA, iter = NA_integer_,
                    stop_reason = NA_character_, objective = NA_real_,
                    seconds = NA_real_, error = NA_character_,
                    fits = NULL, robust = NA, optns = NULL,
                    back = identity, fwd = identity) {
  list(labels = labels, converged = converged, iter = iter,
       stop_reason = stop_reason, objective = objective, seconds = seconds,
       error = error, fits = fits, robust = robust, optns = optns,
       back = back, fwd = fwd)
}

# The one place the result table's shape is defined. run_cell.R builds a
# generation-failure row from this so that every row in a task has the same
# columns -- do.call(rbind, ...) on data frames with different column sets is a
# hard error, so a ragged failure row would kill the whole task at the final
# rbind, after all of its compute had been spent. That is the opposite of what
# recording the failure is for.
empty_scores <- function(method = NA_character_, label = NA_character_,
                         error = NA_character_) {
  data.frame(
    method = method, label = label, error = error,
    seconds = NA_real_, converged = NA, iter = NA_integer_,
    stop_reason = NA_character_, objective = NA_real_,
    ccr = NA_real_, ari = NA_real_, nmi = NA_real_,
    rre_in = NA_real_, rre_true = NA_real_, rre_out = NA_real_,
    n_holdout_excluded = NA_integer_, stringsAsFactors = FALSE)
}

# --- the K-centers family: M1, M2, M3 ---------------------------------------

fit_kcfc <- function(dat, cfg, K, robust, transform = "none") {
  optns <- build_optns(cfg, robust)
  Ly <- dat$Ly; back <- identity; fwd <- identity; shift <- NA_real_

  if (transform == "log") {
    # The transform needs positive data and simulated Y is not. c = 1 - min(Y)
    # is recomputed PER REPLICATE (decided 2026-08-30). Checked across arms: c
    # runs from a median of 14.9 on the Gaussian arm to 48.1 at t nu = 1.5, and
    # log(Y + c) spans 3.3 to 4.6 everywhere, so the transform does not
    # degenerate even where the scores have infinite variance.
    #
    # c DOES vary by roughly 1.5x across replicates within an arm, which is a
    # nuisance source of variation M3 carries and no other method does. It is
    # returned so that a wide M3 interquartile range can be attributed rather
    # than guessed at.
    #
    # `fwd` is returned alongside `back` because the scorers need it: the fitted
    # model lives on the log scale, so a subject's scores must be computed from
    # log(Y + c) rather than from Y. `shift` is captured once here and reused
    # for the held-out refit, so the transform is identical in both.
    shift <- 1 - min(unlist(Ly))
    fwd   <- function(y) log(y + shift)
    Ly    <- lapply(Ly, fwd)
    back  <- function(z) exp(z) - shift
  }

  r <- .timed(RKCFC(Ly = Ly, Lt = dat$Lt, K = K, robust = robust, optns = optns,
                    init.method = cfg$global$init_method,
                    maxIter = gint(cfg, "maxIter"),
                    conv.prop = gnum(cfg, "conv_prop"),
                    LOO = glog(cfg, "LOO"),
                    parallel = FALSE, trace = 0))
  if (inherits(r$value, "mfail"))
    return(.result(seconds = r$seconds, error = as.character(r$value)))
  f <- r$value
  .result(labels = f$res_clust$g, converged = f$conv, iter = f$iter,
          stop_reason = f$stop_reason, objective = f$objective_final,
          seconds = r$seconds, fits = f$FPCA_C, robust = robust,
          optns = optns, back = back, fwd = fwd)
}

# --- score-then-cluster: M4, M5 ---------------------------------------------

# The two-stage ablation. M5 uses the robust estimator WITHOUT the K-centers
# iteration, so the difference between M5 and M1 is what the loop buys. The
# point is that with only two eigenfunctions the simpler competitor may be
# enough; if M5 matches M1, the paper should say so plainly.
fit_kmeans_scores <- function(dat, cfg, K, robust) {
  optns <- build_optns(cfg, robust)
  r <- .timed({
    s <- global_scores(dat, cfg, robust, optns)
    stats::kmeans(s, centers = K, nstart = 10)$cluster
  })
  if (inherits(r$value, "mfail"))
    return(.result(seconds = r$seconds, error = as.character(r$value)))
  .result(labels = as.integer(r$value), converged = TRUE, seconds = r$seconds)
}

# --- M6: latent class mixed model -------------------------------------------

fit_hlme <- function(dat, cfg, K) {
  r <- .timed({
    if (!requireNamespace("lcmm", quietly = TRUE)) stop("lcmm is not installed")
    df <- data.frame(id = rep(seq_along(dat$Ly), lengths(dat$Ly)),
                     time = unlist(dat$Lt), y = unlist(dat$Ly))
    f  <- stats::as.formula("y ~ splines::ns(time, 3)")
    # hlme with ng > 1 will not start from nothing: it requires initial values,
    # and the documented way to supply them is a fitted one-class model. Without
    # B it stops with "Please specify initial values", which is a usage error
    # rather than non-convergence and must not be recorded as the latter.
    #
    # NO RANDOM INTERCEPT, and this is not an oversight -- it is the difference
    # between the comparator working and sitting at chance. Measured on this
    # design at d = 8, n = 80, where every other method attains aRI 1.000:
    #
    #   random = ~ 1    conv 1, aRI -0.009, class means 3.91 / 3.63
    #   random = ~ -1   conv 1, aRI +1.000, class means -0.21 / 7.91  (truth)
    #
    # The mechanism: this design has large within-class score variance
    # (lambda_1 = 9, so a subject-level SD near 3), and a single random
    # intercept with SD around 4.5 fits the bimodal marginal well enough that
    # the likelihood prefers one spread-out class over two separated ones. The
    # random effect absorbs exactly the between-class level difference the
    # mixture exists to find. Both fits CONVERGE, so this would not have shown
    # up in a convergence flag.
    #
    # Reporting an established comparator at chance because of a specification
    # choice of ours would be both unfair and easy to attack, so M6 is given the
    # configuration that works. The random-intercept behaviour is a real finding
    # about applying LCMM to sparse data with large within-class variation and is
    # reported in the manuscript; adding it as a second column is cheap if it is
    # needed. Note that without a random effect this is a latent-class
    # spline model rather than a latent-class MIXED model -- label it accordingly.
    m0 <- lcmm::hlme(fixed = f, random = ~ -1, subject = "id", ng = 1,
                     data = df, maxiter = gint(cfg, "maxIter"), verbose = FALSE)
    m  <- lcmm::hlme(fixed = f, mixture = ~ splines::ns(time, 3),
                     random = ~ -1, subject = "id", ng = K, data = df, B = m0,
                     maxiter = gint(cfg, "maxIter"), verbose = FALSE)
    # VALIDATE THE PARTITION BEFORE RETURNING IT. When hlme fails badly it does
    # not error -- it returns an object whose $pprob is short or carries NA
    # classes, having reported "objective function : -1e+09" and RDM = 1. That
    # NA then reaches ccr(), where K <- max(c(truth, est)) is NA and `if (K >
    # 8L)` stops with "missing value where TRUE/FALSE needed" -- three steps from
    # the cause, in a function that has nothing to do with hlme.
    #
    # Measured 2026-09-07: this killed 72 of 250 core tasks, each losing a whole
    # slice of finished compute for every OTHER method too. Stopping here turns a
    # task-killing crash into one recorded failure row for M6 on that cell.
    cl <- m$pprob$class[order(m$pprob$id)]
    if (length(cl) != length(dat$Ly) || anyNA(cl))
      stop(sprintf(paste("hlme returned %d labels for %d subjects with %d missing",
                         "(conv = %s) -- degenerate fit, not a partition"),
                   length(cl), length(dat$Ly), sum(is.na(cl)),
                   paste(m$conv, collapse = "/")), call. = FALSE)
    list(cl = cl, conv = m$conv == 1)
  })
  if (inherits(r$value, "mfail"))
    return(.result(seconds = r$seconds, error = as.character(r$value)))
  .result(labels = as.integer(r$value$cl), converged = r$value$conv,
          seconds = r$seconds)
}

# --- M7: QDA on the estimated scores ----------------------------------------

# THE ASSIGNMENT-RULE ABLATION.
#
# Equation (12) is the squared distance to the class-specific affine subspace:
# it uses the mean and the eigenfunctions and discards the eigenvalues, and it
# carries no prior term. It is therefore NOT the Bayes rule under Equation (10),
# which would be the quadratic discriminant. M7 supplies that comparison.
#
# Implemented in SCORE space, which is what "QDA on the estimated scores" says
# and which is a Gaussian mixture with unconstrained per-class covariance --
# exactly the QDA rule, fitted rather than assumed.
#
# NOTE, because it retires a piece of planned machinery. The plan (7b, option B)
# specified recovering the eigenvalue scale robustly, as (IQR/1.349)^2 of the
# projected scores, because KFPCA returns no eigenvalues and an OBSERVATION-space
# QDA would need Lambda to form Sigma = Phi Lambda Phi' + sigma^2 I. In score
# space no Lambda is needed: mclust estimates the score covariance directly. The
# recovery step is therefore unnecessary here and is not implemented. Record
# this rather than leaving a reader to wonder where it went.
fit_qda_scores <- function(dat, cfg, K, robust) {
  optns <- build_optns(cfg, robust)
  r <- .timed({
    s <- global_scores(dat, cfg, robust, optns)
    # Mclust() resolves mclustBIC from the search path rather than its own
    # namespace, so it needs mclust ATTACHED -- which is why this file calls
    # library(mclust) at the top. Passing a pre-computed BIC object does not
    # avoid that, because the internal unqualified call happens anyway; the BIC
    # is computed here only to pin G and the model name.
    bic <- mclust::mclustBIC(s, G = K, modelNames = "VVV", verbose = FALSE)
    m <- mclust::Mclust(s, x = bic, verbose = FALSE)
    if (is.null(m)) stop("mclust returned no model")
    m$classification
  })
  if (inherits(r$value, "mfail"))
    return(.result(seconds = r$seconds, error = as.character(r$value)))
  .result(labels = as.integer(r$value), converged = TRUE, seconds = r$seconds)
}

# --- dispatch ----------------------------------------------------------------

run_method <- function(method_row, dat, cfg, K) {
  id <- method_row$method_id
  switch(id,
    M1 = fit_kcfc(dat, cfg, K, robust = TRUE),
    M2 = fit_kcfc(dat, cfg, K, robust = FALSE),
    M3 = fit_kcfc(dat, cfg, K, robust = FALSE, transform = "log"),
    M4 = fit_kmeans_scores(dat, cfg, K, robust = FALSE),
    M5 = fit_kmeans_scores(dat, cfg, K, robust = TRUE),
    M6 = fit_hlme(dat, cfg, K),
    M7 = fit_qda_scores(dat, cfg, K, robust = TRUE),
    stop("unknown method: ", id, call. = FALSE))
}

# --- one method on one replicate, scored --------------------------------------

score_method <- function(method_row, dat, cfg, K, seed) {
  # RE-SEED BEFORE EVERY FIT. Without this the RNG state entering a method is
  # whatever the previous method left behind, so a method's result depends on
  # which methods ran before it. Measured on this harness 2026-09-05, arm 4
  # d = 1.5: M2 returned aRI -0.0082 run first, -0.0064 after M1, and +0.0192
  # after M1 and M4 -- same data, same seed, same package.
  #
  # It matters because `init.method = "two-staged"` draws twice (k-means centres,
  # then the random fallback when k-means splits degenerately) and M4/M5 use
  # kmeans(nstart = 10). Two consequences beyond the noise: a recorded
  # (seed_master, scenario_id, rep) was NOT sufficient to reproduce a fit -- you
  # also needed the method list and its order -- and toggling an `include_*`
  # column or adding an estimator-swap row silently moved every method after it.
  # The pilot runs M1 and M2; the core grid runs seven. Without this the two
  # would not agree on M2 even on an identical cell.
  #
  # ONE COMMON SEED, not one per method, and deliberately. The methods are being
  # COMPARED, so a shared starting state pairs them and reduces the variance of
  # the method-to-method differences the paper reports -- common random numbers,
  # standard variance reduction. A bad draw then hits every method together,
  # which is a property of that replicate rather than an artefact.
  #
  # The "|fit" tag keeps the fitting stream distinct from the one that generated
  # the data, so the algorithm's randomness is not drawn from the stream that
  # produced the sample it is fitting.
  # A hand-built `dat` (the failure-path tests) carries neither field; fall back
  # to the caller's seed rather than silently keying every such call the same.
  sid <- if (is.null(dat$scenario_id)) paste0("adhoc", seed) else dat$scenario_id
  rp  <- if (is.null(dat$rep)) 0L else dat$rep
  set.seed(cell_seed(gint(cfg, "seed_master"), paste0(sid, "|fit"), rp))
  res <- run_method(method_row, dat, cfg, K)
  out <- empty_scores(method_row$method_id, method_row$label, res$error)
  out$seconds <- res$seconds; out$converged <- res$converged
  out$iter <- res$iter; out$stop_reason <- res$stop_reason
  out$objective <- res$objective
  if (!is.na(res$error)) return(out)

  # SCORING IS GUARDED TOO, NOT ONLY FITTING. run_method() catches everything the
  # fit throws, so a method that fails honestly yields a failure row and the
  # other six methods on that cell still get scored. A method that "succeeds"
  # while returning a MALFORMED partition used to escape that: the error came out
  # of ccr(), outside any handler, and killed the whole array task -- taking
  # every cell in its slice, and every other method's finished compute with it.
  # Measured 2026-09-07: 72 of 250 core tasks lost this way to M6/hlme.
  #
  # The failure row keeps the success row's columns, so do.call(rbind, ...) over
  # the task's rows still binds; a short row is a hard error that would take down
  # the task all over again.
  sc <- tryCatch({
    out$ccr <- ccr(dat$labels, res$labels)
    out$ari <- mclust::adjustedRandIndex(dat$labels, res$labels)
    out$nmi <- nmi(dat$labels, res$labels)
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(sc)) {
    out$error <- paste("SCORING:", sc)
    return(out)
  }

  # Reconstruction error only for the methods that produce a fitted trajectory.
  if (isTRUE(method_row$reconstructs) && !is.null(res$fits)) {
    e <- reconstruction_errors(dat, res$labels, res$fits, res$robust, res$optns,
                               back = res$back, fwd = res$fwd)
    out$rre_in <- e$rre_in; out$rre_true <- e$rre_true

    # The held-out split is taken on the ORIGINAL scale so that the withheld
    # target stays raw; the retained part is transformed before refitting, with
    # the same shift the full fit used.
    sp <- holdout_split(dat, gnum(cfg, "holdout_frac"),
                        gint(cfg, "holdout_min_ni"), seed)
    out$n_holdout_excluded <- sp$n_excluded
    ho <- tryCatch({
      fits_ho <- univar_FPCA_k(lapply(sp$fit$Ly, res$fwd), sp$fit$Lt, K,
                               res$labels, res$robust, res$optns)
      holdout_error(sp, res$labels, fits_ho, res$robust, res$optns,
                    back = res$back, fwd = res$fwd)
    }, error = function(e) NA_real_)
    out$rre_out <- ho
  }
  out
}
