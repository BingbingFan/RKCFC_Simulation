# Settings studies -- block S.
#
#   Rscript R/run_settings.R <row> <task_id> [n_tasks] [outdir] [n_reps] [n_over]
#
#   row   init | components | smoothing | estimators
#
# Two of the block's six rows are FREE and have no runner: convergence/cycling
# and computational cost are read off the core grid, because RKCFC() already
# records `objective`, `stop_reason`, `cycle_length` and wall clock on every fit.
# Do not re-run anything to obtain them.
#
# THE `init` ROW IS THE ONE THAT MEASURES A DEFECT, and it is the reason
# RKCFC() gained an `init.labels` argument. Scheme `pace_kmeans` starts from
# k-means on NON-ROBUST PACE scores while fitting robustly -- which is exactly
# what the released 0.1.x package did regardless of the `robust` flag.
# Running it deliberately turns the correction from an assertion into a
# measurement.
#
# The `oracle` scheme bounds what any initialiser could buy: it starts from the
# true partition, so the gap between it and the best real scheme is the whole of
# what initialisation costs. If oracle and two-staged agree, initialisation is
# not the story on that arm and the 0.1.x start's effect is bounded accordingly.

source(file.path("R", "methods.R"))

args    <- commandArgs(trailingOnly = TRUE)
row_id  <- if (length(args) >= 1) args[1] else "init"
task_id <- if (length(args) >= 2) as.integer(args[2]) else 1L
n_tasks <- if (length(args) >= 3) as.integer(args[3]) else 1L
outdir  <- if (length(args) >= 4) args[4] else file.path("results", paste0("S_", row_id))
n_reps  <- if (length(args) >= 5) as.integer(args[5]) else NA_integer_
n_over  <- if (length(args) >= 6) as.integer(args[6]) else NA_integer_

cfg <- load_config()
K   <- gint(cfg, "K_true")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
reps <- seq_len(if (is.na(n_reps)) gint(cfg, "n_reps") else n_reps)

# Which arms each row runs on, from the design table. Arms 1 and 4 for the rows
# about behaviour under tails; arm 2 alone for the components row, which is
# about retention rather than robustness.
# SMOOTHING RUNS FOUR ARMS AS OF 2026-09-07, and the two additions are a matched
# pair rather than more coverage. `mean_robust`, the robust mean option, was measured in
# the estimation block as BETTER on arm 4 and WORSE on arm 3 -- and arm 3 is the
# only arm with nonzero skewness, which is exactly what a Huber M-estimator of
# location is not consistent for. Whether that mean penalty reaches the PARTITION
# could not be inferred from `estim`, because mean_err is a norm and a norm is
# sign-blind: it cannot separate "both class means shifted the same way", which
# the partition cannot see, from "they shifted apart", which it can.
#
# Arm 2 is here to make arm 3 readable, not for its own sake. Fleishman(0, 4) and
# Fleishman(1.5, 4) share a kurtosis and differ only in skewness, so the pair
# isolates skewness the way B10/B12 isolates the observation floor. Arm 3 alone,
# against the Gaussian and t(1.5) arms already present, would confound skewness
# with excess kurtosis and could not settle anything.
ARMS <- list(init = c(1L, 4L), components = 2L, smoothing = c(1L, 2L, 3L, 4L),
             estimators = c(1L, 2L, 4L))
if (!row_id %in% names(ARMS))
  stop("unknown row: ", row_id, ". Use ", paste(names(ARMS), collapse = ", "),
       call. = FALSE)

# The levels of each row. Everything not yet buildable is listed anyway and
# recorded as a stated NOT IMPLEMENTED, so the gap is visible in the combined
# table and the completeness arithmetic stays right, rather than the level being
# silently absent.
LEVELS <- list(
  # two_staged IS THE CONFIGURED DEFAULT AND MUST BE A LEVEL. Without it the row
  # compares five alternatives to each other and never to what the manuscript
  # actually uses, so the baseline has to be borrowed from another block -- which
  # is seeded differently, making the comparison UNPAIRED on exactly the contrast
  # the row exists to measure. The decisions block ran without it on 2026-09-06
  # and the omission cost a defensible answer to "is the default the right one".
  #
  # It needs no branch below: any level that matches none of them leaves init_m
  # at cfg$global$init_method and init_l at NULL, which IS two-staged. Adding a
  # branch for it would be a second definition of the default, free to drift.
  init = data.frame(
    level = c("two_staged", "kendall_kmeans", "pace_kmeans", "random20",
              "ward", "oracle"),
    note = c("the CONFIGURED default: k-means with a random fallback when a class falls below 10 curves",
             "k-means on Kendall scores -- the intended robust start, fallback REMOVED",
             "k-means on PACE scores while fitting robustly -- the 0.1.x start",
             "20 random restarts, best by the objective J",
             "Ward on the global scores",
             "the TRUE partition; bounds what any initialiser could buy"),
    stringsAsFactors = FALSE),
  components = data.frame(
    level = c("M1", "M2", "M3", "M4", "fve80", "fve90", "fve95", "fve99"),
    note = "retained components, fixed or chosen by an FVE threshold",
    stringsAsFactors = FALSE),
  smoothing = data.frame(
    level = c("bw_half", "bw_base", "bw_double", "mean_robust"),
    note = c("half the configured bandwidth", "as configured",
             "twice the configured bandwidth",
             "robust local-linear M-estimator mean -- PACKAGE ITEM 1"),
    stringsAsFactors = FALSE),
  estimators = data.frame(
    level = c("E1_kendall", "E2_sign", "E3_rank"),
    note = c("Kendall's tau -- the proposal",
             "spherical / spatial-sign -- PACKAGE ITEM 2",
             "rank + PACE -- PACKAGE ITEM 3"),
    stringsAsFactors = FALSE))

# mean_robust LEFT THIS LIST on 2026-09-07: the robust mean option is built. KFPCA()
# gained `meanrobust`, MeanEst() gained a Huber IRLS branch, and predict.KFPCA()
# inherits the choice off the fitted object rather than re-deriving it. The
# package's own regression suite asserts the default is byte-for-byte the old
# least-squares call, so every result banked before today still reproduces.
# EMPTY AS OF 2026-09-07: items 2 and 3 are built too, so no level in this block
# is a stated gap any more. Kept as a named vector rather than deleted, because
# the machinery below is the right home for the next unbuilt level and rebuilding
# it from memory is how a gap becomes a silent omission instead of a printed one.
NOT_BUILT <- character(0)

lev  <- LEVELS[[row_id]]
jobs <- expand.grid(arm = ARMS[[row_id]], rep = reps, KEEP.OUT.ATTRS = FALSE)
jobs$d <- gnum(cfg, "d_hard")
jobs$scenario_id <- sprintf("S%s_a%d", toupper(substr(row_id, 1, 3)), jobs$arm)
jobs <- jobs[order(jobs$rep, jobs$arm), ]
mine <- jobs[(seq_len(nrow(jobs)) - 1L) %% n_tasks + 1L == task_id, , drop = FALSE]

cat(sprintf("block S/%s | task %d of %d | %d of %d cells | %d levels\n",
            row_id, task_id, n_tasks, nrow(mine), nrow(jobs), nrow(lev)))
if (!nrow(mine))
  stop(sprintf(paste("task %d of %d was assigned NO cells (the row has %d).",
                     "Pass the ORIGINAL array size as n_tasks when re-running."),
               task_id, n_tasks, nrow(jobs)), call. = FALSE)

empty_srow <- function(level = NA_character_, error = NA_character_) {
  data.frame(row = row_id, level = level, error = error,
             seconds = NA_real_, converged = NA, iter = NA_integer_,
             stop_reason = NA_character_, objective = NA_real_,
             ari = NA_real_, ccr = NA_real_, nmi = NA_real_,
             rre_in = NA_real_, rre_true = NA_real_,
             init_method = NA_character_,
             # restarts only: the share of starts reaching the best partition.
             # This is the direct answer to "how much does initialisation
             # matter" -- a value near 1 says the objective surface has one
             # basin and the 0.1.x start cannot have mattered much on that cell.
             prop_best = NA_real_, n_start_failed = NA_integer_,
             stringsAsFactors = FALSE)
}

# Ward on the global scores, as a starting partition.
ward_start <- function(sc, K) {
  h <- stats::hclust(stats::dist(sc), method = "ward.D2")
  as.integer(stats::cutree(h, k = K))
}

rows <- vector("list", nrow(mine) * nrow(lev))
kk <- 0L
t_start <- Sys.time()

for (i in seq_len(nrow(mine))) {
  jb  <- mine[i, ]
  arm <- cfg$arms[cfg$arms$arm == jb$arm, ]
  dat <- tryCatch(
    generate_replicate(cfg, arm, d = jb$d, rep = jb$rep,
                       scenario_id = jb$scenario_id,
                       n = if (is.na(n_over)) NULL else n_over),
    error = function(e) structure(conditionMessage(e), class = "genfail"))

  if (inherits(dat, "genfail")) {
    for (j in seq_len(nrow(lev))) {
      kk <- kk + 1L
      rows[[kk]] <- cbind(jb, empty_srow(lev$level[j],
                            paste("GENERATION:", as.character(dat))),
                          row.names = NULL)
    }
    next
  }

  for (j in seq_len(nrow(lev))) {
    lv <- lev$level[j]
    r  <- empty_srow(lv)

    if (lv %in% names(NOT_BUILT)) {
      r$error <- paste("NOT IMPLEMENTED:", NOT_BUILT[[lv]])
      kk <- kk + 1L; rows[[kk]] <- cbind(jb, r, row.names = NULL); next
    }

    # Every level starts from the SAME RNG state, so a difference between levels
    # is the level and never the draw.
    set.seed(cell_seed(gint(cfg, "seed_master"),
                       paste0(jb$scenario_id, "|fit"), jb$rep))

    robust <- TRUE                       # the proposal, unless a level says else
    optns  <- build_optns(cfg, robust)
    init_m <- cfg$global$init_method
    init_l <- NULL
    use_restarts <- FALSE
    # The measurements handed to the fit. Only E3_rank changes them; every other
    # level fits the data as generated.
    Ly_fit <- dat$Ly

    if (row_id == "init") {
      if (lv == "kendall_kmeans") { init_m <- "kmeans" }
      else if (lv == "pace_kmeans") {
        # THE 0.1.x START: a start from NON-ROBUST scores, then a robust fit.
        sc0 <- tryCatch(global_scores(dat, cfg, FALSE, build_optns(cfg, FALSE)),
                        error = function(e) NULL)
        if (is.null(sc0)) { r$error <- "PACE scores failed"; kk <- kk + 1L
          rows[[kk]] <- cbind(jb, r, row.names = NULL); next }
        init_l <- stats::kmeans(sc0, centers = K, nstart = 10)$cluster
      }
      else if (lv == "random20") { init_m <- "random"; use_restarts <- TRUE }
      else if (lv == "ward") {
        sc0 <- tryCatch(global_scores(dat, cfg, robust, optns),
                        error = function(e) NULL)
        if (is.null(sc0)) { r$error <- "scores failed"; kk <- kk + 1L
          rows[[kk]] <- cbind(jb, r, row.names = NULL); next }
        init_l <- ward_start(sc0, K)
      }
      else if (lv == "oracle") init_l <- dat$labels
    } else if (row_id == "components") {
      optns$nK <- if (grepl("^M", lv)) as.integer(sub("^M", "", lv)) else optns$nK
      # FVE levels choose nK from the data, per the rule settled 2026-09-05.
      if (grepl("^fve", lv)) {
        thr <- as.numeric(sub("^fve", "", lv)) / 100
        kf <- tryCatch({
          o <- optns; o$nK <- 6L
          f <- KFPCA(Lt = dat$Lt, Ly = dat$Ly, interval = o$interval, nK = 6L,
                     bw = o$bw, bwK = o$bwK, nRegGrid = o$nRegGrid,
                     fdParobj = o$basis)
          v <- pmax(f$values, 0)          # the surface is not PSD; see KFPCA docs
          which(cumsum(v) / sum(v) >= thr)[1]
        }, error = function(e) NA_integer_)
        if (is.na(kf)) { r$error <- "FVE selection failed"; kk <- kk + 1L
          rows[[kk]] <- cbind(jb, r, row.names = NULL); next }
        optns$nK <- max(1L, min(6L, kf))
        r$level <- sprintf("%s(M=%d)", lv, optns$nK)
      }
    } else if (row_id == "smoothing") {
      if (lv == "bw_half")   optns$bw <- optns$bw / 2
      if (lv == "bw_double") optns$bw <- optns$bw * 2
      # PACKAGE ITEM 1. The BANDWIDTH is untouched here on purpose: this level
      # swaps the mean ESTIMATOR at the configured bandwidth, so it sits beside
      # bw_base as a like-for-like contrast rather than confounding the two
      # things this row varies.
      if (lv == "mean_robust") optns$meanrobust <- TRUE
    } else if (row_id == "estimators") {
      # ALL THREE ARE BUILDABLE as of 2026-09-07.
      #
      # E1_kendall  the proposal, and the defaults already set it.
      # E2_sign     spherical / spatial-sign scatter, at the SAME mean, so the
      #             contrast against E1 is the scatter estimator alone.
      # E3_rank     normal-scores marginal transform, then ordinary PACE. It is a
      #             TRANSFORM OF THE DATA, not an option, so it switches the
      #             engine to the non-robust branch and needs that branch's optns
      #             -- passing KFPCA's option names to fdapace::FPCA() is rejected
      #             outright ("Invalid option names").
      #
      # THIS ROW IS WHERE E3 IS JUDGED FAIRLY. A monotone marginal transform
      # changes the estimand, so E3's eigenfunctions and reconstruction are not
      # on the data's scale -- but a PARTITION is invariant to it, and a partition
      # is what this row scores.
      if (lv == "E2_sign") optns$scatter <- "sign"
      if (lv == "E3_rank") {
        robust <- FALSE
        optns  <- build_optns(cfg, FALSE)
        Ly_fit <- normal_scores(dat$Ly)
      }
    }

    t0 <- Sys.time()
    rs <- NULL
    f <- tryCatch({
      if (use_restarts) {
        # The component is $fit, NOT $best -- an earlier version reached for
        # $best, got NULL, and died downstream in adjustedRandIndex() with
        # "arguments must be vectors of the same length", which names the
        # symptom three steps from the cause.
        rs <- RKCFC_restarts(Ly = Ly_fit, Lt = dat$Lt, K = K, robust = robust,
                       optns = optns, nstart = 20L, init.method = "random",
                       maxIter = gint(cfg, "maxIter"),
                       conv.prop = gnum(cfg, "conv_prop"),
                       LOO = glog(cfg, "LOO"), parallel = FALSE, trace = 0)
        rs$fit
      } else {
        RKCFC(Ly = Ly_fit, Lt = dat$Lt, K = K, robust = robust, optns = optns,
              init.method = init_m, init.labels = init_l,
              maxIter = gint(cfg, "maxIter"), conv.prop = gnum(cfg, "conv_prop"),
              LOO = glog(cfg, "LOO"), parallel = FALSE, trace = 0)
      }
    }, error = function(e) structure(conditionMessage(e), class = "sfail"))
    r$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    if (inherits(f, "sfail")) {
      r$error <- as.character(f)
      kk <- kk + 1L; rows[[kk]] <- cbind(jb, r, row.names = NULL); next
    }
    g <- f$res_clust$g
    stopifnot("the fit returned no partition" = length(g) == length(dat$labels))
    r$converged <- f$conv; r$iter <- f$iter; r$stop_reason <- f$stop_reason
    r$objective <- f$objective_final; r$init_method <- f$init.method
    if (!is.null(rs)) {
      r$prop_best <- rs$prop_best; r$n_start_failed <- rs$n_failed
    }
    r$ari <- mclust::adjustedRandIndex(dat$labels, g)
    r$ccr <- ccr(dat$labels, g); r$nmi <- nmi(dat$labels, g)
    # RECONSTRUCTION IS LEFT NA FOR E3_rank, DELIBERATELY. The fit lives on the
    # normal-scores scale and `dat$Ly` does not, so `reconstruction_errors()`
    # would compare two different scales and return a plausible number for a
    # meaningless comparison -- the same defect class as M3's exponentiated-scale
    # error, which passed every test until it was measured at 3.9e13. The
    # transform has no elementwise inverse, so there is no `back` to supply.
    # The PARTITION is scale-invariant and is scored above; that is what this row
    # compares.
    e <- if (identical(lv, "E3_rank")) list(rre_in = NA_real_, rre_true = NA_real_)
         else tryCatch(reconstruction_errors(dat, g, f$FPCA_C, robust, optns),
                       error = function(err) list(rre_in = NA_real_, rre_true = NA_real_))
    r$rre_in <- e$rre_in; r$rre_true <- e$rre_true
    kk <- kk + 1L; rows[[kk]] <- cbind(jb, r, row.names = NULL)
  }
  if (i %% 2L == 0L || i == nrow(mine))
    cat(sprintf("  %d/%d cells, %.1f min elapsed\n", i, nrow(mine),
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}

out <- do.call(rbind, rows[seq_len(kk)])

# PROVENANCE. Two task files in one directory can differ in their slicing or in
# the code that produced them, at exactly 100 percent completeness, and nothing
# else here can see it. The combiner refuses to merge files that disagree.
out$n_tasks   <- n_tasks
out$code_hash <- source_hash(file.path("R", "run_settings.R"))
f <- file.path(outdir, sprintf("S_%s_task%04d.rds", row_id, task_id))
saveRDS(out, f)
cat(sprintf("wrote %s  (%d rows, %.1f min)\n", f, nrow(out),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
