# STUDY F4 -- CLUSTERING STABILITY.
# 36-MONTH VERSION (2026-10-02) of run_realdata_f4.R: the same work list, with the
# window, bandwidths, primary scale and results folder from realdata36_settings.R.
# Combine with the unchanged R/combine_realdata_f4.R <this F4 dir> <36-month F1 dir>.
#
#   Rscript R/run_realdata36_f4.R <task_id> [n_tasks] [outdir] [n_start] [n_sub]
#
# WHAT F4 CAN AND CANNOT ASK. There are no true labels on this cohort, so no
# quantity here measures ACCURACY. That question belongs to the simulation, where
# truth exists. F4 asks the three questions that are answerable without labels:
#
#   1. IS THE PARTITION REPRODUCIBLE?  100 random initialisations, compared to
#      EACH OTHER by adjusted Rand index -- never against a truth. If different
#      starts disagree, the reported partition is one draw from a lottery rather
#      than a finding. The simulation says this is a live risk: `prop_best` was
#      0.25 on the heavy-tailed arm, so three starts in four missed the best
#      objective.
#   2. DOES IT SURVIVE THE SAMPLE?  500 subsamples at 80 percent, scored by
#      per-cluster Jaccard against the full-data partition. A cluster that
#      dissolves under resampling is a property of these particular patients.
#   3. HOW CONFIDENT IS EACH ASSIGNMENT?  The margin between the best and
#      second-best class reconstruction error, per subject. A partition can be
#      perfectly reproducible and still have most subjects on the boundary, which
#      neither of the first two questions can see.
#
# Plus the K selection rerun with the CORRECTED silhouette, because
# the submitted K = 2 rests on a coefficient computed under different smoothing
# from the model it scored.
#
# THE MARGIN NEEDS NO NEW PACKAGE CODE. RKCFC() returns `FPCA_C` and `data`, so
# each subject's reconstruction error under EVERY class is recomputable exactly as
# `.rkcfc_objective` does it. Computing it inside the loop would make it free;
# post hoc costs n*K predictions and touches nothing that produced a banked result.
#
# WORK IS SLICED, not split contiguously: every item is seeded from its own
# identity, so the array can be resized and a failed task re-run alone.

source(file.path("R", "methods.R"))
source(file.path("R", "realdata36_settings.R"))

args    <- commandArgs(trailingOnly = TRUE)
args    <- args[args != "--args"]
task_id <- if (length(args) >= 1) as.integer(args[1]) else 1L
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 1L
outdir  <- if (length(args) >= 3) args[3] else file.path(REAL_RESULTS, "F4")
n_start <- if (length(args) >= 4) as.integer(args[4]) else 100L
n_sub   <- if (length(args) >= 5) as.integer(args[5]) else 500L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cfg <- load_config()
K   <- gint(cfg, "K_true")
# SETTINGS FROM realdata36_settings.R. The first F4 run used bw = 3.6 on [0,36]
# with bwK = "GCV" on the RAW scale and lost 560 of 604 items (92.7%) to "No
# enough points in local window". bwK, not bw, is the binding parameter; set it
# from stage F0 (RKCFC36_BWK) before running this.
BW  <- REAL_BW
BWK <- REAL_BWK
# THE START FOR THE SUBSAMPLES AND THE K-SELECTION FITS -- the one change from
# the 18-month F4 beyond the settings. That F4 used the configured two-staged
# start (k-means on the Kendall scores of a fit to all subjects). On the 36-month
# window it stops at once with "No enough points in local window" at every
# bandwidth tried (8, 10 and 12; local checks 2026-10-02, K = 2): its first
# partition leaves one class with no observations late in follow-up. Random
# starts, the start F1 uses, fitted 7 of 7 at bwK 8 and 10. So the subsamples and
# the K-selection fits start at random here, like the 100 restarts. Stage F0
# measures both starts; RKCFC36_F4_INIT=two-staged restores the 18-month choice.
F4_INIT <- Sys.getenv("RKCFC36_F4_INIT", "random")
if (!F4_INIT %in% c("random", "two-staged"))
  stop("RKCFC36_F4_INIT must be random or two-staged, got '", F4_INIT, "'", call. = FALSE)
# F4 fits on the scale of the PRIMARY partition, so its stability statement is
# about the partition the paper reports: "raw_M1" -> "raw", "log_M1" -> "log".
# The 2026-09-09 run was raw-scale. On that scale the objective's minimum turned
# out to be a 58-patient extreme group (J 402360, reached by 1 start in 94)
# rather than the modal 1199/196 split -- squared error on CEA with excess
# kurtosis 501 is dominated by the few patients in the hundreds -- so if
# REAL_PRIMARY moves to the log scale this block must be re-run there, and the
# combiner checks that the scales agree before it uses F1's partition.
SCALE <- if (exists("REAL_PRIMARY")) sub("_.*$", "", REAL_PRIMARY) else "raw"
d   <- truncate_cea36(load_cea36(scale = SCALE), REAL_TMAX)
iv  <- d$interval
n   <- length(d$Ly)

optns <- realdata_optns(cfg, interval = iv, bw = BW, bwK = BWK)

# --- the work list ------------------------------------------------------------
work <- rbind(
  data.frame(kind = "restart",   idx = seq_len(n_start), stringsAsFactors = FALSE),
  data.frame(kind = "subsample", idx = seq_len(n_sub),   stringsAsFactors = FALSE),
  data.frame(kind = "kselect",   idx = 2:5,              stringsAsFactors = FALSE))
mine <- work[(seq_len(nrow(work)) - 1L) %% n_tasks + 1L == task_id, , drop = FALSE]

cat(sprintf("F4-36 | task %d of %d | %d of %d items | n = %d | subsample/K-selection start: %s\n",
            task_id, n_tasks, nrow(mine), nrow(work), n, F4_INIT))
if (!nrow(mine))
  stop(sprintf(paste("task %d of %d was assigned NO items (there are %d).",
                     "Pass the ORIGINAL array size as n_tasks when re-running one task."),
               task_id, n_tasks, nrow(work)), call. = FALSE)

fit_one <- function(Ly, Lt, Kfit, init, seed) {
  set.seed(seed)
  tryCatch(RKCFC(Ly = Ly, Lt = Lt, K = Kfit, robust = TRUE, optns = optns,
                 init.method = init, maxIter = gint(cfg, "maxIter"),
                 conv.prop = gnum(cfg, "conv_prop"), LOO = glog(cfg, "LOO"),
                 parallel = FALSE, trace = 0),
           error = function(e) structure(conditionMessage(e), class = "f4fail"))
}

# Per-subject reconstruction error under EVERY class, from a fitted object.
class_errors <- function(f, Ly, Lt) {
  Kf <- length(f$FPCA_C)
  m <- matrix(NA_real_, nrow = length(Ly), ncol = Kf)
  for (k in seq_len(Kf)) for (i in seq_along(Ly)) {
    p <- tryCatch(univar_predicted_yhat(f$FPCA_C[[k]], Ly[[i]], Lt[[i]], TRUE, optns),
                  error = function(e) NULL)
    if (!is.null(p)) m[i, k] <- sum((Ly[[i]] - p$fittedX)^2)
  }
  m
}

res <- list(); kk <- 0L
t0 <- Sys.time()

for (i in seq_len(nrow(mine))) {
  kind <- mine$kind[i]; idx <- mine$idx[i]
  seed <- cell_seed(gint(cfg, "seed_master"), paste0("F4|", kind), idx)

  if (kind == "restart") {
    f <- fit_one(d$Ly, d$Lt, K, "random", seed)
    kk <- kk + 1L
    res[[kk]] <- list(kind = kind, idx = idx,
      error = if (inherits(f, "f4fail")) as.character(f) else NA_character_,
      g = if (inherits(f, "f4fail")) NULL else as.integer(f$res_clust$g),
      J = if (inherits(f, "f4fail")) NA_real_ else f$objective_final,
      stop_reason = if (inherits(f, "f4fail")) NA_character_ else f$stop_reason,
      # The margin is computed once, on the FIRST restart only: it is a property
      # of a fitted partition, not of the restart scheme, and recomputing it 100
      # times would cost n*K predictions each for one number.
      margin = if (!inherits(f, "f4fail") && idx == 1L) {
        ce <- class_errors(f, d$Ly, d$Lt)
        s <- t(apply(ce, 1, sort)); (s[, 2] - s[, 1]) / pmax(s[, 1], .Machine$double.eps)
      } else NULL)

  } else if (kind == "subsample") {
    set.seed(seed)
    keep <- sort(sample.int(n, floor(0.8 * n)))
    f <- fit_one(d$Ly[keep], d$Lt[keep], K, F4_INIT, seed)
    kk <- kk + 1L
    res[[kk]] <- list(kind = kind, idx = idx, keep = keep,
      error = if (inherits(f, "f4fail")) as.character(f) else NA_character_,
      g = if (inherits(f, "f4fail")) NULL else as.integer(f$res_clust$g))

  } else {                                   # kselect
    f <- fit_one(d$Ly, d$Lt, idx, F4_INIT, seed)
    r <- list(kind = kind, idx = idx,
              error = if (inherits(f, "f4fail")) as.character(f) else NA_character_)
    if (!inherits(f, "f4fail")) {
      g <- as.integer(f$res_clust$g)
      r$g <- g
      # THE CORRECTED SILHOUETTE: the class mean functions behind the
      # coefficient must be smoothed at the SAME bandwidth as the model being
      # scored. silhouette_score() passes bwK in the robust branch since 0.2.0.
      # silhouette_score() takes a LONG-FORMAT data frame and the option list is
      # named `option`. The 2026-09-10 runs called it as (Ly, Lt, g, optns = ...),
      # which raised "unused argument" inside a tryCatch that returned NA_real_
      # and dropped the message, so BOTH F4 runs printed "-" for every K without
      # saying why. The error is now kept in `sil_error` so a blank is explained.
      long <- data.frame(id = rep(seq_along(d$Ly), lengths(d$Ly)),
                         g = rep(g, lengths(d$Ly)),
                         time = unlist(d$Lt), y = unlist(d$Ly))
      sil_try <- tryCatch(silhouette_score(long, class = "g", id = "id", y = "y",
                                           time = "time", robust = TRUE, option = optns),
                          error = function(e) e)
      if (inherits(sil_try, "error")) {
        r$sil <- NA_real_; r$sil_error <- conditionMessage(sil_try)
      } else {
        r$sil <- as.numeric(sil_try$silhouette); r$sil_error <- NA_character_
      }
      r$J <- f$objective_final
    }
    kk <- kk + 1L; res[[kk]] <- r
  }

  cat(sprintf("  %d/%d  %s %d  (%.1f min elapsed)\n", i, nrow(mine), kind, idx,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

saveRDS(list(items = res, n = n, K = K, n_tasks = n_tasks,
             n_start = n_start, n_sub = n_sub, bw = BW, bwK = BWK, Tmax = REAL_TMAX,
             init_subsample_kselect = F4_INIT,
             scale = SCALE,
             code_hash = source_hash(c(file.path("R", "run_realdata36_f4.R"),
                                       REAL36_SETTINGS_FILES))),
        file.path(outdir, sprintf("F4_task%04d.rds", task_id)))
cat(sprintf("wrote %s (%d items, %.1f min)\n",
            file.path(outdir, sprintf("F4_task%04d.rds", task_id)), kk,
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
