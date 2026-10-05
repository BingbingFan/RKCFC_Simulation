# EVERY DESIGN PARAMETER THAT WAS DECIDED ON THIN EVIDENCE, RE-RUN AT SIZE.
#
#   Rscript R/run_decisions.R <factor> <task_id> [n_tasks] [outdir] [n_reps] [n]
#
#   factor   nkfit | d | delta | convprop | init
#
# WHY THIS BLOCK EXISTS. Five parameters were settled on 2026-09-05 from ONE
# replicate, or three, at n = 60 or n = 500. Each is defensible as a direction
# and none is reportable as a number. A reader who asks "on how many
# replicates?" must not get "one". This block re-runs all five at n = 500 with
# 20 replicates so every choice in the manuscript rests on a measurement.
#
# WHAT EACH FACTOR SETTLES, AND WHAT IT CURRENTLY RESTS ON
#
#   nkfit     nK_fit in {2, 3, 4}. Reverted 4 -> 2 on ONE replicate per cell.
#             M1 was cap-bound at 4 (iter 100/100) and fell 0.968 -> 0.719 on
#             arm 1 and 0.945 -> 0.608 on arm 4. Also re-prices the FVE rule,
#             which selects 4 and DISAGREES with the accuracy evidence -- the
#             disagreement is itself a reportable finding.
#   d         d in {1.5, 2, 3} at the chosen nK_fit. d_easy moved 3 -> 2 on the
#             DELTA sweep alone; the accuracy confirmation was never run. d = 3
#             is believed saturated (M1 and M2 both returned 0.968 on arm 1,
#             identical to three decimals) and that needs 20 replicates, not 1.
#   delta     Delta across all six arms over a d grid. Currently 3 replicates.
#             This is the IQR-matching evidence and the reference-CCR column
#             the manuscript will quote, so it needs the same replication as
#             anything else reported.
#   convprop  conv.prop in {0, 0.001, 0.005, 0.01} with agree_strict. NEVER RUN.
#             The manuscript specifies 0; the question is whether the resulting
#             non-convergence is COSMETIC -- a few boundary subjects oscillating
#             forever -- or whether the partition genuinely has not settled.
#   init      the five initialisation schemes. Run once at n = 60, where the
#             intended robust start scored at CHANCE and the 0.1.x start
#             scored 1.000. That is the opposite of the expected direction and
#             it is one replicate at a tenth of the real n.
#
# TWO RULES THIS DRIVER FOLLOWS, BOTH LEARNED THE HARD WAY
#
# 1. EVERY LEVEL OF A FACTOR RUNS INSIDE ONE TASK. Slicing is on (arm, rep), so
#    all levels see the identical dataset and identical fitting seed and a
#    difference between levels is the level, never the draw. It also keeps
#    agree_strict from crossing a task boundary, which would make it
#    uncomputable.
# 2. THE FAILURE ROW HAS THE SUCCESS ROW'S COLUMNS. do.call(rbind, ...) over
#    data frames with different column sets is a hard error, so a short failure
#    row takes down a whole task after its compute is spent.

source(file.path("R", "methods.R"))
source(file.path("R", "separation.R"))   # brings generate.R and metrics.R with it

args    <- commandArgs(trailingOnly = TRUE)
args    <- args[args != "--args"]
factor_ <- if (length(args) >= 1) args[1] else "nkfit"
task_id <- if (length(args) >= 2) as.integer(args[2]) else 1L
n_tasks <- if (length(args) >= 3) as.integer(args[3]) else 1L
outdir  <- if (length(args) >= 4) args[4] else file.path("results", paste0("D_", factor_))
n_reps  <- if (length(args) >= 5) as.integer(args[5]) else 20L
n_over  <- if (length(args) >= 6) as.integer(args[6]) else NA_integer_

FACTORS <- c("nkfit", "d", "delta", "convprop", "init")
if (!factor_ %in% FACTORS)
  stop("unknown factor: ", factor_, ". Use one of ", paste(FACTORS, collapse = ", "),
       call. = FALSE)

cfg <- load_config()
K   <- gint(cfg, "K_true")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# `delta` needs every arm because it IS the cross-arm comparison. The fit-based
# factors use arms 1 and 4 -- the Gaussian reference and the infinite-variance
# arm the method exists for, which bracket the range and cost a third of the
# full set.
ARMS <- if (factor_ == "delta") 1:6 else c(1L, 4L)
D_GRID <- c(0, 1, 1.5, 2, 2.5, 3, 4)

jobs <- expand.grid(arm = ARMS, rep = seq_len(n_reps), KEEP.OUT.ATTRS = FALSE)
jobs$scenario_id <- sprintf("D%s_a%d", toupper(substr(factor_, 1, 3)), jobs$arm)
jobs <- jobs[order(jobs$rep, jobs$arm), ]
mine <- jobs[(seq_len(nrow(jobs)) - 1L) %% n_tasks + 1L == task_id, , drop = FALSE]

cat(sprintf("decisions/%s | task %d of %d | %d of %d cells | %d reps\n",
            factor_, task_id, n_tasks, nrow(mine), nrow(jobs), n_reps))
if (!is.na(n_over))
  cat(sprintf("*** n OVERRIDDEN to %d. SMOKE RUN -- not reportable. ***\n", n_over))
if (!nrow(mine))
  stop(sprintf(paste("task %d of %d was assigned NO cells (this factor has %d).",
                     "Pass the ORIGINAL array size as n_tasks when re-running one",
                     "task; with n_tasks = 1 every cell is reassigned."),
               task_id, n_tasks, nrow(jobs)), call. = FALSE)

# --- the row shape, defined once ----------------------------------------------

empty_row <- function(level = NA_character_, method = NA_character_,
                      error = NA_character_) {
  data.frame(
    factor = factor_, level = level, method = method, error = error,
    seconds = NA_real_, converged = NA, iter = NA_integer_,
    stop_reason = NA_character_, objective = NA_real_,
    ari = NA_real_, ccr = NA_real_, nmi = NA_real_,
    rre_in = NA_real_, rre_true = NA_real_,
    delta = NA_real_, ref_ccr = NA_real_, band = NA_character_,
    agree_strict = NA_real_,     # convprop only: agreement with the cp = 0 fit
    prop_best = NA_real_,        # init/random20 only
    init_method = NA_character_, stringsAsFactors = FALSE)
}

# One RKCFC fit, with every override this block needs, scored the same way.
one_fit <- function(dat, robust, nk = NULL, cp = NULL, init_m = NULL,
                    init_l = NULL, sid, rp) {
  cfg2 <- cfg
  if (!is.null(nk)) cfg2$global$nK_fit <- as.character(nk)
  optns <- build_optns(cfg2, robust)
  set.seed(cell_seed(gint(cfg, "seed_master"), paste0(sid, "|fit"), rp))
  t0 <- Sys.time()
  f <- tryCatch(
    RKCFC(Ly = dat$Ly, Lt = dat$Lt, K = K, robust = robust, optns = optns,
          init.method = if (is.null(init_m)) cfg$global$init_method else init_m,
          init.labels = init_l,
          maxIter = gint(cfg, "maxIter"),
          conv.prop = if (is.null(cp)) gnum(cfg, "conv_prop") else cp,
          LOO = glog(cfg, "LOO"), parallel = FALSE, trace = 0),
    error = function(e) structure(conditionMessage(e), class = "dfail"))
  list(fit = f, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")),
       optns = optns, robust = robust)
}

score_fit <- function(r, res, dat) {
  r$seconds <- res$seconds
  if (inherits(res$fit, "dfail")) { r$error <- as.character(res$fit); return(r) }
  f <- res$fit; g <- f$res_clust$g
  r$converged <- f$conv; r$iter <- f$iter; r$stop_reason <- f$stop_reason
  r$objective <- f$objective_final; r$init_method <- f$init.method
  r$ari <- mclust::adjustedRandIndex(dat$labels, g)
  r$ccr <- ccr(dat$labels, g); r$nmi <- nmi(dat$labels, g)
  e <- tryCatch(reconstruction_errors(dat, g, f$FPCA_C, res$robust, res$optns),
                error = function(err) list(rre_in = NA_real_, rre_true = NA_real_))
  r$rre_in <- e$rre_in; r$rre_true <- e$rre_true
  attr(r, "g") <- g
  r
}

# --- run ----------------------------------------------------------------------

rows <- list(); kk <- 0L
push <- function(jb, r) { kk <<- kk + 1L
  rows[[kk]] <<- cbind(jb, as.data.frame(r)[, names(empty_row())], row.names = NULL) }
t_start <- Sys.time()

for (i in seq_len(nrow(mine))) {
  jb  <- mine[i, ]
  arm <- cfg$arms[cfg$arms$arm == jb$arm, ]

  # `d` and `delta` vary d themselves; every other factor holds it at d_hard,
  # which is where the methods are closest and a setting can still matter.
  # For the `d` factor: the shape-only FLOOR, the two design levels, and 3.
  #   0   the classes still differ in EIGENFUNCTIONS -- sines against cosines --
  #       so they stay separable with no mean difference at all. This quantifies
  #       how much of the separation d never controlled, and a reader who assumes
  #       d is the whole story will misread every table without it.
  #   3   the RETIRED value, kept deliberately. It is believed saturated (M1 and
  #       M2 both returned 0.968 on arm 1, identical to three decimals) and that
  #       belief currently rests on ONE replicate. Retiring a level on one draw
  #       and never checking it again is how a design acquires folklore.
  d_levels <- switch(factor_,
                     d     = c(0, gnum(cfg, "d_hard"), gnum(cfg, "d_easy"), 3),
                     delta = D_GRID,
                     gnum(cfg, "d_hard"))

  for (dv in d_levels) {
    # THE SCENARIO ID DOES NOT INCLUDE d, AND THAT IS THE POINT. generate_replicate()
    # seeds from (scenario_id, rep), so holding the id fixed across the d levels
    # makes every draw -- labels, scores, observation counts, errors -- IDENTICAL,
    # and the ONLY thing that changes is the class mean. The d levels are then
    # PAIRED: the seed is held fixed across d so only the mean shifts.
    #
    # An earlier version appended _d<value> to the id, which reseeded at every
    # level and made the comparison across d unpaired -- adding draw-to-draw
    # variance to exactly the contrast the factor exists to measure.
    dat <- tryCatch(
      generate_replicate(cfg, arm, d = dv, rep = jb$rep,
                         scenario_id = jb$scenario_id,
                         n = if (is.na(n_over)) NULL else n_over),
      error = function(e) structure(conditionMessage(e), class = "genfail"))
    if (inherits(dat, "genfail")) {
      push(jb, empty_row(sprintf("d=%g", dv), NA_character_,
                         paste("GENERATION:", as.character(dat))))
      next
    }

    if (factor_ == "delta") {
      s <- tryCatch(separation_of(dat, cfg, robust = TRUE),
                    error = function(e) NULL)
      r <- empty_row(sprintf("d=%g", dv), "-")
      if (is.null(s)) r$error <- "separation failed"
      else { r$delta <- s$delta; r$ref_ccr <- s$ref_ccr; r$band <- s$band }
      push(jb, r); next
    }

    if (factor_ == "nkfit") {
      for (nk in 2:4) for (id in c("M1", "M2")) {
        res <- one_fit(dat, robust = id == "M1", nk = nk,
                       sid = jb$scenario_id, rp = jb$rep)
        push(jb, score_fit(empty_row(sprintf("nK=%d", nk), id), res, dat))
      }
    } else if (factor_ == "d") {
      for (id in c("M1", "M2")) {
        res <- one_fit(dat, robust = id == "M1", sid = jb$scenario_id, rp = jb$rep)
        push(jb, score_fit(empty_row(sprintf("d=%g", dv), id), res, dat))
      }
    } else if (factor_ == "convprop") {
      # agree_strict IS THE STATISTIC THIS FACTOR EXISTS FOR. Reporting aRI at
      # each threshold is not enough: two different partitions can score the
      # same aRI, so matching columns would not show the threshold left the
      # ANSWER alone. Every fit is therefore also scored against the cp = 0
      # partition FROM THE SAME REPLICATE.
      #   agree ~ 1   the partition settled long ago and a few boundary
      #               subjects oscillate forever -- non-convergence is COSMETIC
      #   agree << 1  the partition genuinely has not settled and no threshold
      #               disguises that
      for (id in c("M1", "M2")) {
        base_g <- NULL
        for (cp in c(0, 0.001, 0.005, 0.01)) {
          res <- one_fit(dat, robust = id == "M1", cp = cp,
                         sid = jb$scenario_id, rp = jb$rep)
          r <- score_fit(empty_row(sprintf("cp=%g", cp), id), res, dat)
          g <- attr(r, "g")
          if (cp == 0) base_g <- g
          if (!is.null(g) && !is.null(base_g))
            r$agree_strict <- 1 - mclust::classError(base_g, g)$errorRate
          push(jb, r)
        }
      }
    } else if (factor_ == "init") {
      sc_r <- tryCatch(global_scores(dat, cfg, TRUE, build_optns(cfg, TRUE)),
                       error = function(e) NULL)
      sc_n <- tryCatch(global_scores(dat, cfg, FALSE, build_optns(cfg, FALSE)),
                       error = function(e) NULL)
      # two_staged is the CONFIGURED default and belongs in the comparison. It
      # needs no branch: falling through the chain below leaves init_m NULL, and
      # one_fit() then uses cfg$global$init_method, which IS two-staged. The
      # 2026-09-06 run omitted it, so the default's own numbers had to be
      # borrowed from the nkfit/d/convprop baselines -- differently seeded, hence
      # unpaired against the five alternatives measured here.
      for (lv in c("two_staged", "kendall_kmeans", "pace_kmeans", "random20",
                   "ward", "oracle")) {
        r <- empty_row(lv, "M1")
        if (lv == "random20") {
          set.seed(cell_seed(gint(cfg, "seed_master"),
                             paste0(jb$scenario_id, "|fit"), jb$rep))
          t0 <- Sys.time()
          rs <- tryCatch(RKCFC_restarts(
                  Ly = dat$Ly, Lt = dat$Lt, K = K, robust = TRUE,
                  optns = build_optns(cfg, TRUE), nstart = 20L,
                  init.method = "random", maxIter = gint(cfg, "maxIter"),
                  conv.prop = gnum(cfg, "conv_prop"), LOO = glog(cfg, "LOO"),
                  parallel = FALSE, trace = 0),
                error = function(e) structure(conditionMessage(e), class = "dfail"))
          secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
          if (inherits(rs, "dfail")) { r$error <- as.character(rs); r$seconds <- secs }
          else {
            r <- score_fit(r, list(fit = rs$fit, seconds = secs,
                                   optns = build_optns(cfg, TRUE), robust = TRUE), dat)
            r$prop_best <- rs$prop_best
          }
          push(jb, r); next
        }
        init_m <- NULL; init_l <- NULL
        if (lv == "kendall_kmeans") init_m <- "kmeans"
        else if (lv == "pace_kmeans") {
          # THE 0.1.x START: start from NON-ROBUST PACE scores, then fit
          # robustly. That is exactly what the released 0.1.x package did
          # regardless of the `robust` flag, and it is not reachable through
          # init.method because RKCFC() passes its own `robust` to KCFC_init().
          if (is.null(sc_n)) { r$error <- "PACE scores failed"; push(jb, r); next }
          set.seed(cell_seed(gint(cfg, "seed_master"),
                             paste0(jb$scenario_id, "|init"), jb$rep))
          init_l <- stats::kmeans(sc_n, centers = K, nstart = 10)$cluster
        } else if (lv == "ward") {
          if (is.null(sc_r)) { r$error <- "scores failed"; push(jb, r); next }
          init_l <- as.integer(stats::cutree(
            stats::hclust(stats::dist(sc_r), method = "ward.D2"), k = K))
        } else if (lv == "oracle") init_l <- dat$labels
        res <- one_fit(dat, robust = TRUE, init_m = init_m, init_l = init_l,
                       sid = jb$scenario_id, rp = jb$rep)
        push(jb, score_fit(r, res, dat))
      }
    }
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
out$code_hash <- source_hash(file.path("R", "run_decisions.R"))
f <- file.path(outdir, sprintf("D_%s_task%04d.rds", factor_, task_id))
saveRDS(out, f)
cat(sprintf("wrote %s  (%d rows, %.1f min)\n", f, nrow(out),
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
