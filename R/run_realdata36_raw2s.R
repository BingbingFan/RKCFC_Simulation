# 36-MONTH WINDOW, ORIGINAL SCALE, TWO-STAGED START: LONGER RUNS (one fit per task).
#
#   Rscript R/run_realdata36_raw2s.R <task> <stage: 200|500> [root] [n_sub]
#     root   default <REAL_RESULTS>/RAW2S; stage 200 writes root/maxit200, stage 500
#            writes root/maxit500
#     n_sub  fit a random subsample of this many patients -- SMOKE ONLY (no warm start)
#
# WHY. In the F0 sweep (2026-10-02) the two-staged start on the original scale ran only
# at Kendall-surface bandwidths of 10 and 12 months, and every original-scale fit that
# ran, from either start, stopped at the 100-iteration cap. The main run's 20 random
# starts on the original scale ended in 16 different partitions. These runs ask whether
# more iterations let such fits converge, whether they agree with one another, and
# whether a fit started near a good partition stays there.
#
# THE CELLS (task 1..NTASK), from four environment knobs. Each is '+'-separated and must
# never contain a comma, because some job schedulers cut a value at a comma:
#   RKCFC36_RAW2S_TRANSFORM  default "raw"      measurement scale(s): raw = the original
#                            scale; std = (x - median)/IQR, pooled over all measurements
#                            (the package's robust_normalize()); sqrt; log. An affine
#                            map such as std leaves every step of RKCFC unchanged (the
#                            Kendall surface signs pairwise differences, the mean and the
#                            scores are least squares, the assignment rule compares squared
#                            errors), so std must reproduce raw partition for partition;
#                            the reader checks that when both are present.
#   RKCFC36_RAW2S_BWK        default "10+12"    Kendall-surface bandwidths, months
#   RKCFC36_RAW2S_NSEED      default "10"       two-staged starts per cell; seeds follow
#                            the sweep -- cell_seed(seed_master, "F0-36|raw|two-staged|
#                            <bwK>", seed) -- so seeds 1-3 ARE the sweep's fits run on
#                            (and seeds 1 and 3 at 10 months fail again, with the
#                            local-window error, since the start is deterministic given
#                            the seed). The same seeds are used for every transform.
#   RKCFC36_RAW2S_WARM       default "log+raw"  warm starts per cell, run with
#                            init.labels: "log" = the REPORTED log-scale partition, i.e.
#                            the smallest J over the F1 and F4 log-scale starts (the
#                            selection rule of Section 2.3; expect 826/579, J 1077.95);
#                            "raw" = the main run's best original-scale start, the
#                            smallest J over the F1 raw starts (expect 1123/282, J
#                            658545.7, which stopped at the 100-iteration cap); "none"
#                            for no warm start. Both are read from <REAL_RESULTS>/F1 and
#                            /F4.
#   RKCFC36_RAW2S_NCORE      default 1  cores per fit, given to
#                            the package's parallel reassignment step (parLapply over
#                            patients). The Kendall FPCA step stays serial, so five cores
#                            make a fit about 10-15% faster; the partition does not
#                            depend on the core count.
#   NTASK = transforms x bandwidths x (NSEED + warm starts); the submit script's default
#   is 1 x 4 x 12 = 48 (bandwidths 10, 12, 18 and 24 months).
#   Order: transform, then bandwidth, then the NSEED seeds, then the warm starts.
#   hpc/submit_realdata36_raw2s.sh --list prints the layout.
#
# TWO STAGES, ONE AFTER THE OTHER.
#   200  a fresh fit, at most 200 iterations.
#   500  CONTINUES each stage-200 fit that stopped at the cap: RKCFC is restarted with
#        init.labels = its partition at iteration 200 and maxIter 301, which runs
#        iterations 201-500 of the same trajectory (the iteration is deterministic given
#        the partition: LOO is off and nothing in the loop draws random numbers). A fit
#        that converged, cycled, failed or lost a class in stage 200 is not refitted; its
#        stage-200 result is carried into the stage-500 folder unchanged, so that folder
#        holds every cell's result at a maximum of 500 iterations. One difference from
#        an uninterrupted 500-iteration run: the package's cycle check only sees the
#        partitions of the current call, so a cycle straddling iteration 200 is caught
#        one cycle length later, on the same partition; the reader detects it from the
#        saved histories. This saves the 200 repeated iterations of every continued fit.
#
# SETTINGS, as in the sweep: RKCFC (robust), window [0, 36], h' 3.6 months, K, conv_prop
# (0: an exact repeat of the previous partition) and LOO from the configuration. NCORE
# cores per fit (5 in the reported run): the package parallelizes only the
# reassignment step, 18% of an iteration on these data (timed locally 2026-10-02:
# Kendall FPCA 23.5 s, reassignment 5.2 s), so five cores make a fit about 10-15%
# faster. The partition does not depend on the core
# count: each patient's reassignment is a pure function of the current class fits, and
# nothing in the loop draws random numbers.
#
# WHAT IS SAVED (no CEA values, no identifiers): the starting and final partitions, J,
# iterations (counted from the start of stage 200), stop reason, cycle length, the
# package's per-iteration objective trace, the number of patients changing class at
# each iteration, EVERY partition visited (an integer matrix, patients x iterations,
# so that convergence under any tolerance, the oscillating patients and cycles across
# the stage boundary can be read off afterwards without refitting), and the in-sample
# and held-out relative reconstruction errors of the final partition, on the fitted
# scale, scored as F1 scores a fit.

source(file.path("R", "methods.R"))
source(file.path("R", "realdata36_settings.R"))

args    <- commandArgs(trailingOnly = TRUE); args <- args[args != "--args"]
task_id <- if (length(args) >= 1) as.integer(args[1]) else 1L
stage   <- if (length(args) >= 2) as.integer(args[2]) else 200L
# The warm starts come from the LOG-SCALE main run of 2026-10-02 and the sweep, which live
# in results_real36 whatever the settings file's default results folder is (it moved to
# results_real36_raw24 on 2026-10-03, when the original scale became the primary analysis).
RAW2S_SRC <- Sys.getenv("RKCFC36_RAW2S_SRC", if (dir.exists("results_real36")) "results_real36" else REAL_RESULTS)
root    <- if (length(args) >= 3) args[3] else file.path(RAW2S_SRC, "RAW2S")
n_sub   <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_

split_plus <- function(name, default) {
  v <- Sys.getenv(name, ""); if (!nzchar(v)) v <- default
  v <- strsplit(v, "+", fixed = TRUE)[[1]]; v <- trimws(v[nzchar(trimws(v))])
  if (!length(v)) stop(name, " is empty", call. = FALSE)
  v
}
TRANSFORMS <- split_plus("RKCFC36_RAW2S_TRANSFORM", "raw")
BWKS       <- as.numeric(split_plus("RKCFC36_RAW2S_BWK", "10+12"))
NSEED      <- as.integer(Sys.getenv("RKCFC36_RAW2S_NSEED", "10"))
WARMS      <- split_plus("RKCFC36_RAW2S_WARM", "log+raw"); if (identical(WARMS, "none")) WARMS <- character(0)
# TEST ONLY: RKCFC36_RAW2S_TEST_CAPS="3+5" shrinks the two caps for a quick run; files
# written under it are marked as smoke runs, which the reader refuses by default.
CAPS      <- as.integer(split_plus("RKCFC36_RAW2S_TEST_CAPS", "200+500"))
TEST_CAPS <- !identical(CAPS, c(200L, 500L))
NCORE     <- as.integer(Sys.getenv("RKCFC36_RAW2S_NCORE", "1"))
if (is.na(NCORE) || NCORE < 1L) stop("RKCFC36_RAW2S_NCORE must be a positive integer", call. = FALSE)
bad <- setdiff(TRANSFORMS, c("raw", "std", "sqrt", "log")); if (length(bad)) stop("unknown transform: ", paste(bad, collapse = " "), call. = FALSE)
bad <- setdiff(WARMS, c("log", "raw"));                      if (length(bad)) stop("unknown warm start: ", paste(bad, collapse = " "), call. = FALSE)
if (any(!is.finite(BWKS)) || any(BWKS <= 0)) stop("RKCFC36_RAW2S_BWK must be positive numbers", call. = FALSE)
if (is.na(NSEED) || NSEED < 0L) stop("RKCFC36_RAW2S_NSEED must be a non-negative integer", call. = FALSE)
if (length(CAPS) != 2L || any(is.na(CAPS)) || CAPS[1] < 2L || CAPS[2] <= CAPS[1]) stop("caps must be two integers, the second larger", call. = FALSE)
if (!stage %in% c(200L, 500L)) stop("stage must be 200 or 500", call. = FALSE)

# --- the cell table: one row per task ----------------------------------------------
cells <- do.call(rbind, lapply(TRANSFORMS, function(tr) do.call(rbind, lapply(BWKS, function(b) rbind(
  if (NSEED > 0L) data.frame(transform = tr, bwK = b, start = "two-staged", seed_idx = seq_len(NSEED), stringsAsFactors = FALSE),
  if (length(WARMS)) data.frame(transform = tr, bwK = b, start = paste0("warm-", WARMS), seed_idx = 0L, stringsAsFactors = FALSE))))))
if (is.null(cells) || !nrow(cells)) stop("no cells: set NSEED > 0 or a warm start", call. = FALSE)
NTASK <- nrow(cells)
if (is.na(task_id) || task_id < 1L || task_id > NTASK) stop(sprintf("task must be 1..%d", NTASK), call. = FALSE)
cell  <- cells[task_id, ]
TRANSFORM <- cell$transform; BWK <- cell$bwK; START <- cell$start; SEED <- cell$seed_idx
fname <- sprintf("RAW2S_%s_bwK%02d_%s.rds", TRANSFORM, as.integer(BWK),
                 if (START == "two-staged") sprintf("2s_seed%03d", SEED) else sub("-", "", START, fixed = TRUE))
outdir <- file.path(root, sprintf("maxit%d", stage))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cfg <- load_config()
K   <- gint(cfg, "K_true")
d   <- truncate_cea36(load_cea36(scale = "raw"), REAL_TMAX)
if (!is.na(n_sub)) {
  cat(sprintf("*** SUBSAMPLE OF %d SUBJECTS -- SMOKE RUN, NOT REPORTABLE ***\n", n_sub))
  if (START != "two-staged") stop("a warm start needs the full cohort; no n_sub with a warm start", call. = FALSE)
  set.seed(20260909L)
  keep <- sort(sample.int(length(d$Ly), min(n_sub, length(d$Ly))))
  d$Ly <- d$Ly[keep]; d$Lt <- d$Lt[keep]; d$ids <- d$ids[keep]; d$base <- d$base[keep, , drop = FALSE]
}
n <- length(d$Ly)

# --- the measurement scale ---------------------------------------------------------
transform_params <- NULL
d$Ly <- switch(TRANSFORM,
  raw  = d$Ly,
  std  = { z <- RKCFC::robust_normalize(unlist(d$Ly, use.names = FALSE))
           transform_params <- list(centre = attr(z, "centre"), scale = attr(z, "scale"))
           lapply(d$Ly, function(v) (v - transform_params$centre) / transform_params$scale) },
  sqrt = lapply(d$Ly, sqrt),
  log  = lapply(d$Ly, log))
d$scale <- TRANSFORM

optns <- realdata_optns(cfg, interval = d$interval, bw = REAL_BW, bwK = BWK)
# the same seed for every transform: the two-staged start then draws the same fallback
# labels, so transforms are comparable seed by seed
seed  <- cell_seed(gint(cfg, "seed_master"), sprintf("F0-36|%s|%s|%g", "raw", "two-staged", BWK),
                   if (START == "two-staged") SEED else 1L)
out <- list(name = sub("\\.rds$", "", fname), task = task_id, ntask = NTASK, stage = stage,
            transform = TRANSFORM, transform_params = transform_params, scale = TRANSFORM,
            init = START, bwK = BWK, seed_idx = SEED, seed = seed, nseed = NSEED,
            bwks = BWKS, warms = WARMS, transforms = TRANSFORMS,
            n_subjects = n, Tmax = REAL_TMAX, bw = REAL_BW, K = K, caps = CAPS, ncore = NCORE,
            smoke = !is.na(n_sub) || TEST_CAPS)
if (TEST_CAPS) cat(sprintf("*** TEST CAPS %d and %d -- SMOKE RUN, NOT REPORTABLE ***\n", CAPS[1], CAPS[2]))
cat(sprintf("RAW2S | stage %d | task %d of %d | %s | bwK %g | %s%s | n = %d | window [0, %g] | bw %g | %d core%s\n",
            stage, task_id, NTASK, TRANSFORM, BWK, START,
            if (START == "two-staged") sprintf(" seed %d (seed value %d)", SEED, seed) else "",
            n, REAL_TMAX, REAL_BW, NCORE, if (NCORE == 1L) "" else "s"))

save_out <- function(out) {
  out$code_hash <- source_hash(c(file.path("R", "run_realdata36_raw2s.R"), REAL36_SETTINGS_FILES))
  saveRDS(out, file.path(outdir, fname))
  cat(sprintf("wrote %s\n", file.path(outdir, fname)))
}

# --- the warm starts, read from the main run's files on this machine ---------------
no_error <- function(z) is.null(z$error) || isTRUE(is.na(z$error))
warm_partition <- function(kind) {
  f1dir <- file.path(RAW2S_SRC, "F1")
  if (kind == "log") {
    s1 <- lapply(list.files(f1dir, "^F1_log_M1_start[0-9]+\\.rds$", full.names = TRUE), readRDS)
    s1 <- Filter(function(z) isTRUE(z$ok) && length(z$g) == n, s1)
    ff <- list.files(file.path(RAW2S_SRC, "F4"), "^F4_task[0-9]+\\.rds$", full.names = TRUE)
    rs <- list()
    if (length(ff)) {
      z4 <- readRDS(ff[1])
      if (!identical(z4$scale, "log") || !isTRUE(z4$Tmax == REAL_TMAX) || !isTRUE(z4$K == K))
        stop("F4 under ", REAL_RESULTS, " is not the log-scale 36-month run", call. = FALSE)
      rs <- Filter(function(z) identical(z$kind, "restart") && no_error(z) && length(z$g) == n,
                   unlist(lapply(ff, function(f) readRDS(f)$items), recursive = FALSE))
    }
    pool <- c(lapply(s1, function(z) list(source = "F1", start = z$start, g = as.integer(z$g), J = as.numeric(z$J))),
              lapply(rs, function(z) list(source = "F4", start = z$idx, g = as.integer(z$g), J = as.numeric(z$J))))
    rule <- "smallest J over the log-scale RKCFC starts of F1 and F4 (the selection rule of Section 2.3)"
  } else {
    s1 <- lapply(list.files(f1dir, "^F1_raw_M1_start[0-9]+\\.rds$", full.names = TRUE), readRDS)
    s1 <- Filter(function(z) isTRUE(z$ok) && length(z$g) == n, s1)
    pool <- lapply(s1, function(z) list(source = "F1", start = z$start, g = as.integer(z$g), J = as.numeric(z$J)))
    rule <- "smallest J over the original-scale RKCFC starts of F1"
  }
  if (!length(pool)) stop("no ", kind, "-scale starts found under ", REAL_RESULTS, " for the warm start", call. = FALSE)
  J <- vapply(pool, `[[`, 0, "J"); b <- which.min(J)
  list(g = pool[[b]]$g, info = list(kind = kind, rule = rule, n_starts = length(pool), source = pool[[b]]$source,
                                    start = pool[[b]]$start, J = J[b], sizes = as.integer(table(factor(pool[[b]]$g, levels = seq_len(K))))))
}

# --- what to start from --------------------------------------------------------------
init_labels <- NULL; iter_offset <- 0L; max_it <- if (stage == 200L) CAPS[1] else CAPS[2]
if (stage == 200L && START != "two-staged") {
  w <- warm_partition(sub("^warm-", "", START))
  init_labels <- w$g; out$warm <- w$info
  cat(sprintf("  warm start: %s -> %s start %s, J %.4f, sizes %s\n", w$info$rule, w$info$source,
              as.character(w$info$start), w$info$J, paste(w$info$sizes, collapse = "/")))
}
if (stage == 500L) {
  f200 <- file.path(root, "maxit200", fname)
  if (!file.exists(f200)) stop("stage-200 result not found: ", f200, " -- stage 200 must finish first", call. = FALSE)
  r1 <- readRDS(f200)
  if (!isTRUE(r1$ok) || !identical(r1$stop_reason, "maxIter")) {
    r1$stage <- stage; r1$continued <- FALSE
    r1$carried_reason <- if (!isTRUE(r1$ok)) "failed in stage 200" else sprintf("stopped in stage 200 (%s)", r1$stop_reason)
    cat(sprintf("  not continued: %s; stage-200 result carried forward\n", r1$carried_reason))
    save_out(r1)
    quit(save = "no", status = 0)
  }
  init_labels <- r1$g; iter_offset <- r1$iter - 1L; max_it <- CAPS[2] - r1$iter + 1L
  out$g_stage200 <- r1$g; out$warm <- r1$warm
  cat(sprintf("  continuing from the stage-200 partition (%s) for up to %d more iterations\n",
              paste(table(r1$g), collapse = "/"), max_it - 1L))
}

# THE WORKER SOCKET. With NCORE > 1 the package calls parallel::makeCluster(), which
# opens a server socket on the port the parallel package chose when it was loaded (a
# draw from 11000-11999). Tasks of one array land several to a node, so two of them can
# hold the same port, and the second fails at once with "creation of server socket
# failed: port ... cannot be opened". That killed the two warm-raw continuations of the
# first real run (2026-10-03, stage 500, tasks 12 and 24) before any iteration ran. So
# the port is set per task and per attempt here, and that one failure, which costs
# nothing, is retried on another port after a short wait. The RNG is re-seeded before
# each attempt, so a retried two-staged start draws the same labels.
is_socket_failure <- function(f) inherits(f, "rawfail") &&
  grepl("socket|cannot be opened|cannot open the connection", as.character(f), ignore.case = TRUE)
fit_once <- function(attempt) {
  # not exported from parallel, hence :::; it sets the port makeCluster() reads when none is passed
  if (NCORE > 1L) parallel:::setDefaultClusterOptions(port = 11000L + ((task_id * 7L + attempt * 131L) %% 1000L))
  set.seed(seed)
  tryCatch(
    RKCFC(Ly = d$Ly, Lt = d$Lt, K = K, robust = TRUE, optns = optns,
          init.method = "two-staged", init.labels = init_labels, maxIter = max_it,
          conv.prop = gnum(cfg, "conv_prop"), LOO = glog(cfg, "LOO"),
          parallel = NCORE > 1L, nmcores = if (NCORE > 1L) NCORE else NULL, trace = 0),
    error = function(e) structure(conditionMessage(e), class = "rawfail"))
}
t0 <- Sys.time()
for (attempt in 1:6) {
  f <- fit_once(attempt)
  if (!is_socket_failure(f) || attempt == 6L) break
  cat(sprintf("  worker processes could not start (attempt %d of 6): %s -- retrying on another port\n",
              attempt, as.character(f)))
  Sys.sleep(10 + 5 * attempt)
}
out$attempts <- attempt
secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
out$seconds <- secs; out$continued <- stage == 500L; out$iter_offset <- iter_offset

if (inherits(f, "rawfail")) {
  out$ok <- FALSE; out$error <- as.character(f)
  cat(sprintf("  FAILED after %.1f min: %s\n", secs / 60, out$error))
} else {
  g  <- as.integer(f$res_clust$g); cl <- f$clustlist
  out$ok <- TRUE; out$error <- NA_character_
  out$g_start <- as.integer(cl[[1]])
  out$g <- g; out$iter <- iter_offset + f$iter; out$conv <- f$conv; out$stop_reason <- f$stop_reason
  out$cycle_length <- f$cycle_length; out$J <- f$objective_final
  out$sizes <- as.integer(table(factor(g, levels = seq_len(K))))
  ot <- f$objective; if (!is.null(ot) && "iteration" %in% names(ot)) ot$iteration <- ot$iteration + iter_offset
  out$objective_trace <- ot
  # every partition visited in THIS call: column j is iteration iter_offset + j
  out$partitions <- matrix(as.integer(unlist(cl, use.names = FALSE)), nrow = n)
  out$n_changed <- if (ncol(out$partitions) >= 2L)
    vapply(2:ncol(out$partitions), function(j) as.integer(round(
      mclust::classError(out$partitions[, j - 1L], out$partitions[, j])$errorRate * n)), 0L) else integer(0)
  out$rre_in <- tryCatch(reconstruction_errors(d, g, f$FPCA_C, TRUE, optns)$rre_in, error = function(e) NA_real_)
  sp <- holdout_split(d, gnum(cfg, "holdout_frac"), gint(cfg, "holdout_min_ni"), 1L)
  out$rre_out <- tryCatch({
    fh <- univar_FPCA_k(sp$fit$Ly, sp$fit$Lt, K, g, TRUE, optns)
    holdout_error(sp, g, fh, TRUE, optns)
  }, error = function(e) NA_real_)
  cat(sprintf("  ok after %.1f min: %s at iteration %d%s, sizes %s, J %.4f, rre_in %.4f, rre_out %s; patients changing class in the last 10 iterations: %s\n",
              secs / 60, out$stop_reason, out$iter,
              if (is.na(out$cycle_length)) "" else sprintf(" (cycle of %d)", out$cycle_length),
              paste(out$sizes, collapse = "/"), out$J, out$rre_in,
              if (is.finite(out$rre_out)) sprintf("%.4f", out$rre_out) else "NA",
              paste(utils::tail(out$n_changed, 10), collapse = " ")))
}
save_out(out)
