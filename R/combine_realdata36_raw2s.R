# READ THE LONGER ORIGINAL-SCALE RUNS (R/run_realdata36_raw2s.R).
#
#   Rscript R/combine_realdata36_raw2s.R [root] [allow_smoke]
#     root  default <REAL_RESULTS>/RAW2S, holding maxit200/ and (later) maxit500/
#
# Per stage, transform, bandwidth and kind of start: how many fits ran, converged,
# cycled, hit the cap or lost a class; iterations used; distinct final partitions and
# how alike they are; the best J. maxit500 holds every cell's result at a maximum of
# 500 iterations (continued fits refitted from their stage-200 partition, the rest
# carried from stage 200).
#
# FROM THE SAVED PARTITION HISTORIES, without refitting: (a) at which iteration each
# fit would have stopped under a 1% or 0.5% reclassification tolerance (the study's
# rule is 0, an exact repeat); (b) how many patients still change class per iteration
# at the end, and how many distinct patients moved at all in the last 50 iterations
# (the oscillating set); (c) whether a continued fit revisited one of its stage-200
# partitions, which is a cycle straddling the stage boundary that the package could
# not see.
#
# CHECKS: (1) two-staged seeds 1-3 on the original scale at iteration 100 against the F0
# sweep's fits (must be identical); (2) every continued fit started from its stage-200
# partition; (3) std (an affine rescaling) reproduces raw partition for partition when
# both transforms are present; (4) the best fits and the warm-start fits against the
# main run's best original-scale start (F1) and against the reported log-scale
# partition, and how far each warm start moved from where it began. Writes RAW2S.rds
# into root.

source(file.path("R", "realdata36_settings.R"))
if (!requireNamespace("mclust", quietly = TRUE)) stop("package `mclust` is required.", call. = FALSE)
args <- commandArgs(trailingOnly = TRUE); args <- args[args != "--args"]
# The RAW2S run, the F0 sweep and the log-scale F1 live in results_real36 whatever the
# settings file's default results folder is (it moved to results_real36_raw24 on 2026-10-03).
RAW2S_SRC <- if (dir.exists("results_real36")) "results_real36" else REAL_RESULTS
ROOT <- if (length(args) >= 1) args[1] else file.path(RAW2S_SRC, "RAW2S")
ALLOW_SMOKE <- length(args) >= 2 && identical(args[2], "allow_smoke")
ari <- function(a, b) mclust::adjustedRandIndex(a, b)
key <- function(g) paste(match(g, unique(g)), collapse = "")     # partition up to relabelling
fmt <- function(x, d = 2) if (length(x) == 0 || all(is.na(x))) "NA" else formatC(x, digits = d, format = "f")

rd <- function(m) {
  fs <- list.files(file.path(ROOT, sprintf("maxit%d", m)), "^RAW2S_.*\\.rds$", full.names = TRUE)
  r <- lapply(fs, readRDS); names(r) <- vapply(r, function(z) z$name, ""); r
}
runs <- list(`200` = rd(200), `500` = rd(500))
if (!length(runs[["200"]])) stop("no stage-200 files under ", ROOT, call. = FALSE)
if (!ALLOW_SMOKE && any(unlist(lapply(runs, function(r) vapply(r, function(z) isTRUE(z$smoke), TRUE)))))
  stop("smoke-run files found; refusing (pass allow_smoke to read them anyway)", call. = FALSE)
n   <- max(vapply(runs[["200"]], function(r) as.integer(r$n_subjects), 0L))
r0  <- runs[["200"]][[1]]
cat(sprintf("RAW2S under %s: %d stage-200 and %d stage-500 files (expected %d each); transforms %s; bwK %s; %d two-staged seeds; warm starts %s; n = %d; cores per fit %s\n",
            ROOT, length(runs[["200"]]), length(runs[["500"]]), r0$ntask, paste(r0$transforms, collapse = "+"),
            paste(r0$bwks, collapse = "+"), r0$nseed, if (length(r0$warms)) paste(r0$warms, collapse = "+") else "none", n,
            if (is.null(r0$ncore)) "1" else as.character(r0$ncore)))

# --- the full partition history of a fit: stage-200 columns, then the continuation's ---
history <- function(r) {
  h <- r$partitions
  if (is.null(h)) return(NULL)
  if (isTRUE(r$continued) && !is.null(r$g_stage200)) {
    p <- runs[["200"]][[r$name]]
    if (!is.null(p) && !is.null(p$partitions)) h <- cbind(p$partitions, h[, -1, drop = FALSE])
  }
  h
}
# consecutive reclassification counts (labels are consistent within a trajectory)
changes <- function(h) if (ncol(h) < 2) integer(0) else vapply(2:ncol(h), function(j) sum(h[, j] != h[, j - 1L]), 0L)
first_under <- function(ch, tol) { k <- which(ch <= tol * n); if (length(k)) k[1] + 1L else NA_integer_ }   # iteration index
posthoc <- function(r) {
  h <- history(r)
  if (is.null(h) || !isTRUE(r$ok)) return(list(tol1 = NA_integer_, tol05 = NA_integer_, chg_last = NA_real_, osc50 = NA_integer_, straddle = NA_integer_))
  ch <- changes(h)
  last <- h[, max(1L, ncol(h) - 50L):ncol(h), drop = FALSE]
  osc  <- sum(apply(last, 1L, function(v) any(v != v[1L])))
  strad <- NA_integer_
  if (isTRUE(r$continued) && !is.null(r$g_stage200) && !is.null(runs[["200"]][[r$name]]$partitions)) {
    p <- runs[["200"]][[r$name]]$partitions
    k200 <- apply(p, 2L, key); k500 <- apply(r$partitions[, -1, drop = FALSE], 2L, key)
    hit <- match(k500, k200)
    if (any(!is.na(hit))) { j <- which(!is.na(hit))[1]; strad <- as.integer(ncol(p) + j - hit[j]) }
  }
  list(tol1 = first_under(ch, 0.01), tol05 = first_under(ch, 0.005),
       chg_last = stats::median(utils::tail(ch, 20)), osc50 = osc, straddle = strad)
}

row_of <- function(m, r) {
  ph <- posthoc(r)
  data.frame(maxIter = as.integer(m), transform = r$transform, bwK = r$bwK, start = r$init, seed = r$seed_idx, ok = isTRUE(r$ok),
    continued = isTRUE(r$continued), stop = if (isTRUE(r$ok)) r$stop_reason else "FAILED",
    iter = if (isTRUE(r$ok)) r$iter else NA_integer_, cycle = if (isTRUE(r$ok)) r$cycle_length else NA_integer_,
    sizes = if (isTRUE(r$ok)) paste(r$sizes, collapse = "/") else "", J = if (isTRUE(r$ok)) r$J else NA_real_,
    rre_in = if (isTRUE(r$ok)) r$rre_in else NA_real_, rre_out = if (isTRUE(r$ok)) r$rre_out else NA_real_,
    hours = if (is.null(r$seconds)) NA_real_ else r$seconds / 3600,
    chg_last20 = ph$chg_last, osc_last50 = ph$osc50, iter_tol1 = ph$tol1, iter_tol05 = ph$tol05, straddle_cycle = ph$straddle,
    stringsAsFactors = FALSE)
}
tab <- do.call(rbind, unlist(lapply(names(runs), function(m) lapply(runs[[m]], function(r) row_of(m, r))), recursive = FALSE))
tab <- tab[order(tab$maxIter, tab$transform, tab$bwK, tab$start != "two-staged", tab$start, tab$seed), ]
rownames(tab) <- NULL

cat("\n== by stage, transform, bandwidth and start ==\n")
cat("  (iter_tol1 / iter_tol05: the iteration at which a 1% or 0.5% reclassification tolerance would have stopped the fit;\n")
cat("   chg_last20: median patients changing class per iteration over the last 20; osc_last50: distinct patients that moved in the last 50)\n")
groups <- unique(tab[, c("maxIter", "transform", "bwK", "start")])
for (i in seq_len(nrow(groups))) {
  gq <- groups[i, ]
  x <- tab[tab$maxIter == gq$maxIter & tab$transform == gq$transform & tab$bwK == gq$bwK & tab$start == gq$start, ]
  rr <- Filter(function(r) isTRUE(r$ok) && r$transform == gq$transform && r$bwK == gq$bwK && r$init == gq$start, runs[[as.character(gq$maxIter)]])
  cat(sprintf("  max %3d | %-4s | bwK %2g | %-10s: %d fits | ran %d, failed %d | converged %d, cycle %d, at cap %d, small class %d%s\n",
              gq$maxIter, gq$transform, gq$bwK, gq$start, nrow(x), sum(x$ok), sum(!x$ok), sum(x$stop == "converged"),
              sum(x$stop == "cycle"), sum(x$stop == "maxIter"), sum(x$stop == "small_cluster"),
              if (gq$maxIter == 500) sprintf(" | continued: %d", sum(x$continued)) else ""))
  if (length(rr)) {
    fin <- vapply(rr, function(r) key(r$g), "")
    stt <- if (gq$start == "two-staged" && gq$maxIter == 200) vapply(rr, function(r) key(r$g_start), "") else NULL
    pa  <- if (length(rr) > 1) utils::combn(length(rr), 2, function(ij) ari(rr[[ij[1]]]$g, rr[[ij[2]]]$g)) else NA
    best <- rr[[which.min(vapply(rr, `[[`, 0, "J"))]]
    cat(sprintf("        iterations %s | distinct starts %s, distinct final partitions %d | pairwise aRI median %s (%s to %s) | best J %.4f (%s)\n",
                paste(range(x$iter, na.rm = TRUE), collapse = "-"), if (length(stt)) as.character(length(unique(stt))) else "-",
                length(unique(fin)), fmt(stats::median(pa)), fmt(min(pa)), fmt(max(pa)), best$J, paste(best$sizes, collapse = "/")))
    cat(sprintf("        at the cap: %d | would have stopped under 1%%: %d (median iteration %s), under 0.5%%: %d (%s) | chg_last20 median %s | osc_last50 median %s\n",
                sum(x$stop == "maxIter"), sum(!is.na(x$iter_tol1)), fmt(stats::median(x$iter_tol1, na.rm = TRUE), 0),
                sum(!is.na(x$iter_tol05)), fmt(stats::median(x$iter_tol05, na.rm = TRUE), 0),
                fmt(stats::median(x$chg_last20, na.rm = TRUE), 1), fmt(stats::median(x$osc_last50, na.rm = TRUE), 0)))
  }
}

cat("\n== check 1: two-staged seeds 1-3, original scale, at iteration 100 against the F0 sweep (stage 200) ==\n")
f0f <- list.files(file.path(RAW2S_SRC, "F0"), "^F0_task[0-9]+\\.rds$", full.names = TRUE)
if (!length(f0f)) cat("  F0 task files not found; skipped\n") else {
  f0 <- unlist(lapply(f0f, function(f) readRDS(f)$cells), recursive = FALSE)
  any_c <- FALSE
  for (r in runs[["200"]]) if (r$init == "two-staged" && r$transform == "raw" && r$seed_idx <= 3) {
    c0 <- Filter(function(z) z$scale == "raw" && z$init == "two-staged" && z$bwK == r$bwK && z$seed_idx == r$seed_idx, f0)
    if (!length(c0)) next; c0 <- c0[[1]]; any_c <- TRUE
    msg <- if (!isTRUE(r$ok) && !isTRUE(c0$ok)) "both failed (as expected)" else
      if (isTRUE(r$ok) != isTRUE(c0$ok)) "DIFFERENT: one failed, the other ran" else
      if (is.null(r$partitions) || ncol(r$partitions) < 100L) sprintf("stopped before 100 (%s at %d)", r$stop_reason, r$iter) else {
        a <- ari(r$partitions[, 100L], c0$g); sprintf("aRI at iteration 100 = %.3f%s", a, if (a > 0.999) " (identical)" else "  << NOT IDENTICAL") }
    cat(sprintf("  bwK %g seed %d: %s\n", r$bwK, r$seed_idx, msg))
  }
  if (!any_c) cat("  no matching cell (the sweep ran seeds 1-3 at bwK 10 and 12 on the original scale)\n")
}

cat("\n== check 2: every continued fit started from its stage-200 partition ==\n")
cont <- Filter(function(r) isTRUE(r$continued) && !is.null(r$g_stage200), runs[["500"]])
if (!length(cont)) cat("  no continued fits yet\n") else {
  okc <- vapply(cont, function(r) !is.null(r$g_start) && ari(r$g_start, r$g_stage200) > 0.999, TRUE)
  cat(sprintf("  %d of %d continued fits start from their stage-200 partition%s\n", sum(okc), length(okc), if (all(okc)) "" else "  << CHECK"))
  st <- tab$straddle_cycle[tab$maxIter == 500 & !is.na(tab$straddle_cycle)]
  cat(sprintf("  cycles straddling the stage boundary (a stage-500 partition equal to a stage-200 one): %d%s\n", length(st),
              if (length(st)) sprintf(" (lengths %s)", paste(st, collapse = ", ")) else ""))
}

cat("\n== check 3: std (affine rescaling) against raw, same seed and bandwidth ==\n")
if (!all(c("raw", "std") %in% r0$transforms)) cat("  not run (needs both transforms)\n") else {
  for (m in names(runs)) for (r in runs[[m]]) if (r$transform == "std") {
    p <- runs[[m]][[sub("^RAW2S_std_", "RAW2S_raw_", r$name)]]; if (is.null(p)) next
    msg <- if (isTRUE(r$ok) != isTRUE(p$ok)) "DIFFERENT: one failed, the other ran" else if (!isTRUE(r$ok)) "both failed" else {
      hr <- history(p); hs <- history(r)
      same <- !is.null(hr) && !is.null(hs) && ncol(hr) == ncol(hs) && all(vapply(seq_len(ncol(hr)), function(j) ari(hr[, j], hs[, j]) > 0.999, TRUE))
      sprintf("%s; iterations %d vs %d; J ratio raw/std %.4f (IQR^2 = %.4f)%s", if (same) "identical at every iteration" else "NOT identical",
              p$iter, r$iter, p$J / r$J, r$transform_params$scale^2, if (same) "" else "  << CHECK") }
    cat(sprintf("  max %s, %s: %s\n", m, sub("^RAW2S_std_", "", r$name), msg))
  }
}

cat("\n== check 4: against the main run's best original-scale start (F1) and the reported log-scale partition ==\n")
f1dir <- file.path(RAW2S_SRC, "F1")
fr <- lapply(list.files(f1dir, "^F1_raw_M1_start[0-9]+\\.rds$", full.names = TRUE), readRDS)
fr <- Filter(function(z) isTRUE(z$ok) && length(z$g) == n, fr)
g_raw <- if (length(fr)) fr[[which.min(vapply(fr, function(z) as.numeric(z$J), 0))]]$g else NULL
wl <- Filter(function(r) r$init == "warm-log" && !is.null(r$g_start), runs[["200"]])
g_log <- if (length(wl)) wl[[1]]$g_start else NULL
if (is.null(g_raw)) cat("  F1 original-scale starts not found; F1 comparison skipped\n") else
  cat(sprintf("  F1 best original-scale start: sizes %s, J %.1f\n", paste(table(g_raw), collapse = "/"),
              min(vapply(fr, function(z) as.numeric(z$J), 0))))
if (is.null(g_log)) cat("  no warm-log cell; log-partition comparison skipped\n") else
  cat(sprintf("  reported log-scale partition (the warm-log start): sizes %s, J %.4f\n", paste(table(g_log), collapse = "/"), wl[[1]]$warm$J))
for (m in names(runs)) for (i in seq_len(nrow(groups))) {
  gq <- groups[i, ]; if (gq$maxIter != as.integer(m)) next
  rr <- Filter(function(r) isTRUE(r$ok) && r$transform == gq$transform && r$bwK == gq$bwK && r$init == gq$start && length(r$g) == n, runs[[m]])
  if (!length(rr)) next
  best <- rr[[which.min(vapply(rr, `[[`, 0, "J"))]]
  cat(sprintf("  max %s | %-4s | bwK %2g | %-10s best of %d: J %.4f, sizes %s | aRI with F1 best %s, with the log partition %s%s\n",
              m, gq$transform, gq$bwK, gq$start, length(rr), best$J, paste(best$sizes, collapse = "/"),
              if (is.null(g_raw)) "-" else fmt(ari(best$g, g_raw), 3), if (is.null(g_log)) "-" else fmt(ari(best$g, g_log), 3),
              if (gq$start != "two-staged") sprintf(" | moved from its start: aRI(start, final) %s", fmt(ari(best$g_start, best$g), 3)) else ""))
}

cat("\n== fits still at the cap: is J still falling, or oscillating? (last 50 iterations) ==\n")
cat("  J here is the package's trace value after each reassignment (before the class fits are\n")
cat("  re-estimated), so it sits below the final J, which is taken after re-estimation.\n")
for (m in names(runs)) for (i in seq_len(nrow(groups))) {
  gq <- groups[i, ]; if (gq$maxIter != as.integer(m)) next
  rr <- Filter(function(r) isTRUE(r$ok) && r$transform == gq$transform && r$bwK == gq$bwK && r$init == gq$start &&
                 identical(r$stop_reason, "maxIter") && !is.null(r$objective_trace), runs[[m]])
  if (!length(rr)) next
  jj_of <- function(r) { v <- r$objective_trace$J_after; utils::tail(v[is.finite(v)], 50) }
  fall <- vapply(rr, function(r) { jj <- jj_of(r); (jj[length(jj)] - jj[1]) / jj[1] }, 0)
  turn <- vapply(rr, function(r) { jj <- jj_of(r); sum(diff(sign(diff(jj))) != 0) }, 0)
  cat(sprintf("  max %s | %-4s | bwK %2g | %-10s: %d at the cap | change in J over the last 50 iterations, median %+.2f%% (%+.2f%% to %+.2f%%) | direction changes median %g\n",
              m, gq$transform, gq$bwK, gq$start, length(rr), 100 * stats::median(fall), 100 * min(fall), 100 * max(fall), stats::median(turn)))
}

cat("\n== every fit ==\n"); print(tab, row.names = FALSE, digits = 4)
saveRDS(list(table = tab, root = ROOT, nseed = r0$nseed, bwks = r0$bwks, warms = r0$warms, transforms = r0$transforms, date = Sys.Date()),
        file.path(ROOT, "RAW2S.rds"))
cat(sprintf("\nwrote %s\n", file.path(ROOT, "RAW2S.rds")))
