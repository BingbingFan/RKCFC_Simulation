# COMBINE STUDY F4 -- clustering stability. 200 task files, one answer.
#
#   Rscript R/combine_realdata_f4.R [f4dir] [f1dir]
#
# READ IN THIS ORDER, and stop at the first thing that is wrong:
#
#   1. completeness at 100 percent, and NO duplicate items
#   2. the failure table
#   3. REPRODUCIBILITY. If 100 random starts do not agree with each other, the
#      partition F1 reports is one draw from a lottery, and every number that
#      depends on it -- F2's hazard ratio, F3's models, the Abstract -- is
#      describing a particular random seed. Nothing below matters until this has
#      been read.
#   4. only then Jaccard, margin and K selection
#
# WHAT NONE OF THIS MEASURES. There are no true labels on this cohort, so no
# quantity here is an ACCURACY. Reproducibility, survival under resampling and
# assignment confidence are the three things answerable without truth; accuracy
# lives in the simulation and must not be claimed from this file.
#
# WHY THIS CANNOT USE read_task_files(). That helper rbinds data frames and
# checks provenance on columns. F4 tasks write a LIST -- items of three different
# shapes, one carrying a length-n vector -- so the same discipline is applied
# here by hand: the scalar stamps must agree across every file, or the directory
# mixes runs and the item count cannot see it.

# Deliberately sources NOTHING. Every other combiner needs the config tables to
# name its arms; this one reads only what the task files carry, so it runs
# against a downloaded results directory with no code beside it.

args  <- commandArgs(trailingOnly = TRUE)
args  <- args[args != "--args"]
f4dir <- if (length(args) >= 1) args[1] else file.path("results_real", "F4")
f1dir <- if (length(args) >= 2) args[2] else file.path("results_real", "F1")

if (!requireNamespace("mclust", quietly = TRUE))
  stop("package `mclust` is required (adjustedRandIndex).", call. = FALSE)

files <- list.files(f4dir, pattern = "^F4_task[0-9]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no F4 task output in ", f4dir, call. = FALSE)
xs <- lapply(files, readRDS)

# --- 1. provenance and completeness -------------------------------------------
for (fld in c("n", "K", "n_tasks", "n_start", "n_sub", "code_hash")) {
  v <- vapply(xs, function(z) if (is.null(z[[fld]])) NA_character_ else
                              as.character(z[[fld]])[1], character(1))
  u <- unique(v)
  if (length(u) > 1L) {
    tb <- table(v)
    stop(sprintf(paste("task files disagree on `%s`: %s.\n",
                       " This directory mixes runs made under different assumptions,",
                       "and the item count cannot see it.\n  Counts: %s\n",
                       " Clear %s and re-run the array."),
                 fld, paste(u, collapse = " vs "),
                 paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = ", "),
                 f4dir), call. = FALSE)
  }
}

n <- xs[[1]]$n; K <- xs[[1]]$K
n_tasks <- xs[[1]]$n_tasks; n_start <- xs[[1]]$n_start; n_sub <- xs[[1]]$n_sub

items <- unlist(lapply(xs, function(z) z$items), recursive = FALSE)
kind  <- vapply(items, function(z) as.character(z$kind)[1], character(1))
idx   <- vapply(items, function(z) as.integer(z$idx)[1], integer(1))
err   <- vapply(items, function(z) if (is.null(z$error)) NA_character_ else
                                   as.character(z$error)[1], character(1))

# The work list, built EXACTLY as the runner builds it, so the two cannot drift.
work <- rbind(
  data.frame(kind = "restart",   idx = seq_len(n_start), stringsAsFactors = FALSE),
  data.frame(kind = "subsample", idx = seq_len(n_sub),   stringsAsFactors = FALSE),
  data.frame(kind = "kselect",   idx = 2:5,              stringsAsFactors = FALSE))
work$task <- (seq_len(nrow(work)) - 1L) %% n_tasks + 1L

have <- paste(kind, idx); want <- paste(work$kind, work$idx)
miss <- work[!want %in% have, , drop = FALSE]
dupe <- sum(duplicated(have))

cat(sprintf("\n== F4: clustering stability | n = %d, K = %d ==\n", n, K))
cat(sprintf("  %d task files of %d | %d items of %d (%.1f%%)\n",
            length(files), n_tasks, length(items), nrow(work),
            100 * length(items) / nrow(work)))
if (nrow(miss)) {
  cat(sprintf("  *** INCOMPLETE: %d items missing. Do not quote any number below. ***\n",
              nrow(miss)))
  tk <- sort(unique(miss$task))
  cat(sprintf("     missing tasks: %s\n", paste(tk, collapse = ",")))
  cat(sprintf("     re-run each with n_tasks = %d (tasks %s)\n",
              n_tasks, paste(tk, collapse = ",")))
  cat("     n_tasks is not optional: without it a re-run task slices a different\n")
  cat("     work list and writes the wrong items.\n")
}
if (dupe)
  cat(sprintf("  *** %d DUPLICATE items -- stale task files. Clear %s and re-run. ***\n",
              dupe, f4dir))

nfail <- sum(!is.na(err))
cat(sprintf("  %d of %d fits failed (%.1f%%)\n", nfail, length(items),
            100 * nfail / max(1L, length(items))))
if (nfail) {
  tb <- table(kind[!is.na(err)], substr(err[!is.na(err)], 1, 46))
  for (i in seq_len(nrow(tb))) for (j in seq_len(ncol(tb)))
    if (tb[i, j] > 0) cat(sprintf("    %-10s %3d  %s\n", rownames(tb)[i], tb[i, j],
                                  colnames(tb)[j]))
}

# --- the reference partition ---------------------------------------------------
# Jaccard needs a partition to be stable ABOUT. The right reference is the one
# the paper reports -- F1's raw_M1 -- because the question is whether THOSE
# clusters survive resampling, not whether some other partition would. Using the
# modal restart instead is a fallback, and it is labelled as one, because it
# quietly changes the question being answered.
# file.exists() first, rather than letting readRDS warn: R's own "cannot open
# compressed file" goes to stderr and reads like a fault in this script.
f1p <- file.path(f1dir, "F1.rds")
f1  <- if (file.exists(f1p)) tryCatch(readRDS(f1p), error = function(e) NULL) else NULL
ref <- NULL; ref_src <- "none"
if (is.null(f1)) {
  cat(sprintf("\n  %s is not readable -- F4 does not need it, but the Jaccard\n  reference does. Falling back to the modal restart, below.\n", f1p))
} else if (isTRUE(any(f1$table$smoke))) {
  cat("\n  F1.rds is a SMOKE RUN -- refusing to use it as the reference.\n")
} else {
  # The primary F1 was run under, if it recorded one. F4's restarts and
  # subsamples were fitted on the RAW scale, so a log-scale primary cannot be
  # the Jaccard reference here without conflating scale with stability: in
  # that case fall back to raw_M1 and say so, rather than compare across scales.
  # Task files written before 2026-09-10 carry no `scale` field and were raw.
  f4_scale <- if (is.null(xs[[1]]$scale)) "raw" else as.character(xs[[1]]$scale)
  # An F1.rds written before the `primary` field defaults to the M1 fit on F4's
  # own scale -- the one reference that cannot cross scales.
  prim <- if (!is.null(f1$primary)) f1$primary else paste0(f4_scale, "_M1")
  if (sub("_.*$", "", prim) != f4_scale) {
    cat(sprintf("\n  F1's primary is `%s`, but F4's refits are %s-scale. Using %s_M1 as the\n  Jaccard reference; a %s-scale primary needs F4 re-run on that scale.\n",
                prim, f4_scale, f4_scale, sub("_.*$", "", prim)))
    prim <- paste0(f4_scale, "_M1")
  }
  if (!is.null(f1$fits[[prim]]) && length(f1$fits[[prim]]) == n) {
    ref <- as.integer(f1$fits[[prim]]); ref_src <- sprintf("F1 %s, the reported partition", prim)
  } else {
    cat(sprintf("\n  F1.rds has no usable %s partition of length n.\n", prim))
  }
}

# --- 3. REPRODUCIBILITY --------------------------------------------------------
ok_g <- function(z) !is.null(z$g) && length(z$g) == n && !anyNA(z$g)
rs <- items[kind == "restart" & is.na(err)]
rs <- rs[vapply(rs, ok_g, logical(1))]

cat(sprintf("\n== REPRODUCIBILITY: %d usable random starts of %d ==\n",
            length(rs), n_start))
sol_tab <- NULL; ari_pair <- numeric(0); p_best <- NA_real_; p_ref <- NA_real_
if (length(rs) < 2L) {
  cat("  *** fewer than two usable restarts: reproducibility cannot be assessed ***\n")
} else {
  G <- do.call(rbind, lapply(rs, function(z) as.integer(z$g)))
  J <- vapply(rs, function(z) if (is.null(z$J)) NA_real_ else as.numeric(z$J)[1], numeric(1))
  sr <- vapply(rs, function(z) if (is.null(z$stop_reason)) NA_character_ else
                               as.character(z$stop_reason)[1], character(1))

  # Solutions grouped EXACTLY, by relabelling each partition in order of first
  # appearance. ARI == 1 gives the same grouping but is not transitive in
  # general; an exact canonical form is, so the classes below are well defined.
  canon <- apply(G, 1, function(g) paste(match(g, unique(g)), collapse = ","))
  cls   <- match(canon, unique(canon))

  m <- nrow(G); ari_pair <- numeric(m * (m - 1) / 2); p <- 0L
  for (i in seq_len(m - 1L)) for (j in (i + 1L):m) {
    p <- p + 1L; ari_pair[p] <- mclust::adjustedRandIndex(G[i, ], G[j, ])
  }
  qa <- stats::quantile(ari_pair, c(0, .25, .5, .75))
  cat(sprintf("  pairwise ARI over %d pairs: median %.4f, IQR %.4f-%.4f, min %.4f\n",
              length(ari_pair), qa[3], qa[2], qa[4], qa[1]))
  cat(sprintf("  %.1f%% of pairs identical up to relabelling\n",
              100 * mean(ari_pair >= 1 - 1e-9)))

  Jb <- min(J, na.rm = TRUE)
  p_best <- mean(J <= Jb * (1 + 1e-6), na.rm = TRUE)
  cat(sprintf("\n  %d DISTINCT solutions among %d starts\n", length(unique(cls)), m))
  cat(sprintf("  %-4s %7s %8s %14s %14s %s\n", "sol", "starts", "share", "median J",
              "best J", "sizes"))
  sol_tab <- do.call(rbind, lapply(sort(unique(cls)), function(c0) {
    w <- which(cls == c0)
    data.frame(solution = c0, starts = length(w), share = length(w) / m,
               J_median = stats::median(J[w], na.rm = TRUE),
               J_best = suppressWarnings(min(J[w], na.rm = TRUE)),
               sizes = paste(as.integer(table(factor(G[w[1], ], levels = seq_len(K)))),
                             collapse = "/"),
               stringsAsFactors = FALSE) }))
  sol_tab <- sol_tab[order(-sol_tab$starts), , drop = FALSE]
  # THE BEST-J ROW IS ALWAYS SHOWN, even when it is not among the most frequent.
  # The first version printed the top six BY START COUNT and nothing else, and on
  # the real run that hid the best solution entirely: the winning objective was
  # 402360 against the modal 577064, and its row was inside "... and 1 more".
  # Since the rule this file recommends is "report the best-J partition", hiding
  # that partition made the recommendation unusable.
  ib   <- which.min(sol_tab$J_best)
  show <- unique(c(seq_len(min(6L, nrow(sol_tab))), ib))
  for (i in show)
    cat(sprintf("  %-4d %7d %7.1f%% %14.4f %14.4f %s%s\n", sol_tab$solution[i],
                sol_tab$starts[i], 100 * sol_tab$share[i], sol_tab$J_median[i],
                sol_tab$J_best[i], sol_tab$sizes[i],
                if (i == ib) "   <- BEST J" else ""))
  if (nrow(sol_tab) > length(show))
    cat(sprintf("  ... and %d more (all with a worse objective than the row marked BEST J)\n",
                nrow(sol_tab) - length(show)))

  cat(sprintf("\n  prop_best: %.0f%% of starts reach the best objective found\n",
              100 * p_best))
  cat("  The simulation puts this at 0.25 on the heavy-tailed arm -- three starts\n")
  cat("  in four miss the best J -- so a low value here is expected behaviour of\n")
  cat("  the algorithm being measured, not evidence the code is wrong. What it\n")
  cat("  means for the paper is that a SINGLE start is not a defensible protocol.\n")

  # Does the MOST FREQUENT solution also have the best objective? If not, the
  # partition a single run most often returns is not the one the criterion
  # prefers, and multi-start selection by J is doing real work rather than
  # tidying.
  top <- sol_tab[1, ]
  if (is.finite(top$J_best) && top$J_best > Jb * (1 + 1e-6)) {
    cat(sprintf("  *** the MODAL solution (%.0f%% of starts) does NOT attain the best J\n",
                100 * top$share))
    cat(sprintf("      (%.4f against %.4f). Report the best-J partition, not the\n",
                top$J_best, Jb))
    cat("      most common one, and say in the Methods which rule was used. ***\n")
  }

  cat("\n  how the restarts ended:\n")
  for (nm in names(sort(table(sr), decreasing = TRUE)))
    cat(sprintf("    %-16s %d\n", nm, sum(sr == nm, na.rm = TRUE)))

  if (!is.null(ref)) {
    a_ref <- apply(G, 1, function(g) mclust::adjustedRandIndex(ref, g))
    p_ref <- mean(a_ref >= 1 - 1e-9)
    cat(sprintf("\n  AGAINST THE REPORTED PARTITION (%s)\n", ref_src))
    cat(sprintf("    median ARI %.4f; %.0f%% of random starts recover it exactly\n",
                stats::median(a_ref), 100 * p_ref))
    cat("    This is not circular: the restarts are independent fits, and the\n")
    cat("    question is how often a start lands on the answer the paper reports.\n")
  }
}

if (is.null(ref) && length(rs)) {
  G0 <- do.call(rbind, lapply(rs, function(z) as.integer(z$g)))
  cn <- apply(G0, 1, function(g) paste(match(g, unique(g)), collapse = ","))
  ref <- as.integer(G0[which.max(as.integer(table(cn)[cn])), ])
  ref_src <- "the MODAL restart -- F1 was not usable, so this is a fallback"
  cat(sprintf("\n  reference for Jaccard: %s\n", ref_src))
}

# --- 4a. survival under resampling --------------------------------------------
ss <- items[kind == "subsample" & is.na(err)]
ss <- ss[vapply(ss, function(z) !is.null(z$g) && !is.null(z$keep) &&
                                length(z$g) == length(z$keep), logical(1))]
jac_tab <- NULL; ari_sub <- numeric(0)
cat(sprintf("\n== SURVIVAL UNDER RESAMPLING: %d usable subsamples of %d, at 80%% ==\n",
            length(ss), n_sub))
if (!length(ss) || is.null(ref)) {
  cat("  not computable (no usable subsamples, or no reference partition)\n")
} else {
  Kref <- max(ref)
  jac <- matrix(NA_real_, nrow = length(ss), ncol = Kref)
  ari_sub <- numeric(length(ss))
  for (b in seq_along(ss)) {
    keep <- as.integer(ss[[b]]$keep); gb <- as.integer(ss[[b]]$g); rk <- ref[keep]
    ari_sub[b] <- mclust::adjustedRandIndex(rk, gb)
    for (k in seq_len(Kref)) {
      A <- rk == k
      if (!any(A)) next
      # Hennig's clusterwise Jaccard: the reference cluster RESTRICTED to the
      # resampled subjects, against its best match among the refit's clusters.
      jac[b, k] <- max(vapply(sort(unique(gb)), function(j) {
        B <- gb == j; sum(A & B) / sum(A | B) }, numeric(1)))
    }
  }
  qs <- stats::quantile(ari_sub, c(.25, .5, .75))
  cat("  ARI of each refit against the reference, restricted to that refit's own\n")
  cat(sprintf("  subsample: median %.4f, IQR %.4f-%.4f\n", qs[2], qs[1], qs[3]))
  cat(sprintf("\n  %-8s %8s %10s %10s %12s %12s %12s\n", "cluster", "n_ref",
              "mean J", "median J", "J >= 0.85", "J >= 0.75", "J < 0.60"))
  jac_tab <- do.call(rbind, lapply(seq_len(Kref), function(k) {
    v <- jac[, k]; v <- v[is.finite(v)]
    data.frame(cluster = k, n_ref = sum(ref == k),
               jaccard_mean = mean(v), jaccard_median = stats::median(v),
               p_ge_085 = mean(v >= 0.85), p_ge_075 = mean(v >= 0.75),
               p_lt_060 = mean(v < 0.60), n_b = length(v), stringsAsFactors = FALSE) }))
  for (i in seq_len(nrow(jac_tab)))
    cat(sprintf("  %-8d %8d %10.3f %10.3f %11.0f%% %11.0f%% %11.0f%%\n",
                jac_tab$cluster[i], jac_tab$n_ref[i], jac_tab$jaccard_mean[i],
                jac_tab$jaccard_median[i], 100 * jac_tab$p_ge_085[i],
                100 * jac_tab$p_ge_075[i], 100 * jac_tab$p_lt_060[i]))
  cat("\n  Hennig's reading, and it is a convention rather than a test: a mean\n")
  cat("  Jaccard of 0.75 or more is a stable cluster and 0.85 or more a highly stable\n")
  cat("  one; 0.60-0.75 indicates a pattern in the data, and below 0.60 the cluster\n")
  cat("  should not be interpreted at all (Hennig, 2007).\n")
  cat(sprintf("  Reference: %s\n", ref_src))
}

# --- 4b. how confident is each assignment -------------------------------------
mi  <- which(kind == "restart" & idx == 1L)
mrg <- if (length(mi)) items[[mi[1]]]$margin else NULL
gm  <- if (length(mi)) items[[mi[1]]]$g else NULL
cat("\n== ASSIGNMENT MARGIN: (second best - best) / best reconstruction error ==\n")
if (!is.numeric(mrg) || length(mrg) != n) {
  cat("  NOT AVAILABLE. The margin is computed on restart idx == 1 only, which\n")
  cat(sprintf("  lives in task %d.\n", (1L - 1L) %% n_tasks + 1L))
  cat("  If that task failed, re-run it alone:\n")
  cat(sprintf("    re-run task 1 with n_tasks = %d\n", n_tasks))
  cat("  A partition can be perfectly reproducible and still have most subjects\n")
  cat("  on the boundary, so this is not an optional panel.\n")
} else {
  v  <- mrg[is.finite(mrg)]
  qm <- stats::quantile(v, c(.1, .25, .5, .75, .9))
  cat(sprintf("  %d of %d subjects scored\n", length(v), n))
  cat(sprintf("  median %.3f, IQR %.3f-%.3f, deciles %.3f / %.3f\n",
              qm[3], qm[2], qm[4], qm[1], qm[5]))
  cat(sprintf("  margin < 0.05: %.0f%%   < 0.10: %.0f%%   < 0.25: %.0f%%\n",
              100 * mean(v < 0.05), 100 * mean(v < 0.10), 100 * mean(v < 0.25)))
  if (!is.null(gm) && length(gm) == n) {
    cat(sprintf("\n  %-8s %8s %10s %12s\n", "cluster", "n", "median", "share < 0.10"))
    for (k in sort(unique(as.integer(gm)))) {
      w <- mrg[as.integer(gm) == k]; w <- w[is.finite(w)]
      cat(sprintf("  %-8d %8d %10.3f %11.0f%%\n", k, length(w),
                  stats::median(w), 100 * mean(w < 0.10)))
    }
  }
  cat("\n  A margin of 0.10 means the second-choice class reconstructs that\n")
  cat("  subject only 10 percent worse. Those subjects are the ones that move\n")
  cat("  between solutions, and they are where the clinical story is weakest.\n")
}

# --- 4c. K selection, with the CORRECTED silhouette ---------------------------
ks <- items[kind == "kselect"]
cat("\n== K SELECTION, corrected silhouette ==\n")
k_tab <- NULL
if (!length(ks)) {
  cat("  no K-selection items present\n")
} else {
  k_tab <- do.call(rbind, lapply(ks, function(z) data.frame(
    K = as.integer(z$idx)[1],
    error = if (is.null(z$error)) NA_character_ else as.character(z$error)[1],
    silhouette = if (is.null(z$sil)) NA_real_ else as.numeric(z$sil)[1],
    sil_error = if (is.null(z$sil_error)) NA_character_ else as.character(z$sil_error)[1],
    J = if (is.null(z$J)) NA_real_ else as.numeric(z$J)[1],
    sizes = if (is.null(z$g)) NA_character_ else
            paste(as.integer(table(z$g)), collapse = "/"),
    stringsAsFactors = FALSE)))
  k_tab <- k_tab[order(k_tab$K), , drop = FALSE]
  cat(sprintf("  %-4s %14s %16s %s\n", "K", "silhouette", "objective", "sizes"))
  for (i in seq_len(nrow(k_tab)))
    cat(sprintf("  %-4d %14s %16s %s%s\n", k_tab$K[i],
                if (is.finite(k_tab$silhouette[i])) sprintf("%.4f", k_tab$silhouette[i]) else "-",
                if (is.finite(k_tab$J[i])) sprintf("%.4f", k_tab$J[i]) else "-",
                if (is.na(k_tab$sizes[i])) "-" else k_tab$sizes[i],
                if (!is.na(k_tab$error[i])) "   FAILED" else ""))
  # A blank silhouette on a fit that succeeded needs its reason printed: the
  # 2026-09-10 runs had one (a wrong call signature) and printed nothing.
  for (i in seq_len(nrow(k_tab)))
    if (is.na(k_tab$error[i]) && !is.finite(k_tab$silhouette[i]))
      cat(sprintf("  K = %d fitted but no silhouette: %s\n", k_tab$K[i],
                  if (is.na(k_tab$sil_error[i])) "no error recorded (runner predates sil_error)"
                  else k_tab$sil_error[i]))
  if (any(is.finite(k_tab$silhouette))) {
    sel <- k_tab$K[which.max(replace(k_tab$silhouette,
                                     !is.finite(k_tab$silhouette), -Inf))]
    cat(sprintf("\n  selected K = %d on the corrected coefficient%s\n", sel,
                if (sel != 2L) "   *** NOT the submitted K = 2 ***" else
                               "   (the submitted K)"))
  }
  cat("  The submitted coefficient smoothed the class means at a different\n")
  cat("  bandwidth from the model it scored, so this rerun is the number that\n")
  cat("  can be quoted; the old one cannot, whichever way it comes out.\n")
}

# --- the sentences this block licenses ----------------------------------------
cat("\n== WHAT THIS LICENSES ==\n")
if (length(ari_pair))
  cat(sprintf("   Across %d random initialisations the pairwise ARI had median %.3f,\n   and %.0f%% of starts reached the best objective.\n",
              length(rs), stats::median(ari_pair), 100 * p_best))
if (!is.null(jac_tab))
  cat(sprintf("   Under %d subsamples at 80%%, clusterwise Jaccard was %s.\n",
              length(ss), paste(sprintf("%.2f", jac_tab$jaccard_mean), collapse = " and ")))
if (is.numeric(mrg) && length(mrg) == n)
  cat(sprintf("   %.0f%% of subjects had an assignment margin below 0.10.\n",
              100 * mean(mrg[is.finite(mrg)] < 0.10)))
cat("   None of these is an accuracy. They measure STABILITY.\n")

saveRDS(list(n = n, K = K, n_tasks = n_tasks, n_start = n_start, n_sub = n_sub,
             n_items = length(items), n_expected = nrow(work),
             n_missing = nrow(miss), n_failed = nfail,
             reference = ref_src, solutions = sol_tab,
             ari_pairwise = ari_pair, prop_best = p_best, prop_recover_ref = p_ref,
             ari_subsample = ari_sub, jaccard = jac_tab,
             margin = if (is.numeric(mrg) && length(mrg) == n) mrg else NULL,
             margin_group = if (!is.null(gm) && length(gm) == n) as.integer(gm) else NULL,
             kselect = k_tab, code_hash = xs[[1]]$code_hash),
        file.path(f4dir, "F4.rds"))
cat(sprintf("\nwrote %s\n", file.path(f4dir, "F4.rds")))
