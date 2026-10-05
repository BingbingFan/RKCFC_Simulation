# READ STUDY F1 -- the common-pipeline refit.
#
#   Rscript R/read_realdata_f1.R [f1dir]
#
# READ IN THIS ORDER:
#
#   1. the run is REAL, not a smoke subsample, and no fit failed
#   2. the cohort
#   3. WITHIN EACH SCALE, the two methods compared on one footing. This is the
#      only comparison the submitted paper needed to make and did not.
#   4. all four fits placed on the RAW scale, which is the only honest way to ask
#      whether the log transform helps
#   5. the cluster sizes, because a 1399/6 split is not a finding
#
# WHY THE SUBMITTED RATIO WAS WITHDRAWN. The submitted 67.2 percent improvement divided an error
# measured on one scale by an error measured on another. It is not a small
# reporting problem: the ratio has no units and no meaning, and its SIZE is an
# artefact of which pairing was chosen. Panel 3 below shows how far the number
# moves when the pairing changes: the figure was not slightly off, it was never
# a quantity.

# Sources nothing on purpose: F1.rds carries everything this needs, so the panels
# can be read on a laptop from a downloaded results directory.

args  <- commandArgs(trailingOnly = TRUE)
args  <- args[args != "--args"]
f1dir <- if (length(args) >= 1) args[1] else file.path("results_real", "F1")

fp <- file.path(f1dir, "F1.rds")
if (!file.exists(fp)) stop("not found: ", fp, "\n  Has F1 finished?", call. = FALSE)
f1 <- readRDS(fp)
tb <- f1$table
d  <- f1$descriptives

fmt <- function(x, k = 4) if (is.finite(x)) sprintf(paste0("%.", k, "f"), x) else "-"

# --- 1. is this reportable at all ---------------------------------------------
cat("\n== F1: both methods under one pipeline ==\n")
if (isTRUE(any(tb$smoke))) {
  cat("\n  *** THIS IS A SMOKE RUN ON A SUBSAMPLE. Nothing below is reportable, and\n")
  cat("  *** no number from it may be reported.\n")
  cat("  *** Re-run F1 without the n_sub argument. ***\n")
}
bad <- !is.na(tb$error)
if (any(bad)) {
  cat(sprintf("\n  *** %d of %d fits FAILED ***\n", sum(bad), nrow(tb)))
  for (i in which(bad))
    cat(sprintf("      %s / %s: %s\n", tb$scale[i], tb$method[i],
                substr(tb$error[i], 1, 80)))
  cat("  A missing fit is not a missing row: the comparison it belonged to cannot\n")
  cat("  be made at all, so read the panels below as incomplete.\n")
}
cat(sprintf("\n  K = %d, subject-level bandwidth %.2f months, %d subjects\n",
            f1$K, f1$bw_months, d$n_subjects))
ch <- unique(tb$code_hash[!is.na(tb$code_hash)])
cat(sprintf("  code_hash %s\n", if (length(ch) == 1L) ch else
            paste("*** DISAGREES:", paste(ch, collapse = " vs "), "***")))

# --- 2. the cohort ----------------------------------------------------------
cat("\nCOHORT\n")
cat(sprintf("  %d subjects of %d in the source, %d excluded for no baseline record\n",
            d$n_subjects, d$n_source, d$excluded))
cat(sprintf("  %d CEA observations\n", d$n_obs))
cat(sprintf("  measurements per subject: median %.0f, IQR %.0f-%.0f, range %.0f-%.0f\n",
            d$ni_median, d$ni_iqr[1], d$ni_iqr[2], d$ni_range[1], d$ni_range[2]))
cat(sprintf("  observation times: median %.1f, IQR %.1f-%.1f, range %.0f-%.0f months\n",
            d$time_median, d$time_iqr[1], d$time_iqr[2],
            d$time_range[1], d$time_range[2]))

# --- 3. WITHIN SCALE: the comparison that was missing -------------------------
own <- tb[tb$metric_scale == tb$scale, , drop = FALSE]
cat("\n== PANEL 1. WITHIN EACH SCALE -- both methods, one footing ==\n")
cat(sprintf("  %-6s %-6s %10s %10s %7s %-16s %12s %9s\n", "scale", "method",
            "rre_in", "rre_out", "iter", "stopped", "sizes", "seconds"))
for (sc in c("raw", "log")) for (m in c("M1", "M2")) {
  r <- own[own$scale == sc & own$method == m, , drop = FALSE]
  if (!nrow(r)) next
  cat(sprintf("  %-6s %-6s %10s %10s %7s %-16s %12s %9.0f\n",
              if (m == "M1") sc else "", m, fmt(r$rre_in[1]), fmt(r$rre_out[1]),
              if (is.na(r$iter[1])) "-" else as.character(r$iter[1]),
              if (is.na(r$stop_reason[1])) "-" else r$stop_reason[1],
              if (is.na(r$n1[1])) "-" else sprintf("%d/%d", r$n1[1], r$n2[1]),
              r$seconds[1]))
}

imp <- function(a, b) if (is.finite(a) && is.finite(b) && b != 0)
                        100 * (b - a) / b else NA_real_
cat("\n  IMPROVEMENT OF M1 OVER M2, computed within a scale and nowhere else:\n")
imp_tab <- NULL
for (sc in c("raw", "log")) {
  a1 <- own$rre_in[own$scale == sc & own$method == "M1"]
  a2 <- own$rre_in[own$scale == sc & own$method == "M2"]
  o1 <- own$rre_out[own$scale == sc & own$method == "M1"]
  o2 <- own$rre_out[own$scale == sc & own$method == "M2"]
  if (!length(a1) || !length(a2)) next
  ii <- imp(a1[1], a2[1]); oo <- if (length(o1) && length(o2)) imp(o1[1], o2[1]) else NA_real_
  cat(sprintf("    %-6s in-sample %s%%   out-of-sample %s%%\n", sc,
              fmt(ii, 1), fmt(oo, 1)))
  imp_tab <- rbind(imp_tab, data.frame(scale = sc, improvement_in = ii,
                                       improvement_out = oo, stringsAsFactors = FALSE))
}
cat("  OUT-OF-SAMPLE IS THE ONE TO QUOTE. In-sample reconstruction error falls\n")
cat("  whenever a method fits harder, so it rewards flexibility rather than\n")
cat("  accuracy; the held-out figure is the one to report.\n")

# --- 4. all four fits on the RAW scale ----------------------------------------
raw_panel <- tb[grepl("^raw", tb$metric_scale), , drop = FALSE]
cat("\n== PANEL 2. EVERY FIT SCORED ON THE RAW SCALE ==\n")
cat("  Log-scale fits are back-transformed with exp() before scoring, so all four\n")
cat("  numbers are errors in the same units. This is the only place a raw fit and\n")
cat("  a log fit may be compared, and the comparison it supports is whether the\n")
cat("  transform helps -- not which method is better.\n\n")
cat(sprintf("  %-8s %-6s %-26s %10s\n", "fitted", "method", "scored on", "rre_in"))
for (i in seq_len(nrow(raw_panel)))
  cat(sprintf("  %-8s %-6s %-26s %10s\n", raw_panel$scale[i], raw_panel$method[i],
              raw_panel$metric_scale[i], fmt(raw_panel$rre_in[i])))
best <- raw_panel[which.min(replace(raw_panel$rre_in,
                                    !is.finite(raw_panel$rre_in), Inf)), , drop = FALSE]
if (nrow(best))
  cat(sprintf("\n  lowest raw-scale error: %s fitted on the %s scale (%s)\n",
              best$method[1], best$scale[1], fmt(best$rre_in[1])))

# --- 5. the defect, demonstrated ----------------------------------------------
cat("\n== PANEL 3. WHY THE SUBMITTED 67.2 PERCENT IS NOT A QUANTITY ==\n")
cat("  Below are the improvement figures obtainable by pairing a method on one\n")
cat("  scale against the other method on the other scale -- the arithmetic the\n")
cat("  submitted analysis performed. NONE OF THESE IS A RESULT. They are printed\n")
cat("  only to show that the figure moves with an arbitrary\n")
cat("  choice, which is a stronger and more honest concession than calling it\n")
cat("  an error of a few percentage points.\n\n")
gv <- function(sc, m) { v <- own$rre_in[own$scale == sc & own$method == m]
                        if (length(v)) v[1] else NA_real_ }
mix <- NULL
for (p in list(c("raw", "M1", "log", "M2"), c("log", "M1", "raw", "M2"),
               c("raw", "M2", "log", "M1"), c("log", "M2", "raw", "M1"))) {
  v <- imp(gv(p[1], p[2]), gv(p[3], p[4]))
  cat(sprintf("    %s/%s against %s/%s : %s%%\n", p[1], p[2], p[3], p[4], fmt(v, 1)))
  mix <- c(mix, v)
}
mix <- mix[is.finite(mix)]
if (length(mix) > 1L)
  cat(sprintf("\n  Range across four arbitrary pairings: %.1f to %.1f percent.\n",
              min(mix), max(mix)))
cat("  A quantity whose value depends on which of two units you divide by is not\n")
cat("  a measurement of anything.\n")

# --- 6. the partitions --------------------------------------------------------
cat("\n== PANEL 4. THE PARTITIONS ==\n")
cat(sprintf("  %-10s %14s %14s\n", "fit", "sizes", "smallest"))
disagree <- character(0)
for (k in names(f1$fits)) {
  g <- f1$fits[[k]]
  sz <- as.integer(table(factor(g, levels = seq_len(f1$K))))
  sm <- min(sz) / sum(sz)
  cat(sprintf("  %-10s %14s %13.1f%%%s\n", k, paste(sz, collapse = "/"), 100 * sm,
              if (sm < 0.05) "   *** near-degenerate ***" else ""))

  # THE SAME FACT, RECORDED TWICE. `n1`/`n2` in the table and the stored
  # partition are written by separate lines of the runner, so a disagreement
  # means the object is internally inconsistent -- a half-written file, or two
  # runs merged. Sorted, because the class labels carry no order.
  p <- strsplit(k, "_", fixed = TRUE)[[1]]
  r <- own[own$scale == p[1] & own$method == p[2], , drop = FALSE]
  if (nrow(r) && !is.na(r$n1[1]) &&
      !identical(sort(sz), sort(as.integer(c(r$n1[1], r$n2[1])))))
    disagree <- c(disagree, sprintf("%s: table says %d/%d, partition says %s",
                                    k, r$n1[1], r$n2[1], paste(sz, collapse = "/")))
}
if (length(disagree)) {
  cat("\n  *** THE TABLE AND THE STORED PARTITIONS DISAGREE ON CLUSTER SIZES ***\n")
  for (s in disagree) cat(sprintf("      %s\n", s))
  cat("  *** These are two records of one fit. F1.rds is internally inconsistent,\n")
  cat("  *** so it was half-written or two runs were merged. Re-run F1; do not\n")
  cat("  *** reconcile them by hand, and do not let F2/F3/F4 read this file. ***\n")
}
cat("\n  A partition whose smaller class holds under 5 percent of subjects is a\n")
cat("  detected outlier group, not a subtype, and F2's hazard ratio on it will\n")
cat("  have an interval too wide to interpret. Say which it is in the Results.\n")

cat("\n== WHAT THIS LICENSES ==\n")
if (!is.null(imp_tab)) {
  for (i in seq_len(nrow(imp_tab)))
    cat(sprintf("   On the %s scale, with both methods given identical data and identical\n   metrics, the robust method reduced held-out reconstruction error by %s%%.\n",
                imp_tab$scale[i], fmt(imp_tab$improvement_out[i], 1)))
}
cat("   Report the two scales AS A PAIR. Reporting whichever is larger, without\n")
cat("   the other, reintroduces exactly that selection.\n")
