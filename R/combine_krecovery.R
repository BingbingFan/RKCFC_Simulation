# Combine block K -- K recovery and the silhouette criterion.
#
#   Rscript R/combine_krecovery.R <resultsdir> [n_reps]
#
# THE QUESTION THIS BLOCK EXISTS TO ANSWER. On the submitted real data both
# methods select K = 2, but KCFC scores 0.6694 against RKCFC's 0.0534 on the
# silhouette -- twelve and a half times higher, on the criterion this paper
# proposes. Simulation is the only place the true K is known, so it is the only
# place that gap can be explained rather than defended.
#
# THREE REPORTING RULES, ALL OF WHICH THE SUBMITTED PAPER BROKE:
#   1. the full DISTRIBUTION of selected K, never the mode
#   2. the silhouette VALUE, not only its argmax
#   3. BOTH criteria, silhouette and Calinski-Harabasz
#
# AND THE SILHOUETTE IS REPORTED TWICE ON PURPOSE. sil_manuscript is the
# package's curve-space variant -- squared, self-inclusive, distance to the
# cluster MEAN FUNCTION -- which is the quantity the submitted table quotes.
# sil_pairwise is Rousseeuw's on the scores. The squared form inflates the
# coefficient relative to a true distance, so the two
# are NOT interchangeable and the paper must say which it is quoting.
#
# K = 1 IS A CONVENTION, NOT A MEASUREMENT. The silhouette is undefined at
# K = 1; the 0 printed there is a convention. An argmax that can land on K = 1 is
# comparing a convention against a measurement, so selection below is taken over
# K >= 2 and the K = 1 column is shown separately.

source(file.path("config", "validate_config.R"))
source(file.path("R", "metrics.R"))

args   <- commandArgs(trailingOnly = TRUE)
args   <- args[args != "--args"]
resdir <- if (length(args) >= 1) args[1] else "results/krecovery"
n_reps <- if (length(args) >= 2) as.integer(args[2]) else NA_integer_

files <- list.files(resdir, pattern = "_task[0-9]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no task output in ", resdir, call. = FALSE)
raw <- read_task_files(files)
cfg <- load_config()

scens <- sort(unique(raw$scenario_id))
meths <- sort(unique(raw$method))
kfits <- sort(unique(raw$K_fit))
reps  <- length(unique(raw$rep))
if (is.na(n_reps)) n_reps <- reps

cat(sprintf("\n%d task files, %d rows | block krecovery | scenarios %s | %d replicates\n",
            length(files), nrow(raw), paste(scens, collapse = ", "), reps))

expected <- length(scens) * n_reps * length(meths) * length(kfits)
cat(sprintf("  %d rows expected (%d scenarios x %d reps x %d methods x %d K_fit), %d present (%.1f%%)\n",
            expected, length(scens), n_reps, length(meths), length(kfits),
            nrow(raw), 100 * nrow(raw) / expected))
if (nrow(raw) < expected)
  cat("  *** INCOMPLETE. Do not quote any number until this reads 100 percent. ***\n")
# EXCESS IS AS BAD AS SHORTFALL, and it is the quieter of the two. Every combiner
# globs _task*.rds and rbinds whatever it finds, so STALE FILES from an earlier
# run are merged in silently and the tables below mix two designs. Caught on
# 2026-09-06 in a smoke directory still holding the previous day's task files;
# the only symptom was a completeness figure ABOVE 100 percent, which reads as
# reassuring rather than alarming.
if (nrow(raw) > expected) {
  cat("  *** MORE ROWS THAN THE DESIGN HAS -- almost always STALE TASK FILES.\n")
  cat("  *** Check the file dates in this directory, clear it, and re-run.\n")
  cat("  *** Do not quote any number until this reads exactly 100 percent. ***\n")
}
cat(sprintf("  K_fit present: %s\n", paste(kfits, collapse = ", ")))

nfail <- sum(!is.na(raw$error))
cat(sprintf("  %d of %d fits failed (%.1f%%)\n", nfail, nrow(raw),
            100 * nfail / nrow(raw)))
if (nfail) {
  bad <- raw[!is.na(raw$error), ]
  tb <- table(paste0(bad$scenario_id, "/", bad$method, " K=", bad$K_fit),
              substr(bad$error, 1, 40))
  for (i in seq_len(nrow(tb))) for (j in seq_len(ncol(tb)))
    if (tb[i, j] > 0) cat(sprintf("    %-22s %3d  %s\n",
                                  rownames(tb)[i], tb[i, j], colnames(tb)[j]))
}

med  <- function(x) { x <- x[is.finite(x)]; if (length(x)) stats::median(x) else NA_real_ }
rate <- function(x) { x <- x[!is.na(x)]; if (length(x)) mean(x) else NA_real_ }

# Selection over K >= 2 ONLY, for the reason in the header.
select_over <- function(df, col) {
  d <- df[df$K_fit >= 2L & is.finite(df[[col]]), ]
  if (!nrow(d)) return(NA_integer_)
  d$K_fit[which.max(d[[col]])]
}

# --- 1. THE FULL DISTRIBUTION OF SELECTED K, never the mode -------------------
cat("\n== SELECTED K, FULL DISTRIBUTION over replicates (K >= 2 only) ==\n")
cat("   The mode alone hides whether the runner-up was 3 or 5, which is the\n")
cat("   whole question when the true K is 2.\n\n")

for (sc in scens) {
  kt <- unique(raw$K_true[raw$scenario_id == sc])[1]
  cat(sprintf("  %s   true K = %d\n", sc, kt))
  for (crit in c("sil_manuscript", "sil_pairwise", "CH")) {
    for (m in meths) {
      sel <- vapply(sort(unique(raw$rep)), function(rp) {
        d <- raw[raw$scenario_id == sc & raw$method == m & raw$rep == rp, ]
        as.integer(select_over(d, crit))
      }, integer(1))
      sel <- sel[!is.na(sel)]
      tb <- table(factor(sel, levels = kfits[kfits >= 2L]))
      cat(sprintf("     %-15s %-3s  ", crit, m))
      for (k in names(tb)) cat(sprintf("K=%s:%-3d ", k, tb[[k]]))
      # WHEN THE TRUTH IS K = 1 THERE IS NO CORRECTNESS RATE TO REPORT. Selection
      # runs over K >= 2, because the silhouette at K = 1 is a CONVENTION and an
      # argmax that can land there compares a convention against a measurement.
      # The true K is then outside the candidate set by construction, so a rate
      # would be 0.00 always -- which reads as catastrophic failure when it is
      # only a definitional exclusion. What K01 measures instead is how badly the
      # criteria OVER-SPLIT a single homogeneous population, which is the
      # distribution printed to the left.
      if (kt < 2L) cat("  correct n/a (true K = 1 is outside K >= 2)\n")
      else if (length(sel)) cat(sprintf("  correct %.2f\n", mean(sel == kt)))
      else cat("  correct n/a (no replicate selected)\n")
    }
  }
  cat("\n")
}

# --- 2. THE SILHOUETTE VALUES, which is what the real-data gap is about -------
cat("== CRITERION VALUES, median over replicates ==\n")
cat("   sil_manuscript is the SQUARED, self-inclusive curve-space variant the\n")
cat("   submitted table quotes; sil_pairwise is Rousseeuw's. They are not\n")
cat("   interchangeable -- quote one and name it.\n\n")

for (sc in scens) {
  kt <- unique(raw$K_true[raw$scenario_id == sc])[1]
  cat(sprintf("  %s   true K = %d\n", sc, kt))
  cat("     method  K_fit   sil_manu   sil_pair         CH      aRI   conv   n_ok\n")
  for (m in meths) for (k in kfits) {
    s <- raw[raw$scenario_id == sc & raw$method == m & raw$K_fit == k, ]
    if (!nrow(s)) next
    tag <- if (k == 1L) " (convention)" else ""
    cat(sprintf("     %-6s  %5d  %9.4f  %9.4f  %9.2f  %7.3f  %5.2f  %5d%s\n",
                m, k, med(s$sil_manuscript), med(s$sil_pairwise), med(s$CH),
                med(s$ari), rate(s$converged), sum(is.na(s$error)), tag))
  }
  cat("\n")
}

cat("   THE COMPARISON THAT MATTERS is KCFC's criterion value against RKCFC's\n")
cat("   AT THE SAME K, beside which of them actually recovers the true K. If\n")
cat("   KCFC scores higher while recovering K no better, the criterion rewards\n")
cat("   clusters separated in LEVEL over clusters separated in SHAPE -- which\n")
cat("   is the argument the paper needs.\n")
