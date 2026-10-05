# Combine block E -- estimation accuracy on the TRUE partition.
#
#   Rscript R/combine_estim.R <resultsdir> [n_reps]
#
# This block judges the estimators more directly than any other, because it
# removes the confound of the clustering loop: on the heavy-tailed arms
# CLUSTERING is dominated by initialisation failure, so an estimator judged
# through a clustering endpoint is being judged through a failure mode that has
# nothing to do with it. Fitting under the true partition strips all of that out.
#
# READ IN THIS ORDER, and stop at the first thing that is wrong:
#   1. completeness at 100 percent
#   2. the failure table
#   3. THE ANGLE RANGE. Angles are DEGREES and must lie in [0, 90]. An earlier
#      runner converted twice and produced 550 to 1450, which rbinds cleanly and
#      takes a perfectly well-behaved median. A range violation invalidates the
#      block outright and no amount of replication repairs it.
#   4. only then any accuracy number

source(file.path("config", "validate_config.R"))
source(file.path("R", "metrics.R"))

args   <- commandArgs(trailingOnly = TRUE)
args   <- args[args != "--args"]
resdir <- if (length(args) >= 1) args[1] else "results/estim"
n_reps <- if (length(args) >= 2) as.integer(args[2]) else NA_integer_

files <- list.files(resdir, pattern = "_task[0-9]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no task output in ", resdir, call. = FALSE)
raw <- read_task_files(files)
cfg <- load_config()

arms <- sort(unique(raw$arm))
ests <- sort(unique(raw$estimator))
cls  <- sort(unique(raw$class[!is.na(raw$class)]))
reps <- length(unique(raw$rep))
if (is.na(n_reps)) n_reps <- reps

cat(sprintf("\n%d task files, %d rows | block estim | arms %s | %d replicates\n",
            length(files), nrow(raw), paste(arms, collapse = ", "), reps))

expected <- length(arms) * n_reps * length(ests) * length(cls)
cat(sprintf("  %d rows expected (%d arms x %d reps x %d estimators x %d classes), %d present (%.1f%%)\n",
            expected, length(arms), n_reps, length(ests), length(cls),
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

# NOTE ON THE ARITHMETIC ABOVE. The estimator and class counts come from the ROWS
# PRESENT, so a level absent from every task cannot be detected here -- both
# sides of the ratio would shrink together. It catches a missing TASK, not a
# missing LEVEL. Check the estimator list printed above against config.
cat(sprintf("  estimators present: %s | classes present: %s\n",
            paste(ests, collapse = ", "), paste(cls, collapse = ", ")))

# NOT IMPLEMENTED is counted APART from real failures. E2 and E3 are package
# items 2 and 3 and are recorded deliberately rather than omitted, so that the
# gap stays visible and the completeness arithmetic stays right. Conflating the
# two would either hide a broken run inside an expected gap, or make an expected
# gap look like a broken run -- and with two of four estimators unbuilt the
# combined figure reads "50 percent failed", which is alarming and wrong.
is_ni   <- !is.na(raw$error) & grepl("^NOT IMPLEMENTED", raw$error)
is_fail <- !is.na(raw$error) & !is_ni
cat(sprintf("  %d of %d rows are a stated NOT IMPLEMENTED (%.1f%%)\n",
            sum(is_ni), nrow(raw), 100 * sum(is_ni) / nrow(raw)))
if (any(is_ni)) {
  tb <- table(raw$estimator[is_ni], substr(raw$error[is_ni], 1, 52))
  for (i in seq_len(nrow(tb))) for (j in seq_len(ncol(tb)))
    if (tb[i, j] > 0) cat(sprintf("    %-14s %3d  %s\n",
                                  rownames(tb)[i], tb[i, j], colnames(tb)[j]))
}
cat(sprintf("  %d of %d rows actually failed (%.1f%%)\n",
            sum(is_fail), nrow(raw), 100 * sum(is_fail) / nrow(raw)))
if (any(is_fail)) {
  tb <- table(raw$estimator[is_fail], substr(raw$error[is_fail], 1, 46))
  for (i in seq_len(nrow(tb))) for (j in seq_len(ncol(tb)))
    if (tb[i, j] > 0) cat(sprintf("    %-14s %3d  %s\n",
                                  rownames(tb)[i], tb[i, j], colnames(tb)[j]))
}

# --- 3. THE ANGLE RANGE GUARD, before any accuracy number is printed ----------
ang <- unlist(raw[, c("angle_1", "angle_2", "angle_sub")], use.names = FALSE)
ang <- ang[is.finite(ang)]
cat("\n== ANGLE RANGE (degrees; MUST lie in [0, 90]) ==\n")
if (!length(ang)) {
  cat("  *** NO FINITE ANGLES AT ALL. The block produced nothing readable. ***\n")
} else {
  cat(sprintf("  observed %.2f to %.2f over %d values\n",
              min(ang), max(ang), length(ang)))
  if (min(ang) < -1e-8 || max(ang) > 90 + 1e-8) {
    cat("  *** OUT OF RANGE. This is the double-conversion defect, not a result.\n")
    cat("  *** Every angle number below is meaningless. Fix run_estim.R first. ***\n")
  } else cat("  in range\n")
}

med  <- function(x) { x <- x[is.finite(x)]; if (length(x)) stats::median(x) else NA_real_ }
mcse <- function(x) { x <- x[is.finite(x)]; if (length(x) >= 5L) mcse_median(x) else NA_real_ }

arm_name <- function(a) {
  r <- cfg$arms[cfg$arms$arm == a, ]
  if (nrow(r) && !is.null(r$label)) as.character(r$label)[1] else paste("arm", a)
}

# --- 4. the endpoints ---------------------------------------------------------
cat("\n== ESTIMATION ENDPOINTS, median over replicates and classes ==\n")
cat("   mean_err  relative L2 error of mu_hat -- the NON-robust local-linear\n")
cat("             mean\n")
cat("   angle_m   degrees between phi_hat_m and the true phi_m\n")
cat("   angle_sub largest principal angle between the estimated and true spans\n")
cat("   cancor    mean canonical correlation of the scores\n\n")

for (a in arms) {
  cat(sprintf("  arm %d  %s\n", a, arm_name(a)))
  cat("     estimator     mean_err +/- mcse    angle_1   angle_2  angle_sub   cancor   n_ok\n")
  for (e in ests) {
    s <- raw[raw$arm == a & raw$estimator == e, ]
    ok <- sum(is.na(s$error))
    # Width on the MCSE field is load-bearing: sprintf pads NA to the field
    # width, but a width-less %.4f gives it two characters and shifts every
    # column to the right of it.
    cat(sprintf("     %-12s  %7.4f +/- %6.4f  %8.2f  %8.2f  %9.2f  %7.4f  %5d\n",
                e, med(s$mean_err), mcse(s$mean_err),
                med(s$angle_1), med(s$angle_2), med(s$angle_sub),
                med(s$cancor), ok))
  }
  cat("\n")
}

# --- the mean, per class ------------------------------------------------------
# Pooled over classes the mean error hides the case that matters: a single
# contaminated class drags its own mu_hat and leaves the other alone, so a pooled
# median can look healthy while one class is badly estimated.
cat("== MEAN-FUNCTION ERROR BY CLASS (unpooled) ==\n")
cat("     arm  estimator     ")
for (k in cls) cat(sprintf("  class %d ", k))
cat("\n")
for (a in arms) for (e in ests) {
  cat(sprintf("     %3d  %-12s  ", a, e))
  for (k in cls) {
    s <- raw[raw$arm == a & raw$estimator == e & raw$class == k, ]
    cat(sprintf("  %7.4f ", med(s$mean_err)))
  }
  cat("\n")
}

cat("\n   The mean is local-linear least squares and is NOT robust. If mean_err\n")
cat("   rises across the tail arms while the ANGLES hold, that is exactly the\n")
cat("   story to tell: Kendall's tau protects the eigenspace, and the mean is\n")
cat("   the part that still needs the robust mean option.\n")
