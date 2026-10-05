# Combine one row of block S -- the settings studies.
#
#   Rscript R/combine_settings.R <resultsdir> [n_reps]
#
# ONE ROW AT A TIME. The four rows measure different things on different arms and
# a shared table would have to drop whatever the others do not have.
#
# READ IN THIS ORDER, and stop at the first thing that is wrong:
#   1. completeness at 100 percent
#   2. the failure table -- WHICH FOR THIS BLOCK INCLUDES THE STATED
#      "NOT IMPLEMENTED" LEVELS. Those are not failures of the run; they are
#      levels the package cannot do yet, recorded deliberately so the gap is
#      visible and the completeness arithmetic stays right. They are counted
#      separately below for exactly that reason.
#   3. convergence and iter -- a fit pinned at maxIter is an artefact of where
#      the loop was cut off, not a result
#   4. only then any accuracy number

source(file.path("config", "validate_config.R"))
source(file.path("R", "metrics.R"))

args   <- commandArgs(trailingOnly = TRUE)
args   <- args[args != "--args"]
resdir <- if (length(args) >= 1) args[1] else "results/S_init"
n_reps <- if (length(args) >= 2) as.integer(args[2]) else NA_integer_

files <- list.files(resdir, pattern = "_task[0-9]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no task output in ", resdir, call. = FALSE)
raw <- read_task_files(files)
cfg <- load_config()

row_id <- unique(raw$row)
if (length(row_id) != 1L)
  stop("results mix rows: ", paste(row_id, collapse = ", "), call. = FALSE)

arms <- sort(unique(raw$arm))
levs <- unique(raw$level)
reps <- length(unique(raw$rep))
if (is.na(n_reps)) n_reps <- reps

cat(sprintf("\n%d task files, %d rows | block S/%s | arms %s | %d replicates\n",
            length(files), nrow(raw), row_id, paste(arms, collapse = ", "), reps))

# ROWS PER CELL, NOT DISTINCT LEVEL STRINGS. The components row RENAMES its FVE
# levels with the M they selected -- run_settings.R sets level to "fve95(M=5)" --
# so a threshold that lands on different M in different replicates contributes
# TWO strings while still producing ONE row per replicate. Counting strings then
# inflates the denominator and reports a complete run as incomplete: on
# 2026-09-07 S_components read 88.9 percent with all 400 rows present, because
# fve95 resolved to M=5 on 48 replicates and M=6 on 2.
#
# The rows a single (arm, rep) cell produces is the quantity that is actually
# fixed by the design, so use that -- and assert it is constant, since a cell
# with too few rows is the real failure this check exists to catch.
per_cell  <- table(paste(raw$arm, raw$rep))
n_percell <- as.integer(names(sort(table(per_cell), decreasing = TRUE))[1])
expected  <- length(arms) * n_reps * n_percell
cat(sprintf("  %d rows expected (%d arms x %d reps x %d rows/cell), %d present (%.1f%%)\n",
            expected, length(arms), n_reps, n_percell, nrow(raw),
            100 * nrow(raw) / expected))
if (length(unique(per_cell)) > 1L) {
  cat(sprintf("  *** CELLS DISAGREE ON ROW COUNT: %s. A short cell is a task that\n",
              paste(sort(unique(per_cell)), collapse = ", ")))
  cat("  *** died part way, and the ratio above cannot see it. ***\n")
}
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

# The level count comes from the ROWS PRESENT, so a level missing from every
# task cannot be detected by the ratio above -- both sides shrink together. Print
# the list so it can be checked against LEVELS in run_settings.R by eye.
cat(sprintf("  levels present: %s\n", paste(levs, collapse = ", ")))

# NOT IMPLEMENTED is counted apart from real failures. Conflating them would
# either hide a broken run inside an expected gap, or make an expected gap look
# like a broken run.
is_ni   <- !is.na(raw$error) & grepl("^NOT IMPLEMENTED", raw$error)
is_fail <- !is.na(raw$error) & !is_ni
cat(sprintf("  %d of %d rows are a stated NOT IMPLEMENTED (%.1f%%)\n",
            sum(is_ni), nrow(raw), 100 * sum(is_ni) / nrow(raw)))
if (any(is_ni)) {
  tb <- table(raw$level[is_ni], substr(raw$error[is_ni], 1, 52))
  for (i in seq_len(nrow(tb))) for (j in seq_len(ncol(tb)))
    if (tb[i, j] > 0) cat(sprintf("    %-18s %3d  %s\n",
                                  rownames(tb)[i], tb[i, j], colnames(tb)[j]))
}
cat(sprintf("  %d of %d fits actually failed (%.1f%%)\n",
            sum(is_fail), nrow(raw), 100 * sum(is_fail) / nrow(raw)))
if (any(is_fail)) {
  tb <- table(raw$level[is_fail], substr(raw$error[is_fail], 1, 46))
  for (i in seq_len(nrow(tb))) for (j in seq_len(ncol(tb)))
    if (tb[i, j] > 0) cat(sprintf("    %-18s %3d  %s\n",
                                  rownames(tb)[i], tb[i, j], colnames(tb)[j]))
}

med  <- function(x) { x <- x[is.finite(x)]; if (length(x)) stats::median(x) else NA_real_ }
mcse <- function(x) { x <- x[is.finite(x)]; if (length(x) >= 5L) mcse_median(x) else NA_real_ }
rate <- function(x) { x <- x[!is.na(x)]; if (length(x)) mean(x) else NA_real_ }

arm_name <- function(a) {
  r <- cfg$arms[cfg$arms$arm == a, ]
  if (nrow(r) && !is.null(r$label)) as.character(r$label)[1] else paste("arm", a)
}

HEAD <- list(
  init = c("INITIALISATION. `oracle` bounds what any initialiser could buy and",
           "`pace_kmeans` reproduces the 0.1.x start deliberately. `two_staged` is the",
           "CONFIGURED default and is the row every other level is judged against."),
  components = c("RETAINED COMPONENTS. Read `conv` and `iter` FIRST: at nK = 4 the",
                 "univariate loop has been measured CAP-BOUND at 100 of 100, and a",
                 "cap-bound fit is unreportable however its accuracy reads."),
  smoothing = c("BANDWIDTH. The mean is local-linear least squares and is NOT",
                "robust; `mean_robust` is the robust mean option and is the level that",
                "would change that."),
  estimators = c("SCATTER ESTIMATORS. E1 is the proposal; E2 and E3 are package",
                 "items 2 and 3. Until those exist this row measures E1 alone and",
                 "is NOT a comparison."))

cat("\n== ")
cat(paste(HEAD[[row_id]], collapse = "\n   "))
cat(" ==\n\n")

for (a in arms) {
  cat(sprintf("  arm %d  %s\n", a, arm_name(a)))
  cat("     level               aRI med +/- mcse   CCR med +/- mcse     conv   iter   n_ok rre_true\n")
  for (lv in levs) {
    s <- raw[raw$arm == a & raw$level == lv, ]
    if (!nrow(s)) next
    ok <- sum(is.na(s$error))
    if (!ok) {
      msg <- unique(substr(s$error, 1, 40))[1]
      cat(sprintf("     %-18s %s\n", lv, msg))
      next
    }
    # WIDTHS ON THE MCSE FIELDS ARE LOAD-BEARING. sprintf pads NA to the field
    # width, but a width-less %.4f gives it two characters and every column to
    # the right shifts -- which on a table whose whole purpose is comparing a
    # number against twice an MCSE is a good way to misread one.
    cat(sprintf("     %-18s %6.3f +/- %6.4f   %6.3f +/- %6.4f   %4.2f %6.1f %6d %8.3f\n",
                lv, med(s$ari), mcse(s$ari), med(s$ccr), mcse(s$ccr),
                rate(s$converged), med(as.numeric(s$iter)), ok, med(s$rre_true)))
  }
  # prop_best belongs only to the restart level and is the direct answer to
  # "how much does initialisation matter": near 1 means one basin, so the 0.1.x start
  # cannot have cost much on that arm.
  pb <- raw[raw$arm == a & is.finite(raw$prop_best), ]
  if (nrow(pb))
    cat(sprintf("     random20 prop_best (share of starts reaching the best partition): %.2f\n",
                med(pb$prop_best)))
  cat("\n")
}

cat("  A difference smaller than twice the LARGER of the two MCSEs is not a\n")
cat("  difference. The MCSE is a bootstrap, valid at a ceiling and on a bimodal\n")
cat("  cell where the normal-theory form is not.\n")

if (row_id == "init")
  cat("\n  Judge every scheme against two_staged, NOT against each other. The\n  alternatives strip the random fallback, and on the heavy-tailed arms that\n  fallback is the difference between a fit and a refusal.\n")
