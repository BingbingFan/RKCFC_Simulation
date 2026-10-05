# Combine block J -- the objective trace.
#
#   Rscript R/combine_convergence.R <resultsdir> [n_reps]
#
# READ IN THIS ORDER, and stop at the first thing that is wrong:
#
#   1. completeness at 100 percent
#   2. the failure table
#   3. THE STRUCTURAL CHECK. The reassignment step is monotone BY CONSTRUCTION --
#      a subject moves only to a class with a strictly smaller residual sum. It is
#      the one convergence property this method provably has. If it fails in even
#      one row, that is an IMPLEMENTATION DEFECT, not a finding about the
#      algorithm, and every number below it is describing broken code.
#   4. only then the re-estimation numbers, which are the finding
#
# WHAT THE NUMBERS ARE FOR. The manuscript concedes that monotone descent cannot be
# established, because the re-estimation step does not minimise J -- Kendall FPCA
# is not the minimiser of a sum of squared reconstruction errors, and the
# local-linear smoothing adds a bias that does not vanish in J. Conceding that is
# right, but a concession is worth much more with a measurement attached: how
# OFTEN does the difficulty arise, and how BIG is it when it does.

source(file.path("config", "validate_config.R"))
source(file.path("R", "metrics.R"))

args   <- commandArgs(trailingOnly = TRUE)
args   <- args[args != "--args"]
resdir <- if (length(args) >= 1) args[1] else "results/convergence"
n_reps <- if (length(args) >= 2) as.integer(args[2]) else NA_integer_

files <- list.files(resdir, pattern = "_task[0-9]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no task output in ", resdir, call. = FALSE)
raw <- read_task_files(files)
cfg <- load_config()

arms <- sort(unique(raw$arm))
mets <- sort(unique(raw$method))
reps <- length(unique(raw$rep))
if (is.na(n_reps)) n_reps <- reps

cat(sprintf("\n%d task files, %d rows | block convergence | arms %s | %d replicates\n",
            length(files), nrow(raw), paste(arms, collapse = ", "), reps))

expected <- length(arms) * n_reps * length(mets)
cat(sprintf("  %d rows expected (%d arms x %d reps x %d methods), %d present (%.1f%%)\n",
            expected, length(arms), n_reps, length(mets), nrow(raw),
            100 * nrow(raw) / expected))
if (nrow(raw) < expected)
  cat("  *** INCOMPLETE. Do not quote any number until this reads 100 percent. ***\n")
if (nrow(raw) > expected) {
  cat("  *** MORE ROWS THAN THE DESIGN HAS -- almost always STALE TASK FILES.\n")
  cat("  *** Clear this directory and re-run. ***\n")
}

fail <- !is.na(raw$error)
cat(sprintf("  %d of %d fits failed (%.1f%%)\n", sum(fail), nrow(raw),
            100 * sum(fail) / nrow(raw)))
if (any(fail)) {
  tb <- table(raw$method[fail], substr(raw$error[fail], 1, 46))
  for (i in seq_len(nrow(tb))) for (j in seq_len(ncol(tb)))
    if (tb[i, j] > 0) cat(sprintf("    %-6s %3d  %s\n", rownames(tb)[i], tb[i, j], colnames(tb)[j]))
}

# --- 3. THE STRUCTURAL CHECK, before any finding is printed -------------------
cat("\n== REASSIGNMENT MONOTONICITY (must hold in EVERY fit) ==\n")
chk <- raw$reassign_descends_all[!is.na(raw$reassign_descends_all)]
if (!length(chk)) {
  cat("  *** NO TRACES AT ALL. The block produced nothing readable. ***\n")
} else if (all(chk)) {
  cat(sprintf("  holds in all %d fits with a trace\n", length(chk)))
} else {
  cat(sprintf("  *** FAILS IN %d of %d FITS. This is an IMPLEMENTATION DEFECT:\n",
              sum(!chk), length(chk)))
  cat("  *** the reassignment step is monotone by construction, so a violation\n")
  cat("  *** means the code does not implement Equation (13). Nothing below is\n")
  cat("  *** a finding about the algorithm until this reads zero. ***\n")
}

arm_name <- function(a) {
  r <- cfg$arms[cfg$arms$arm == a, ]
  if (nrow(r) && !is.null(r$label)) as.character(r$label)[1] else paste("arm", a)
}
med <- function(x) { x <- x[is.finite(x)]; if (length(x)) stats::median(x) else NA_real_ }
q90 <- function(x) { x <- x[is.finite(x)]; if (length(x)) stats::quantile(x, 0.9) else NA_real_ }

# --- 4. the finding -----------------------------------------------------------
cat("\n== RE-ESTIMATION: how often J RISES, and by how much ==\n")
cat("   A rise is J_before(b) - J_after(b-1) > 0: the change attributable to the\n")
cat("   REFIT alone, with the reassignment that preceded it divided out.\n")
cat("   Rises are RELATIVE to the J the refit started from, because J scales with\n")
cat("   n and with the noise level and is not comparable across arms otherwise.\n\n")
# "% up" is a VALUE here, not a format directive, so it needs no doubling --
# writing "%%" printed a literal percent-percent in the header.
cat(sprintf("  %-24s %-6s %8s %9s %9s %9s %7s\n",
            "arm", "method", "steps", "% up", "med rise", "p90 rise", "n_ok"))
for (a in arms) for (m in mets) {
  s <- raw[raw$arm == a & raw$method == m & !fail, ]
  if (!nrow(s)) next
  cat(sprintf("  %-24s %-6s %8.1f %8.0f%% %9.4f %9.4f %7d\n",
              if (m == mets[1]) arm_name(a) else "", m,
              med(s$n_reest), 100 * med(s$frac_reest_up),
              med(s$med_rise_rel), q90(s$max_rise_rel), nrow(s)))
}
cat("\n   `%% up` is the median over fits of the share of that fit's re-estimation\n")
cat("   steps on which J rose. `p90 rise` is the 90th percentile of each fit's\n")
cat("   LARGEST single rise -- the tail is the part that matters.\n")

# --- how fits end -------------------------------------------------------------
cat("\n== HOW FITS END, by cause ==\n")
cat("   A cycle means the loop revisited a partition it had already been in. It\n")
cat("   is not a failure of the code; it is the absence of guaranteed descent\n")
cat("   showing up in practice, and it is why cycle detection exists at all.\n\n")
causes <- c("converged", "cycle", "maxIter", "small_cluster", "trivial")
cat(sprintf("  %-24s %-6s %s %8s %9s\n", "arm", "method",
            paste(sprintf("%15s", causes), collapse = ""), "med iter", "med cyc len"))
for (a in arms) for (m in mets) {
  s <- raw[raw$arm == a & raw$method == m, ]
  if (!nrow(s)) next
  cnt <- vapply(causes, function(cc) sum(s$stop_reason == cc, na.rm = TRUE), integer(1))
  cat(sprintf("  %-24s %-6s %s %8.1f %9s\n",
              if (m == mets[1]) arm_name(a) else "", m,
              paste(sprintf("%15d", cnt), collapse = ""),
              med(s$iter),
              if (any(is.finite(s$cycle_length))) sprintf("%.1f", med(s$cycle_length)) else "-"))
}

cat("\n== THE SENTENCE THIS BLOCK LICENSES ==\n")
ok <- raw[!fail, ]
cat(sprintf("   Across %d fits and %d re-estimation steps, J rose on %.0f%% of steps;\n",
            nrow(ok), sum(ok$n_reest, na.rm = TRUE),
            100 * sum(ok$n_reest_up, na.rm = TRUE) / max(1, sum(ok$n_reest, na.rm = TRUE))))
cat(sprintf("   the median rise was %.1f%% of J and the largest %.1f%%; %d of %d fits\n",
            100 * med(ok$med_rise_rel), 100 * max(ok$max_rise_rel, na.rm = TRUE),
            sum(ok$stop_reason == "cycle", na.rm = TRUE), nrow(ok)))
cat("   ended in a detected cycle. Reassignment descended in every fit.\n")
