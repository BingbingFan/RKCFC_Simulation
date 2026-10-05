# Combine one decisions factor and print the table the decision rests on.
#
#   Rscript R/combine_decisions.R <resultsdir> [n_reps]
#
# One table per factor, laid out so the choice is readable rather than derivable.
# Medians with the bootstrap Monte Carlo SE beside them, convergence and
# iteration counts where they bear on whether a number is reportable at all, and
# the failure count so a median over 17 replicates is never mistaken for one over
# 20.

source(file.path("config", "validate_config.R"))
source(file.path("R", "metrics.R"))

args   <- commandArgs(trailingOnly = TRUE)
args   <- args[args != "--args"]
resdir <- if (length(args) >= 1) args[1] else "results/D_nkfit"
n_reps <- if (length(args) >= 2) as.integer(args[2]) else 20L

files <- list.files(resdir, pattern = "_task[0-9]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no task output in ", resdir, call. = FALSE)
raw <- read_task_files(files)
cfg <- load_config()
fac <- unique(raw$factor)
if (length(fac) != 1L) stop("results mix factors: ", paste(fac, collapse = ", "))

arms <- sort(unique(raw$arm))
cat(sprintf("\n%d task files, %d rows | factor %s | arms %s | %d replicates\n",
            length(files), nrow(raw), fac, paste(arms, collapse = ", "),
            length(unique(raw$rep))))

# --- completeness, from the config rather than from the rows present ----------
n_lev  <- length(unique(raw$level))
n_meth <- length(unique(raw$method))
expected <- length(arms) * n_reps * n_lev * n_meth
cat(sprintf("  %d rows expected (%d arms x %d reps x %d levels x %d methods), %d present (%.1f%%)\n",
            expected, length(arms), n_reps, n_lev, n_meth, nrow(raw),
            100 * nrow(raw) / expected))
if (nrow(raw) < expected)
  cat("  *** INCOMPLETE. Do not quote any number until this reads 100 percent. ***\n")
# EXCESS IS AS BAD AS SHORTFALL, and it is the quieter of the two. Every combiner
# globs _task*.rds and rbinds whatever it finds, so STALE FILES from an earlier
# run are merged in silently and the tables below mix two designs. Caught on
# 2026-09-06 in a smoke directory that still held task files from the day before;
# the only symptom was a completeness figure above 100 percent, which read as
# reassuring rather than alarming.
if (nrow(raw) > expected) {
  cat("  *** MORE ROWS THAN THE DESIGN HAS -- almost always STALE TASK FILES.\n")
  cat("  *** Check the file dates in this directory, clear it, and re-run.\n")
  cat("  *** Do not quote any number until this reads exactly 100 percent. ***\n")
}

nfail <- sum(!is.na(raw$error))
cat(sprintf("  %d of %d fits failed (%.1f%%)\n", nfail, nrow(raw),
            100 * nfail / nrow(raw)))
if (nfail) {
  tb <- table(raw$level[!is.na(raw$error)], substr(raw$error[!is.na(raw$error)], 1, 46))
  for (i in seq_len(nrow(tb))) for (j in seq_len(ncol(tb)))
    if (tb[i, j] > 0) cat(sprintf("    %-16s %3d  %s\n",
                                  rownames(tb)[i], tb[i, j], colnames(tb)[j]))
}

med  <- function(x) { x <- x[is.finite(x)]; if (length(x)) stats::median(x) else NA_real_ }
mcse <- function(x) { x <- x[is.finite(x)]; if (length(x) >= 5L) mcse_median(x) else NA_real_ }
rate <- function(x) { x <- x[!is.na(x)]; if (length(x)) mean(x) else NA_real_ }

# --- the delta factor needs no fitting and gets its own layout ----------------

if (fac == "delta") {
  cat("\n== Delta by arm and d, median over replicates ==\n")
  cat("   Near-constant DOWN each column is what the IQR matching promises.\n\n")
  lev <- unique(raw$level)
  lev <- lev[order(as.numeric(sub("d=", "", lev)))]
  cat(sprintf("   %-24s %s\n", "arm", paste(sprintf("%9s", lev), collapse = "")))
  tab <- matrix(NA_real_, length(arms), length(lev))
  for (a in seq_along(arms)) {
    for (j in seq_along(lev)) {
      s <- raw[raw$arm == arms[a] & raw$level == lev[j], ]
      tab[a, j] <- med(s$delta)
    }
    # Look the label up BY ARM NUMBER, not by row position. They coincide today
    # because arms.csv happens to be in order, and would stop coinciding the
    # moment a row is reordered or an arm is dropped from a run.
    cat(sprintf("   %-24s %s\n",
                substr(cfg$arms$label[cfg$arms$arm == arms[a]], 1, 22),
                paste(sprintf("%9s", sprintf("%.2f", tab[a, ])), collapse = "")))
  }
  sp <- apply(tab, 2, function(z) max(z, na.rm = TRUE) / min(z, na.rm = TRUE))
  cat(sprintf("\n   %-24s %s\n", "cross-arm spread",
              paste(sprintf("%9s", sprintf("%.2fx", sp)), collapse = "")))
  cat(sprintf("   %-24s %s\n", "median reference CCR",
              paste(sprintf("%9s", sprintf("%.3f", apply(tab, 2, function(z)
                stats::pnorm(stats::median(z, na.rm = TRUE) / 2)))), collapse = "")))
  cat("\n   A spread above about 1.3x means the arms are NOT exchangeable at that\n")
  cat("   d, so accuracy must be reported RELATIVE to each arm's own reference\n")
  cat("   rather than compared across arms directly.\n\n")
  quit(status = 0L)
}

# --- the fit-based factors ----------------------------------------------------

lab <- switch(fac,
  nkfit = "RETAINED COMPONENTS. Read `conv` and `iter` FIRST: a fit pinned on\n   maxIter is an artefact of where the loop was cut off, not a result.",
  d = "SEPARATION. A level where the methods return the SAME number cannot rank\n   them, however high that number is.",
  convprop = "CONVERGENCE THRESHOLD. `agree` is the statistic this exists for:\n   near 1 means the partition settled and only boundary subjects oscillate,\n   so the non-convergence is COSMETIC. Do NOT relax cp until convergence\n   looks good and then report that rate.",
  init = "INITIALISATION. `oracle` bounds what any initialiser could buy;\n   `pace_kmeans` reproduces the 0.1.x start deliberately.")
cat(sprintf("\n== %s ==\n\n", lab))

for (a in arms) {
  cat(sprintf("  arm %d  %s\n", a, cfg$arms$label[cfg$arms$arm == a]))
  cat(sprintf("     %-16s %-6s %-18s %-18s %6s %6s %6s %8s\n",
              "level", "method", "aRI med +/- mcse", "CCR med +/- mcse",
              "conv", "iter", "n_ok", if (fac == "convprop") "agree" else "rre_true"))
  s0 <- raw[raw$arm == a, ]
  lv <- unique(s0$level)
  if (fac %in% c("nkfit", "d"))
    lv <- lv[order(as.numeric(gsub("[^0-9.]", "", lv)))]
  for (l in lv) for (m in sort(unique(s0$method))) {
    s <- s0[s0$level == l & s0$method == m, ]
    if (!nrow(s)) next
    extra <- if (fac == "convprop") med(s$agree_strict) else med(s$rre_true)
    cat(sprintf("     %-16s %-6s %7.3f +/- %-6.4f %7.3f +/- %-6.4f %6.2f %6s %6d %8s\n",
                l, m, med(s$ari), mcse(s$ari), med(s$ccr), mcse(s$ccr),
                rate(s$converged),
                ifelse(is.na(med(s$iter)), "-", format(med(s$iter))),
                sum(is.na(s$error)),
                ifelse(is.na(extra), "-", sprintf("%.3f", extra))))
  }
  if (fac == "init" && any(!is.na(s0$prop_best)))
    cat(sprintf("     random20 prop_best (share of 20 starts reaching the best partition): %.2f\n",
                med(s0$prop_best)))
  cat("\n")
}

cat("  A difference smaller than twice the LARGER of the two MCSEs is not a\n")
cat("  difference. The MCSE is a bootstrap, valid at a ceiling and on a bimodal\n")
cat("  cell where the normal-theory form is not.\n\n")
