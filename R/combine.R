# Combine the per-task outputs of one block and summarise them.
#
#   Rscript R/combine.R <resultsdir> [outstem]
#
# Reads every *_task*.rds under resultsdir, writes the raw long table and a
# per-scenario summary, and prints the three things worth checking before any
# number is quoted:
#
#   1. COMPLETENESS. A block is not finished because the array reported success.
#      Tasks can die silently; the expected cell count is what says otherwise.
#   2. FAILURES, by method and message. A method failing is a result in itself,
#      not an inconvenience to be hidden.
#   3. NON-CONVERGENCE SPLIT BY CAUSE. `iter == maxIter` means the cap bound the
#      fit; `iter < maxIter` with converged FALSE means the small-class rule
#      fired. These need opposite remedies and pooling them destroys the
#      information.

source(file.path("config", "validate_config.R"))
source(file.path("R", "metrics.R"))
source(file.path("R", "jobs.R"))

args   <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results/core"
stem   <- if (length(args) >= 2) args[2] else file.path(resdir, "combined")
n_reps <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_

files <- list.files(resdir, pattern = "_task[0-9]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no task output found in ", resdir, call. = FALSE)
raw <- read_task_files(files)
cfg <- load_config()

cat(sprintf("\n%d task files, %d rows\n", length(files), nrow(raw)))

# --- 1. completeness ----------------------------------------------------------

cat("\n== completeness ==\n")
blk <- unique(raw$block)
if (length(blk) != 1L) stop("results mix blocks: ", paste(blk, collapse = ", "))

# Expected comes from the CONFIG, never from the rows present. Deriving it from
# the data is circular and reports 100 percent on one task out of a hundred --
# which is exactly the failure this check exists to catch, and which an earlier
# version of this file did.
expected <- expected_rows(cfg, blk, n_reps)
jobs_all <- build_jobs(cfg, blk, n_reps)
cat(sprintf("  block %s: %d cells x %d methods = %d rows expected, %d present (%.1f%%)\n",
            blk, nrow(jobs_all), expected / nrow(jobs_all), expected,
            nrow(raw), 100 * nrow(raw) / expected))
if (nrow(raw) < expected) {
  missing <- setdiff(paste(jobs_all$scenario_id, jobs_all$rep),
                     paste(raw$scenario_id, raw$rep))
  cat(sprintf("  *** INCOMPLETE: %d of %d cells absent. ***\n",
              length(missing), nrow(jobs_all)))
  cat("  *** Do not quote any number until this reads 100 percent. ***\n")
  if (length(missing) <= 12)
    cat("  missing cells:", paste(utils::head(missing, 12), collapse = ", "), "\n")
} else if (nrow(raw) > expected) {
  cat("  *** MORE rows than expected. Duplicate task output in this directory? ***\n")
}

# --- 2. failures --------------------------------------------------------------

cat("\n== failures ==\n")
bad <- raw[!is.na(raw$error), , drop = FALSE]
if (!nrow(bad)) {
  cat("  none\n")
} else {
  cat(sprintf("  %d of %d fits failed (%.2f%%)\n", nrow(bad), nrow(raw),
              100 * nrow(bad) / nrow(raw)))
  tb <- table(bad$method, substr(bad$error, 1, 44))
  for (i in seq_len(nrow(tb))) for (j in seq_len(ncol(tb)))
    if (tb[i, j] > 0)
      cat(sprintf("    %-4s %4d  %s\n", rownames(tb)[i], tb[i, j], colnames(tb)[j]))
}

# --- 3. convergence, split by cause -------------------------------------------

cat("\n== convergence, by cause ==\n")
ok <- raw[is.na(raw$error) & !is.na(raw$converged), , drop = FALSE]
if (nrow(ok)) {
  for (m in sort(unique(ok$method))) {
    s <- ok[ok$method == m, ]
    if (all(is.na(s$stop_reason))) next
    tb <- table(factor(s$stop_reason,
                       levels = c("converged", "cycle", "maxIter",
                                  "small_cluster", "trivial")))
    cat(sprintf("  %-4s  %s   median iter %s\n", m,
                paste(sprintf("%s %d", names(tb), as.integer(tb)), collapse = "  "),
                ifelse(all(is.na(s$iter)), "-", format(stats::median(s$iter, na.rm = TRUE)))))
  }
} else cat("  no convergence information (score-only methods)\n")

# --- summary ------------------------------------------------------------------
#
# THE POOLED MEDIAN IS THE WRONG STATISTIC ON A BIMODAL CELL, and at least one
# cell here is expected to be bimodal. On the infinite-variance arm (arm 4) the
# median can move from about 0.5 to near 0 between the two halves of ONE
# 100-replicate run, same method, same configuration, because the outcome is
# bimodal with the failure fraction near one half, so the median reports which
# side the sample fell on rather than typical performance.
#
# So every cell also reports P(aRI > success_thr) and the median GIVEN success.
# "A method produces the best partitions of any method when it works and works
# least often" is a coherent finding that is completely invisible in the pooled
# median. This is NOT conditioning on convergence, which would be selection
# bias: the success RATE is itself reported, so nothing is hidden.
SUCCESS_THR <- 0.05

num <- c("ccr", "ari", "nmi", "rre_in", "rre_true", "rre_out", "seconds")
agg <- function(df, v) {
  x <- df[[v]]; x <- x[is.finite(x)]
  if (!length(x)) return(c(NA, NA, NA, NA))
  c(stats::median(x), stats::quantile(x, .25, names = FALSE),
    stats::quantile(x, .75, names = FALSE), mcse_median(x))
}
key <- unique(raw[, c("block", "scenario_id", "arm", "d", "method")])
key <- key[!is.na(key$method), , drop = FALSE]
sm <- do.call(rbind, lapply(seq_len(nrow(key)), function(i) {
  s <- raw[raw$scenario_id == key$scenario_id[i] &
             !is.na(raw$method) & raw$method == key$method[i], ]
  row <- key[i, ]
  # SURVIVORSHIP. Every median below is taken over the replicates on which that
  # method RETURNED a partition. n_ok is what makes a cell built from 89 fits
  # distinguishable from one built from 100; without it a comparator can hold a
  # win on the replicates it happened to be able to fit.
  row$n_fits <- nrow(s); row$n_failed <- sum(!is.na(s$error))
  row$n_ok   <- sum(is.na(s$error))
  row$conv_rate <- mean(s$converged, na.rm = TRUE)

  a <- s$ari[is.finite(s$ari)]
  row$p_success <- if (length(a)) mean(a > SUCCESS_THR) else NA_real_
  row$ari_med_given_success <-
    if (sum(a > SUCCESS_THR)) stats::median(a[a > SUCCESS_THR]) else NA_real_
  # A crude bimodality flag, not a test: the share of replicates landing in the
  # middle third of the observed range. A low value with a mid-range median is
  # the signature that makes the pooled median unstable.
  row$mid_share <- if (length(a) >= 10L && diff(range(a)) > 1e-8) {
    rg <- range(a); mean(a > rg[1] + diff(rg)/3 & a < rg[2] - diff(rg)/3)
  } else NA_real_

  for (v in num) {
    a2 <- agg(s, v)
    row[[paste0(v, "_med")]]  <- a2[1]; row[[paste0(v, "_q25")]] <- a2[2]
    row[[paste0(v, "_q75")]]  <- a2[3]; row[[paste0(v, "_mcse")]] <- a2[4]
  }
  row
}))

saveRDS(raw, paste0(stem, "_raw.rds"))
utils::write.csv(sm, paste0(stem, "_summary.csv"), row.names = FALSE)
cat(sprintf("\nwrote %s_raw.rds and %s_summary.csv\n", stem, stem))

# --- headline table -----------------------------------------------------------

cat("\n== CCR and aRI by scenario and method ==\n")
cat("   n     replicates the method RETURNED a partition on. A median over 89\n")
cat("         fits is not comparable to one over 100 -- read this column first.\n")
cat(sprintf("   P(s)  share of replicates with aRI > %.2f, and the median GIVEN that.\n",
            SUCCESS_THR))
cat("         On a bimodal cell the pooled median reports which side the sample\n")
cat("         fell on, so P(s) and the conditional median carry the conclusion.\n")
cat("   MCSE  bootstrap SE of the median, in the summary csv. A difference\n")
cat("         smaller than twice the LARGER of the two MCSEs is not a difference.\n")
cat("         Bootstrap, not the normal-theory form, which was measured 84 percent\n")
cat("         too SMALL on a bimodal cell and 340 percent too large at a ceiling.\n\n")
for (sc in unique(sm$scenario_id)) {
  s <- sm[sm$scenario_id == sc, ]
  lab <- cfg$arms$label[cfg$arms$arm == s$arm[1]]
  cat(sprintf("  %-10s arm %-22s d = %-4g\n", sc, lab, s$d[1]))
  cat(sprintf("     %-4s  %-22s %-22s %5s  %5s %6s  %s\n",
              "", "CCR med [IQR]", "aRI med [IQR]", "n", "P(s)", "aRI|s", "conv"))
  for (i in seq_len(nrow(s))) {
    flag <- if (!is.na(s$mid_share[i]) && s$mid_share[i] < 0.20 &&
                s$ari_med[i] > 0.05 && s$ari_med[i] < 0.95) "  <- BIMODAL?" else ""
    cat(sprintf("     %-4s  %.3f [%.3f, %.3f]      %+.3f [%+.3f, %+.3f]   %3d  %5.2f %6s  %.2f%s%s\n",
                s$method[i], s$ccr_med[i], s$ccr_q25[i], s$ccr_q75[i],
                s$ari_med[i], s$ari_q25[i], s$ari_q75[i],
                s$n_ok[i], s$p_success[i],
                ifelse(is.na(s$ari_med_given_success[i]), "-",
                       sprintf("%.3f", s$ari_med_given_success[i])),
                s$conv_rate[i],
                if (s$n_failed[i] > 0) sprintf("  (%d failed)", s$n_failed[i]) else "",
                flag))
  }
}
cat("\n   A cell flagged BIMODAL has under 20 percent of its replicates in the\n")
cat("   middle third of its own range while the median sits mid-scale. That is a\n")
cat("   HEURISTIC, not a test: confirm it by splitting the replicates in half and\n")
cat("   checking whether the median moves. Report the rate and the conditional\n")
cat("   median for such a cell, and cite the split as the justification.\n\n")
