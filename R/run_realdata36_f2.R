# STUDY F2 -- LANDMARK SURVIVAL.
# 36-MONTH VERSION (2026-10-02) of run_realdata_f2.R: the same analysis, with the
# window, landmark and results folder from realdata36_settings.R.
#
#   Rscript R/run_realdata36_f2.R [f1dir] [outdir]
#
# THE DEFECT THIS FIXES, and it is the most serious comment in either report.
# Group membership is built from CEA over months 0-36, then used as a BASELINE
# exposure for survival measured from surgery. That reverses the temporal
# ordering: the exposure is not known at time zero. Worse, it confers GUARANTEE
# TIME -- entering the cohort at all required surviving long enough to accrue
# four measurements, so a patient who died at month 8 could never be assigned to
# a "slowly rising" group. Any survival difference is then partly an artefact of
# who lived long enough to be classified.
#
# THE FIX: a LANDMARK at the END OF THE TRAJECTORY WINDOW, MONTH 36 here (the
# clustering window is [0,36], as in the submitted analysis). REAL_LANDMARK ties
# the two together in realdata36_settings.R.
#
#   * time origin moves to month 36
#   * eligible set is patients ALIVE AND UNDER FOLLOW-UP at month 36; patients
#     who died or were censored during the window still contribute their
#     measurements to the clustering, and are COUNTED here rather than silently
#     dropped
#   * the estimand becomes survival from month 36 onward among those alive at
#     month 36, which is what must be written in the Abstract, Results and
#     Discussion in identical words.
#
# THE COST OF THE LONGER WINDOW: 1259 of 1405 subjects reach a month-36
# landmark (1341 of 1395 at month 18), and 127 deaths follow it (187 after month
# 18). The exclusion is 10.4% (3.9% at month 18).
#
# Survival_Time is in MONTHS (verified: median 70.2, max 143.3), so the landmark
# subtraction is direct and no unit conversion is involved.

source(file.path("R", "methods.R"))
source(file.path("R", "realdata36_settings.R"))

args   <- commandArgs(trailingOnly = TRUE)
args   <- args[args != "--args"]
f1dir  <- if (length(args) >= 1) args[1] else file.path(REAL_RESULTS, "F1")
outdir <- if (length(args) >= 2) args[2] else file.path(REAL_RESULTS, "F2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("survival", quietly = TRUE))
  stop("package `survival` is required for F2.", call. = FALSE)

# THE LANDMARK SITS AT THE END OF THE 36-MONTH WINDOW -- the grouping is not
# known before then. Both come from REAL_LANDMARK in realdata36_settings.R, which
# ties them so one cannot be changed without the other.
LANDMARK <- REAL_LANDMARK

f1 <- readRDS(file.path(f1dir, "F1.rds"))
d  <- truncate_cea36(load_cea36(), LANDMARK)
b  <- d$base
pk <- peak_cea(d)

# THE PARTITION F2 USES. F1 fits four ways; survival is reported against the
# PRIMARY analysis -- the robust method on the scale F1 designates -- and the
# others are a sensitivity. Taking whichever partition happens to be first would
# make the survival result depend on an ordering nobody chose.
# From realdata_load.R, so F2 and F3 cannot disagree; the fallback keeps this
# runnable against a loader uploaded before the constant existed.
prim <- if (exists("REAL_PRIMARY")) REAL_PRIMARY else "raw_M1"
if (is.null(f1$fits[[prim]]))
  stop("F1 has no `", prim, "` partition; run F1 first.", call. = FALSE)

results <- list()

for (key in names(f1$fits)) {
  g <- f1$fits[[key]]
  if (length(g) != nrow(b))
    stop(sprintf("partition `%s` has %d labels for %d subjects.", key, length(g), nrow(b)),
         call. = FALSE)

  # --- the landmark, and the accounting it demands ---------------------------
  at_risk <- b$Survival_Time >= LANDMARK
  n_all   <- nrow(b)
  n_lm    <- sum(at_risk)

  # Time and status measured FROM the landmark, for those who reached it.
  tt <- b$Survival_Time[at_risk] - LANDMARK
  ev <- b$Death[at_risk]
  gg <- factor(g[at_risk])

  sd_fit <- survival::survdiff(survival::Surv(tt, ev) ~ gg)
  p_lr   <- stats::pchisq(sd_fit$chisq, df = length(sd_fit$n) - 1, lower.tail = FALSE)
  km     <- survival::survfit(survival::Surv(tt, ev) ~ gg)
  cx     <- survival::coxph(survival::Surv(tt, ev) ~ gg)
  sm     <- summary(cx)

  # THE NAIVE ANALYSIS, computed alongside so that the size of the bias can be
  # stated rather than only its direction. This is survival from SURGERY with
  # the same grouping -- exactly what the submitted paper did.
  cx_naive <- survival::coxph(survival::Surv(b$Survival_Time, b$Death) ~ factor(g))
  sm_naive <- summary(cx_naive)

  results[[key]] <- list(
    partition = key,
    n_source = n_all, n_landmark = n_lm, n_excluded = n_all - n_lm,
    excluded_pct = 100 * (n_all - n_lm) / n_all,
    group_sizes_all = as.integer(table(factor(g))),
    group_sizes_landmark = as.integer(table(gg)),
    events_landmark = as.integer(tapply(ev, gg, sum)),
    followup_median = stats::median(tt),
    logrank_chisq = unname(sd_fit$chisq), logrank_p = p_lr,
    hr = unname(sm$conf.int[1, 1]),
    hr_lo = unname(sm$conf.int[1, 3]), hr_hi = unname(sm$conf.int[1, 4]),
    hr_p = unname(sm$coefficients[1, 5]),
    hr_naive = unname(sm_naive$conf.int[1, 1]),
    hr_naive_lo = unname(sm_naive$conf.int[1, 3]),
    hr_naive_hi = unname(sm_naive$conf.int[1, 4]),
    hr_naive_p = unname(sm_naive$coefficients[1, 5]),
    km = summary(km, times = c(0, 12, 24, 36, 48, 60), extend = TRUE)[
           c("time", "strata", "surv", "lower", "upper", "n.risk")])
}

# --- report -------------------------------------------------------------------
cat(sprintf("\n== F2: landmark survival at %d months ==\n", LANDMARK))
r <- results[[prim]]
cat(sprintf("\nEXCLUSION ACCOUNTING (primary analysis, %s)\n", prim))
cat(sprintf("  %d subjects clustered\n", r$n_source))
cat(sprintf("  %d died or were censored before month %d and are NOT at risk (%.1f%%)\n",
            r$n_excluded, LANDMARK, r$excluded_pct))
cat(sprintf("  %d analysed from the landmark, median follow-up %.1f months after it\n",
            r$n_landmark, r$followup_median))
cat(sprintf("  group sizes: all %s | at landmark %s | events %s\n",
            paste(r$group_sizes_all, collapse = "/"),
            paste(r$group_sizes_landmark, collapse = "/"),
            paste(r$events_landmark, collapse = "/")))

cat("\nHAZARD RATIO, group 2 against group 1\n")
cat(sprintf("  LANDMARK (correct) : HR %.2f (%.2f-%.2f), p = %.4g\n",
            r$hr, r$hr_lo, r$hr_hi, r$hr_p))
cat(sprintf("  from surgery (the submitted design, biased): HR %.2f (%.2f-%.2f), p = %.4g\n",
            r$hr_naive, r$hr_naive_lo, r$hr_naive_hi, r$hr_naive_p))
cat(sprintf("  log-rank at the landmark: chisq %.2f, p = %.4g\n", r$logrank_chisq, r$logrank_p))
cat("\n  The second line is reported so that the SIZE of the guarantee-time\n")
cat("  bias can be stated, not only its direction.\n")

# --- DOES THE CLUSTERING BEAT ONE NUMBER? -------------------------------------
# F0C found that peak CEA reproduces the partition at AUC 0.961. That does not
# make the clusters fake -- their mean curves differ in shape, one rebounding to
# a peak at month 13 while the other stays flat -- but it does raise the only
# question that matters for a clinical claim: does cluster membership predict
# survival BETTER THAN A THRESHOLD ON PEAK CEA? If it does not, the honest
# reading is that the trajectory machinery recovered something a single lab value
# already carried, and the paper should say so.
#
# The threshold is the one that best separates the CLUSTERS (not survival), so
# this is the fairest available stand-in for "what a clinician would have done
# instead", and it is not tuned on the outcome.
cat("\n== THE COMPARATOR F0C MADE NECESSARY: cluster against peak CEA ==\n")
# Rebuilt from `prim` explicitly. `tt`, `ev` and `gg` above are LOOP variables
# and hold whichever partition the loop happened to end on, which is not the
# primary one -- reusing them here would have compared the wrong grouping.
gp      <- f1$fits[[prim]]
lab_min <- gp == which.min(table(factor(gp)))
ar      <- b$Survival_Time >= LANDMARK
tt      <- b$Survival_Time[ar] - LANDMARK
ev      <- b$Death[ar]
gg      <- factor(gp[ar])
o <- order(pk); cum <- c(0, cumsum(lab_min[o])); tot <- cum[length(pk) + 1L]
idx <- c(0L, which(pk[o][-length(pk)] != pk[o][-1]), length(pk))
accs <- ((idx - cum[idx + 1L]) + (tot - cum[idx + 1L])) / length(pk)
cut_at <- pk[o][max(1L, idx[which.max(pmax(accs, 1 - accs))])]
pk_hi  <- pk > cut_at

cx_pk <- survival::coxph(survival::Surv(tt, ev) ~ factor(pk_hi[ar]))
cx_cl <- survival::coxph(survival::Surv(tt, ev) ~ gg)
cx_bo <- survival::coxph(survival::Surv(tt, ev) ~ gg + factor(pk_hi[ar]))
cmp <- data.frame(
  model = c("cluster only", "peak CEA only", "both"),
  hr_first = c(summary(cx_cl)$conf.int[1, 1], summary(cx_pk)$conf.int[1, 1],
               summary(cx_bo)$conf.int[1, 1]),
  p_first = c(summary(cx_cl)$coefficients[1, 5], summary(cx_pk)$coefficients[1, 5],
              summary(cx_bo)$coefficients[1, 5]),
  concordance = c(summary(cx_cl)$concordance[1], summary(cx_pk)$concordance[1],
                  summary(cx_bo)$concordance[1]),
  loglik = c(stats::logLik(cx_cl), stats::logLik(cx_pk), stats::logLik(cx_bo)),
  stringsAsFactors = FALSE)
cat(sprintf("  peak-CEA threshold that best reproduces the clusters: %.2f\n", cut_at))
cat(sprintf("  %-16s %12s %10s %14s %12s\n", "model", "HR (1st term)", "p",
            "concordance", "logLik"))
for (i in seq_len(nrow(cmp)))
  cat(sprintf("  %-16s %12.2f %10.4g %14.4f %12.1f\n", cmp$model[i], cmp$hr_first[i],
              cmp$p_first[i], cmp$concordance[i], cmp$loglik[i]))
lr <- 2 * (as.numeric(stats::logLik(cx_bo)) - as.numeric(stats::logLik(cx_pk)))
cat(sprintf("\n  Does the CLUSTER add anything once peak CEA is in the model?\n"))
cat(sprintf("    likelihood-ratio chisq %.2f on 1 df, p = %.4g\n", lr,
            stats::pchisq(lr, df = 1, lower.tail = FALSE)))
cat("  A non-significant value here means the partition carries no survival\n")
cat("  information beyond a single threshold on peak CEA. That is reportable and\n")
cat("  it is not a disaster -- but it must not be written up as a trajectory\n")
cat("  subtype effect, and the Discussion has to say which way this came out.\n")

cat("\nSENSITIVITY: the other partitions\n")
for (k in names(results)) if (k != prim) {
  s <- results[[k]]
  cat(sprintf("  %-8s HR %.2f (%.2f-%.2f) p=%.3g | at landmark %d, groups %s\n",
              k, s$hr, s$hr_lo, s$hr_hi, s$hr_p, s$n_landmark,
              paste(s$group_sizes_landmark, collapse = "/")))
}

saveRDS(list(landmark = LANDMARK, primary = prim, results = results,
             peak_cut = cut_at, comparator = cmp, lr_chisq = lr,
             code_hash = source_hash(c(file.path("R", "run_realdata36_f2.R"),
                                       REAL36_SETTINGS_FILES))),
        file.path(outdir, "F2.rds"))
cat(sprintf("\nwrote %s\n", file.path(outdir, "F2.rds")))
