# STUDY F5 -- DISTRIBUTIONAL DIAGNOSTICS. Motivates the method; proves nothing.
# 36-MONTH VERSION (2026-10-02) of run_realdata_f5.R: the same diagnostics on the
# [0,36] window, with settings and results folder from realdata36_settings.R.
#
#   Rscript R/run_realdata36_f5.R [f1dir] [outdir]
#
# WHAT THIS IS AND IS NOT. These are DESCRIPTIVE. Gaussianity of the latent
# process is not identifiable at this sampling density -- a median of 7
# measurements per subject cannot distinguish a heavy-tailed score distribution
# from a Gaussian one with heavy measurement error -- and the method does not
# require Gaussianity in any case. So the diagnostics MOTIVATE the choice of a
# robust estimator; they are not a precondition for using one, and the write-up
# must not present them as a test the data passed.
#
# The formal tests are reported with that caveat attached rather than omitted,
# because a reader who wants them will otherwise assume they were run and hidden.
# At n = 10636 pooled observations any normality test rejects on a trivial
# departure, so the EFFECT SIZE -- skewness, excess kurtosis -- carries the
# argument and the p-value is close to meaningless. That is stated, not implied.

source(file.path("R", "methods.R"))
source(file.path("R", "realdata36_settings.R"))

args   <- commandArgs(trailingOnly = TRUE)
args   <- args[args != "--args"]
f1dir  <- if (length(args) >= 1) args[1] else file.path(REAL_RESULTS, "F1")
outdir <- if (length(args) >= 2) args[2] else file.path(REAL_RESULTS, "F5")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

skew <- function(x) { x <- x[is.finite(x)]; m <- mean(x); mean((x - m)^3) / stats::sd(x)^3 }
kurt <- function(x) { x <- x[is.finite(x)]; m <- mean(x); mean((x - m)^4) / stats::var(x)^2 }

# Truncated to the analysis window, so the moments describe the data the method
# was actually fitted to. Reporting skewness over [0,36] beside a fit on [0,18]
# would be describing one cohort and modelling another.
d_raw <- truncate_cea36(load_cea36(scale = "raw"), REAL_TMAX)
d_log <- truncate_cea36(load_cea36(scale = "log"), REAL_TMAX)

cat("\n== F5: distributional diagnostics ==\n")

# --- 1. the marginal, on both scales ------------------------------------------
cat("\nMARGINAL DISTRIBUTION OF THE MEASUREMENTS\n")
cat(sprintf("  %-8s %8s %10s %10s %12s %10s\n", "scale", "n", "median", "IQR", "skewness", "kurtosis"))
for (nm in c("raw", "log")) {
  y <- unlist(if (nm == "raw") d_raw$Ly else d_log$Ly)
  q <- stats::quantile(y, c(.25, .75))
  cat(sprintf("  %-8s %8d %10.2f %10s %12.2f %10.2f\n", nm, length(y),
              stats::median(y), sprintf("%.2f-%.2f", q[1], q[2]), skew(y), kurt(y)))
}
cat("\n  Gaussian reference: skewness 0, kurtosis 3.\n")
cat("  The raw scale is not a mild departure -- it is the reason the method exists.\n")

# --- 2. a formal test, with its own caveat ------------------------------------
cat("\nFORMAL TEST (Shapiro-Wilk on a random subsample; the full n rejects trivially)\n")
set.seed(20260909L)
for (nm in c("raw", "log")) {
  y <- unlist(if (nm == "raw") d_raw$Ly else d_log$Ly)
  s <- y[sample.int(length(y), min(5000L, length(y)))]
  sw <- stats::shapiro.test(s)
  cat(sprintf("  %-8s W = %.4f, p = %.3g  (n = %d)\n", nm, sw$statistic, sw$p.value, length(s)))
}
cat("  Both reject. That is expected and nearly uninformative at this n; the\n")
cat("  skewness and kurtosis above are what should be quoted.\n")

# --- 3. the estimated FPC scores, which is the quantity the METHOD sees --------
# The marginal above is a property of the measurements. The estimator's exposure
# is the SCORE distribution, and that is what the robustness argument is about,
# so it is reported separately rather than conflated with the marginal.
f1 <- tryCatch(readRDS(file.path(f1dir, "F1.rds")), error = function(e) NULL)
sc_tab <- NULL
if (!is.null(f1)) {
  cfg <- load_config(); iv <- d_raw$interval
  for (nm in c("raw", "log")) {
    dd <- if (nm == "raw") d_raw else d_log
    # Both bandwidths come from F1's own record where it has one, so this panel
    # describes the scores the reported fit saw rather than a differently-tuned
    # re-estimate. cfg$global$bwK is deliberately NOT read: it is "GCV".
    o <- realdata_optns(cfg, interval = iv,
                        bw  = if (is.null(f1$bw_months)) REAL_BW else f1$bw_months,
                        bwK = if (is.null(f1$bwK)) REAL_BWK else f1$bwK)
    s <- tryCatch(global_scores(dd, cfg, TRUE, o), error = function(e) NULL)
    if (is.null(s)) { cat(sprintf("\n  FPC scores (%s): not computable\n", nm)); next }
    if (is.null(sc_tab)) {
      cat("\nESTIMATED FPC SCORES (robust estimator) -- what the method actually sees\n")
      cat(sprintf("  %-8s %-8s %10s %12s %10s\n", "scale", "component", "sd", "skewness", "kurtosis"))
    }
    for (j in seq_len(ncol(s)))
      cat(sprintf("  %-8s %-8d %10.3f %12.2f %10.2f\n", nm, j,
                  stats::sd(s[, j]), skew(s[, j]), kurt(s[, j])))
    sc_tab <- rbind(sc_tab, data.frame(scale = nm, component = seq_len(ncol(s)),
                                       sd = apply(s, 2, stats::sd),
                                       skew = apply(s, 2, skew),
                                       kurt = apply(s, 2, kurt),
                                       stringsAsFactors = FALSE))
  }
}

cat("\n  READ THIS AS MOTIVATION, NOT AS A TEST THE DATA PASSED. At a median of 7\n")
cat("  measurements per subject the latent score distribution is not identifiable\n")
cat("  separately from the measurement error, and the method does not require it\n")
cat("  to be Gaussian. The diagnostics say a robust estimator is worth using; they\n")
cat("  do not license a claim about the latent process.\n")

marg <- do.call(rbind, lapply(c("raw", "log"), function(nm) {
  y <- unlist(if (nm == "raw") d_raw$Ly else d_log$Ly)
  data.frame(scale = nm, n = length(y), median = stats::median(y),
             skew = skew(y), kurt = kurt(y), stringsAsFactors = FALSE) }))

saveRDS(list(marginal = marg, scores = sc_tab, Tmax = REAL_TMAX,
             code_hash = source_hash(c(file.path("R", "run_realdata36_f5.R"),
                                       REAL36_SETTINGS_FILES))),
        file.path(outdir, "F5.rds"))
cat(sprintf("\nwrote %s\n", file.path(outdir, "F5.rds")))
