# READ STUDY F5 -- distributional diagnostics. Motivates the method; proves nothing.
#
#   Rscript R/read_realdata_f5.R [f5dir] [f1dir]
#
# WHAT THIS IS. Descriptive evidence that a robust estimator is worth using on
# this cohort. It is NOT a test the data passed, and it is not a premise the
# method requires: at a median of 7 measurements per subject the latent score
# distribution is not identifiable separately from the measurement error, so
# nothing here licenses a claim about the latent process. Write "motivates";
# never write "confirms".
#
# THE OBJECTION THIS FILE MUST MEET, and it is the strongest one available to a
# reader: if the log transform makes the data approximately Gaussian, then the
# honest recommendation is to take logs and use the ordinary estimator, and the
# robust machinery is unnecessary. Panel 2 puts the two scales side by side so
# that comparison is made here rather than left for someone else to make. If the
# log scale IS well behaved, that must be conceded in the Discussion and the
# argument moved to where it still holds -- F1's raw-scale panel, and the
# simulation arms where contamination is not a transform away.

args  <- commandArgs(trailingOnly = TRUE)
args  <- args[args != "--args"]
f5dir <- if (length(args) >= 1) args[1] else file.path("results_real", "F5")
f1dir <- if (length(args) >= 2) args[2] else file.path("results_real", "F1")

fp <- file.path(f5dir, "F5.rds")
if (!file.exists(fp)) stop("not found: ", fp, "\n  Has F5 finished?", call. = FALSE)
f5 <- readRDS(fp)
mg <- f5$marginal
sc <- f5$scores

fmt <- function(x, k = 2) if (is.finite(x)) sprintf(paste0("%.", k, "f"), x) else "-"
get <- function(tabl, s, col) { v <- tabl[[col]][tabl$scale == s]
                                if (length(v)) v[1] else NA_real_ }

cat("\n== F5: distributional diagnostics (descriptive) ==\n")
cat(sprintf("  code_hash %s\n", if (is.null(f5$code_hash)) "-" else f5$code_hash))

# --- 1. the marginal ----------------------------------------------------------
cat("\nPANEL 1. THE MEASUREMENTS\n")
cat(sprintf("  %-8s %10s %12s %12s %14s %12s\n", "scale", "n", "median",
            "skewness", "excess kurt", "kurtosis"))
for (i in seq_len(nrow(mg)))
  cat(sprintf("  %-8s %10d %12s %12s %14s %12s\n", mg$scale[i], mg$n[i],
              fmt(mg$median[i]), fmt(mg$skew[i]), fmt(mg$kurt[i] - 3),
              fmt(mg$kurt[i])))
cat("\n  Gaussian reference: skewness 0, excess kurtosis 0 (kurtosis 3).\n")
cat("  EXCESS kurtosis is the column to quote -- the raw fourth moment is often\n")
cat("  reported against the wrong reference, and 3 versus 0 is exactly the kind\n")
cat("  of slip a reader reads as carelessness about the rest.\n")

# --- 2. does the transform already fix it -------------------------------------
cat("\nPANEL 2. DOES THE LOG TRANSFORM ALREADY SOLVE THIS\n")
sr <- get(mg, "raw", "skew"); sl <- get(mg, "log", "skew")
kr <- get(mg, "raw", "kurt"); kl <- get(mg, "log", "kurt")
cat(sprintf("  skewness        %8s  ->  %8s\n", fmt(sr), fmt(sl)))
cat(sprintf("  excess kurtosis %8s  ->  %8s\n", fmt(kr - 3), fmt(kl - 3)))
tame <- is.finite(sl) && is.finite(kl) && abs(sl) < 1 && abs(kl - 3) < 1
if (tame) {
  cat("\n  *** THE LOG SCALE IS CLOSE TO GAUSSIAN on these two moments. A reader\n")
  cat("  *** is entitled to ask why a robust estimator is needed at all once the\n")
  cat("  *** transform is applied. CONCEDE IT, and make the argument where it\n")
  cat("  *** still stands: on the raw scale, where the transform is a modelling\n")
  cat("  *** choice a reader may not want to make, and in the simulation, where\n")
  cat("  *** contamination is not something a transform removes. ***\n")
} else {
  cat("\n  The log scale remains non-Gaussian on at least one moment, so the\n")
  cat("  transform is not a substitute for a robust estimator here. State it in\n")
  cat("  exactly those terms -- the transform helps and does not suffice -- and\n")
  cat("  put F1's raw-scale panel beside this one.\n")
}

# --- 3. what the estimator actually sees --------------------------------------
cat("\nPANEL 3. THE ESTIMATED FPC SCORES\n")
cat("  The marginal above describes the measurements. The estimator's exposure is\n")
cat("  the SCORE distribution, and the robustness argument is about that, so the\n")
cat("  two are reported separately rather than conflated.\n\n")
if (is.null(sc) || !nrow(sc)) {
  cat("  NOT AVAILABLE. This panel needs F1's fitted bandwidth, so F5 was run\n")
  cat(sprintf("  before %s existed or could not read it.\n", file.path(f1dir, "F1.rds")))
  cat("  Re-run F5 after F1 succeeds:\n")
  cat("    Rscript R/run_realdata36_f5.R <f1dir> <f5dir>\n")
} else {
  cat(sprintf("  %-8s %-10s %10s %12s %14s\n", "scale", "component", "sd",
              "skewness", "excess kurt"))
  for (i in seq_len(nrow(sc)))
    cat(sprintf("  %-8s %-10d %10s %12s %14s\n",
                if (i == 1L || sc$scale[i] != sc$scale[i - 1L]) sc$scale[i] else "",
                sc$component[i], fmt(sc$sd[i], 3), fmt(sc$skew[i]),
                fmt(sc$kurt[i] - 3)))
  # Counted WITHIN a scale, not pooled: a pooled count mixes the two answers
  # Panel 2 just separated, and the whole question there is whether the log
  # scale behaves differently from the raw one.
  cat("\n")
  for (s in unique(sc$scale)) {
    z <- sc[sc$scale == s & is.finite(sc$kurt), , drop = FALSE]
    if (!nrow(z)) next
    cat(sprintf("  %-4s: %d of %d components with excess kurtosis above 1\n",
                s, sum(z$kurt - 3 > 1), nrow(z)))
  }
  zr <- sc[sc$scale == "raw" & is.finite(sc$kurt), , drop = FALSE]
  zl <- sc[sc$scale == "log" & is.finite(sc$kurt), , drop = FALSE]
  if (nrow(zr) && any(zr$kurt - 3 > 1))
    cat("  Heavy tails in the raw-scale scores are the DIRECT form of the argument:\n  they sit in the quantity the estimator is exposed to, not merely in the\n  measurements.\n")
  else if (nrow(zr))
    cat("  No raw-scale score component carries excess kurtosis above 1. The heavy\n  tails are in the measurements rather than in the estimated scores, which is\n  a weaker version of the argument and must be written as such.\n")
  if (nrow(zl) && !any(zl$kurt - 3 > 1))
    cat("  The log-scale scores are not heavy-tailed. Read this beside Panel 2 --\n  together they say the transform does most of the work on THIS cohort.\n")
}

# --- 4. the formal tests, and why they are nearly useless here ----------------
cat("\nPANEL 4. THE FORMAL TESTS\n")
cat("  The Shapiro-Wilk results are printed by the RUNNER, not banked here, and\n")
cat("  that is deliberate: at this n any normality test rejects on a trivial\n")
cat("  departure, so the p-value carries no information about whether the\n")
cat("  departure matters. They are reported in the log so a reader knows they\n")
cat("  were run rather than hidden, and the moments above are what gets quoted.\n")

cat("\n== WHAT THIS LICENSES ==\n")
if (is.finite(sr) && is.finite(kr))
  cat(sprintf("   Pooled CEA measurements had skewness %s and excess kurtosis %s on the\n   raw scale (%s and %s after log transformation), which motivates a robust\n   estimator.\n",
              fmt(sr), fmt(kr - 3), fmt(sl), fmt(kl - 3)))
cat("   MOTIVATES. Not requires, not confirms, not justifies. The method assumes\n")
cat("   nothing about Gaussianity, and at a median of 7 observations per subject\n")
cat("   this cohort could not settle the question if it did.\n")
