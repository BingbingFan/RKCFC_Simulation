# Checks the robust Mahalanobis separation, and then RUNS the cross-arm check
# the design's IQR matching depends on.
#
#   Rscript R/test_separation.R [n_reps]
#
# The second half is not a unit test. It is the measurement the pilot exists to
# produce: Delta should be near-constant ACROSS arms at fixed d, because the arms
# are matched on IQR precisely so that any difference can be attributed to
# signal-to-noise rather than to tails. If it is not, the arms are not
# exchangeable and tail effects are confounded with separation.

source(file.path("R", "separation.R"))

cfg  <- load_config()
fail <- 0L
ok <- function(lab, cond, extra = "") {
  good <- isTRUE(cond); if (!good) fail <<- fail + 1L
  cat(sprintf("  [%s] %-54s %s\n", if (good) "ok" else "FAIL", lab, extra))
}

# A construction with a separation we know, so the measure can be checked
# against a number rather than against itself.
mk <- function(shift, n = 400L, p = 2L, seed = 1L) {
  set.seed(seed)
  X <- rbind(matrix(stats::rnorm(n * p), n, p),
             matrix(stats::rnorm(n * p), n, p) +
               matrix(rep(c(shift, rep(0, p - 1L)), each = n), n, p))
  list(X = X, lab = rep(1:2, each = n))
}

cat("\n== recovers a known separation ==\n")
for (s in c(2, 4)) {
  r <- robust_separation(mk(s)$X, mk(s)$lab)
  ok(sprintf("shift %g gives Delta about %g", s, s),
     abs(r$min - s) < 0.45, sprintf("Delta = %.3f", r$min))
}

cat("\n== behaves as a distance ==\n")
r <- robust_separation(mk(3)$X, mk(3)$lab)
ok("symmetric", isTRUE(all.equal(r$matrix, t(r$matrix))))
ok("zero on the diagonal", all(abs(diag(r$matrix)) < 1e-8))
ok("no separation gives Delta near 0",
   robust_separation(mk(0)$X, mk(0)$lab)$min < 0.35)

cat("\n== scale invariance: why Mahalanobis and not Euclidean ==\n")
# Stated COMPARATIVELY, because that is the actual claim and because an absolute
# tolerance here is tunable and therefore worthless. The spatial-sign scatter is
# orthogonally equivariant but NOT affine equivariant, so rescaling one
# coordinate does move Delta a little; what matters is that it moves it far less
# than a Euclidean centroid distance, which is swamped outright.
a <- mk(3); b <- a; b$X[, 2] <- b$X[, 2] * 100
euc <- function(z) {
  m <- vapply(sort(unique(z$lab)),
              function(k) colMeans(z$X[z$lab == k, , drop = FALSE]),
              numeric(ncol(z$X)))
  sqrt(sum((m[, 1] - m[, 2])^2))
}
mah_rel <- abs(robust_separation(b$X, b$lab)$min /
               robust_separation(a$X, a$lab)$min - 1)
euc_rel <- abs(euc(b) / euc(a) - 1)
ok("Mahalanobis moves far less than Euclidean under a 100x rescale",
   mah_rel < 0.25 && mah_rel < euc_rel / 5,
   sprintf("Mahalanobis %+.1f%% against Euclidean %+.1f%%", 100 * mah_rel, 100 * euc_rel))

cat("\n== defined without moments, which is the point on arm 4 ==\n")
set.seed(7)
n <- 400L
Xt <- rbind(matrix(stats::rt(n * 2, df = 1.5), n, 2),
            matrix(stats::rt(n * 2, df = 1.5), n, 2) +
              matrix(rep(c(3, 0), each = n), n, 2))
rt_ <- robust_separation(Xt, rep(1:2, each = n))
ok("finite Delta under an infinite-variance law", is.finite(rt_$min),
   sprintf("Delta = %.3f -- no covariance exists here", rt_$min))

cat("\n== the reference CCR, which is NOT a bound ==\n")
ok("Delta = 0 is chance", abs(reference_ccr(0) - 0.5) < 1e-9)
ok("Delta = 4 implies about 0.977", abs(reference_ccr(4) - 0.9772) < 1e-3)
ok("monotone in Delta", reference_ccr(5) > reference_ccr(3))
ok("NA for K > 2 rather than a number quietly wrong",
   is.na(reference_ccr(3, K = 3L)))
ok("bands land where the design says",
   separation_band(0.5) == "negligible" && separation_band(1.5) == "weak" &&
     separation_band(2.5) == "medium" && separation_band(4.0) == "strong" &&
     separation_band(6.0) == "saturated",
   "cutoffs are set on attainable accuracy, then inverted")

cat("\n== more classes: the minimum pair is what binds ==\n")
set.seed(3)
X3 <- rbind(matrix(stats::rnorm(200 * 2), 200, 2),
            matrix(stats::rnorm(200 * 2), 200, 2) + matrix(c(2, 0), 200, 2, byrow = TRUE),
            matrix(stats::rnorm(200 * 2), 200, 2) + matrix(c(9, 0), 200, 2, byrow = TRUE))
r3 <- robust_separation(X3, rep(1:3, each = 200))
ok("min is the hardest pair, not the average", r3$min < r3$mean,
   sprintf("min %.2f against mean %.2f", r3$min, r3$mean))

# --- the measurement, not a test ----------------------------------------------

nrep <- as.integer(c(commandArgs(TRUE), 3L)[1])
cat(sprintf("\n== Delta ACROSS ARMS at fixed d -- the IQR-matching check, %d reps ==\n",
            nrep))
cat("   Near-constant down each column is what the matching promises. A spread\n")
cat("   much above about 1.3x means the arms are NOT exchangeable at that d and\n")
cat("   tail effects are confounded with signal-to-noise.\n\n")
cat(sprintf("   %-24s %-22s %-22s\n", "arm", "d = 1.5", "d = 3"))
tab <- matrix(NA_real_, 6, 2)
for (arm in 1:6) {
  for (j in seq_along(c(1.5, 3))) {
    d <- c(1.5, 3)[j]
    v <- vapply(seq_len(nrep), function(r) {
      dat <- generate_replicate(cfg, cfg$arms[cfg$arms$arm == arm, ], d = d,
                                rep = r, scenario_id = sprintf("SEP%d_%g", arm, d))
      separation_of(dat, cfg, robust = TRUE)$delta
    }, numeric(1))
    tab[arm, j] <- stats::median(v, na.rm = TRUE)
  }
  cat(sprintf("   %-24s %-22s %-22s\n",
              substr(cfg$arms$label[cfg$arms$arm == arm], 1, 22),
              sprintf("%.2f  ref %.3f  %s", tab[arm, 1],
                      reference_ccr(tab[arm, 1]), separation_band(tab[arm, 1])),
              sprintf("%.2f  ref %.3f  %s", tab[arm, 2],
                      reference_ccr(tab[arm, 2]), separation_band(tab[arm, 2]))))
}
# REPORTED, NOT ASSERTED, and the distinction is deliberate. A test suite asserts
# that CODE is correct. Whether `d` was well chosen is a property of the DESIGN,
# and failing the suite on it would say the harness is broken when the harness is
# working exactly as intended -- it is telling us something we did not want to
# hear. The file header said so before the number was known; this keeps faith
# with that. The flag below is loud so it cannot be skimmed past.
cat("\n")
not_exchangeable <- character(0)
for (j in 1:2) {
  sp <- max(tab[, j], na.rm = TRUE) / min(tab[, j], na.rm = TRUE)
  d_lab <- c("1.5", "3")[j]
  lo <- which.min(tab[, j]); hi <- which.max(tab[, j])
  if (sp >= 1.4) not_exchangeable <- c(not_exchangeable, d_lab)
  cat(sprintf("   [%s] cross-arm spread at d = %-4s %.2fx   (%s at %.2f, %s at %.2f)\n",
              if (sp < 1.4) "ok" else "!!", d_lab, sp,
              cfg$arms$label[lo], tab[lo, j],
              cfg$arms$label[hi], tab[hi, j]))
  if (sp >= 1.4) {
    cat(sprintf("        THE ARMS ARE NOT EXCHANGEABLE AT d = %s. Reference CCR runs\n", d_lab))
    cat(sprintf("        %.3f to %.3f, so %d points of attainable accuracy come from\n",
                reference_ccr(tab[lo, j]), reference_ccr(tab[hi, j]),
                round(100 * (reference_ccr(tab[hi, j]) - reference_ccr(tab[lo, j])))))
    cat("        nothing but WHICH ARM you are in. That is the confound the IQR\n")
    cat("        matching exists to foreclose, and the manuscript asserts it holds.\n")
    cat("        Either report accuracy RELATIVE to each arm's own reference, or\n")
    cat("        recalibrate d per arm. Do not assert the matching unqualified.\n")
  }
}

# The summary must not contradict the report above it. An earlier version printed
# "THE ARMS ARE EXCHANGEABLE" directly beneath a flag saying they are not, because
# it summarised only the code assertions. A reader who skims to the last line is
# the one most in need of the finding.
cat(sprintf("\n%s  %d check%s failed",
            if (fail == 0L) "THE MEASURE IS CORRECT." else "THE MEASURE IS BROKEN.",
            fail, if (fail == 1L) "" else "s"))
if (length(not_exchangeable)) {
  cat(sprintf("\nBUT THE ARMS ARE NOT EXCHANGEABLE AT d = %s -- see the flag above.\n",
              paste(not_exchangeable, collapse = " and ")))
  cat("That is a finding about the DESIGN, not a failure of this code.\n\n")
} else {
  cat("  The arms are exchangeable at every d tested.\n\n")
}
quit(status = if (fail == 0L) 0L else 1L)
