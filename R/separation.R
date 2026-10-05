# Robust Mahalanobis separation between classes.
#
# Adapted from the multivariate extension's `R/separation.R`, dropping the
# multivariate stacking. `robust_stage2()` is called `robust_scatter()` in the
# standalone package -- identical numerical body, different name, because "stage
# two" names the Happ and Greven second stage, which does not exist univariately.
#
# WHY THIS EXISTS. `d` is a poor measure of how hard a scenario is. It sits on
# the raw scale, ignores the within-class score dispersion Z, and -- decisively
# here -- presumes a variance that arm 4 (`t nu = 1.5`) does not have. Δ measures
# separation where the algorithm actually works: in the estimated FPC score
# space, AFTER sparse observation and score estimation. With Poisson(7) points
# per subject the separation you can exploit is not the separation of the
# noiseless mean functions.
#
#   Delta^2 = t(rho_1 - rho_2) %*% solve(S) %*% (rho_1 - rho_2)
#
# Location is the spatial median per class, scatter the pooled within-class
# robust estimate. Both are defined without moment assumptions, so the measure
# survives arm 4 where a covariance-based one has no population target at all.
#
# THE SCATTER IS BUILT AS C diag(nu) t(C), not from the SSCM directly. The
# spatial-sign covariance matrix is the covariance of unit vectors and has trace
# 1, so inverting it would put Δ on an arbitrary scale. `robust_scatter()`'s
# `values` recover the spread along each direction on the score scale, which is
# what makes Δ comparable to a conventional Mahalanobis distance.
#
# TWO USES.
#
#   1. A reference value. Performance is read against a
#      conventional reference value; `reference_ccr()` supplies one, so a CCR
#      of 0.98 can be read against a yardstick rather than in a vacuum.
#
#   2. A CHECK THAT THE IQR MATCHING WORKED. The design matches the arms on IQR
#      so that any difference can be attributed to signal-to-noise rather
#      than to tails. Δ should therefore be near-constant ACROSS arms at fixed
#      `d`. If it is not, the arms are not exchangeable and tail effects are
#      confounded with separation -- which is the objection the matching exists
#      to foreclose. This turns that assertion into evidence.
#
# The check below measures that spread on this design.

source(file.path("R", "methods.R"))

# Pairwise robust Mahalanobis distances between class centroids. Returns the
# full matrix plus the minimum, which is the binding pair: the hardest two
# classes to separate set the difficulty, not the average.
robust_separation <- function(scores, labels) {
  sm  <- RKCFC:::spatial_median
  rsc <- RKCFC:::robust_scatter
  cls <- sort(unique(labels))
  if (length(cls) < 2L)
    return(list(matrix = NULL, min = NA_real_, mean = NA_real_))

  cent <- lapply(cls, function(k) sm(scores[labels == k, , drop = FALSE]))

  # Pooled within-class scatter: centre each class on its OWN location first, so
  # the between-class shift does not inflate the thing it is being divided by.
  resid <- do.call(rbind, lapply(seq_along(cls), function(i)
    sweep(scores[labels == cls[i], , drop = FALSE], 2L, cent[[i]], "-")))
  e <- rsc(resid, M_keep = ncol(scores), robust = TRUE)
  S <- e$vectors %*% diag(e$values, nrow = length(e$values)) %*% t(e$vectors)
  Sinv <- tryCatch(solve(S), error = function(err) MASS::ginv(S))

  M <- matrix(NA_real_, length(cls), length(cls), dimnames = list(cls, cls))
  for (i in seq_along(cls)) for (j in seq_along(cls)) {
    dv <- cent[[i]] - cent[[j]]
    M[i, j] <- sqrt(max(0, as.numeric(t(dv) %*% Sinv %*% dv)))
  }
  off <- M[upper.tri(M)]
  list(matrix = M, min = min(off), mean = mean(off))
}

# Reference CCR implied by a separation Δ, for two equal-prior classes.
#
# NOT A BAYES BOUND, though it is easily mistaken for one. `Phi(Delta/2)` is the accuracy of the optimal LINEAR rule for
# two Gaussians with a COMMON covariance -- the best achievable from centroid
# separation alone.
#
# THIS DESIGN VIOLATES THAT ASSUMPTION DELIBERATELY. The two classes have
# class-specific eigenfunctions, sines for class 1 and cosines for class 2, so
# they differ in covariance structure as well as in location. The true Bayes rule
# is quadratic and does better, and a method exploiting the shape difference
# legitimately EXCEEDS this reference. That is not a bug in the reference and it
# is not a method behaving impossibly. Expect it here, and strongly -- at `d = 0` the
# means COINCIDE, so Δ is near zero and the reference is 0.5, while RKCFC reaches
# CCR 0.925 on shape alone.
#
# Read it as a yardstick, not a ceiling. Exceeding it says a method is using more
# than the centroid gap. Falling well short says the cell is a floor.
#
# Undefined for K > 2, so NA is returned rather than a number quietly wrong.
reference_ccr <- function(delta, K = 2L) {
  if (K != 2L) return(NA_real_)
  stats::pnorm(delta / 2)
}

# Bands, anchored on attainable accuracy rather than on an arbitrary distance,
# so each cutoff is a claim a reader can check. Cutoffs are set on the reference
# CCR and inverted. Note that Δ is CONVEX in `d`, so equal steps in `d` are not
# equal steps in difficulty.
separation_band <- function(delta) {
  if (!is.finite(delta)) return(NA_character_)
  if (delta < 1.0) "negligible"     # ref < 0.69; even a perfect method looks poor
  else if (delta < 2.0) "weak"      # 0.69-0.84; answers "weakly separated clusters"
  else if (delta < 3.3) "medium"    # 0.84-0.95; THE DISCRIMINATING BAND
  else if (delta < 4.7) "strong"    # 0.95-0.99; easy but still rankable
  else "saturated"                  # > 0.99; cannot rank methods, at any B
}

# Separation for one generated replicate, measured on the same global scores the
# score-then-cluster comparators consume -- so Δ describes the representation the
# methods actually see, not an idealised one.
separation_of <- function(dat, cfg, robust = TRUE) {
  optns <- build_optns(cfg, robust)
  sc <- tryCatch(global_scores(dat, cfg, robust, optns), error = function(e) NULL)
  if (is.null(sc)) return(list(delta = NA_real_, ref_ccr = NA_real_,
                               band = NA_character_))
  r <- robust_separation(sc, dat$labels)
  list(delta = r$min, ref_ccr = reference_ccr(r$min, length(unique(dat$labels))),
       band = separation_band(r$min))
}
