# Reads the configuration tables and checks they are internally consistent.
# Run from the study/ directory:
#
#   Rscript config/validate_config.R
#
# This is the ONLY place the tables are interpreted. Simulation code calls
# load_config() rather than hard-coding any of these values, so that changing a
# parameter never means editing code.
#
# The checks below are not decoration. Three classes of error in this design are
# silent -- they produce plausible numbers rather than an error -- and each has
# a check here:
#
#   1. The design's arm numbering is NOT the package's gen_scores numbering.
#      Design arm 4 (t nu = 1.5) is gen_scores arm 6. Getting this wrong would
#      simulate a contaminated normal and label it infinite-variance.
#   2. The basis must be orthonormal ON THE GRID, not merely analytically. A
#      coarse grid breaks orthogonality even when the closed forms are exact.
#   3. The overlap between each class mean and its own eigenfunctions is
#      REPORTED rather than asserted. It is not zero in the published design:
#      class 1 is orthogonal, class 2 is not. See the diagnostic below.

CFG <- "config"

load_config <- function(dir = CFG) {
  g <- utils::read.csv(file.path(dir, "params_global.csv"), stringsAsFactors = FALSE)
  vals <- stats::setNames(as.list(g$value), g$param)
  rd <- function(f) utils::read.csv(file.path(dir, f), stringsAsFactors = FALSE)
  list(
    global      = vals,
    meta        = g,
    arms        = rd("arms.csv"),
    basis       = rd("basis.csv"),
    Z           = rd("Z.csv"),
    methods     = rd("methods.csv"),
    estimators  = rd("estimators.csv"),
    pilot       = rd("scenarios_p.csv"),
    core        = rd("scenarios_a.csv"),
    variations  = rd("scenarios_b.csv"),
    krecovery   = rd("scenarios_k.csv")
  )
}

# Numeric scalar from the global table, which stores everything as character.
gnum <- function(cfg, key) as.numeric(cfg$global[[key]])
gint <- function(cfg, key) as.integer(cfg$global[[key]])
glog <- function(cfg, key) as.logical(cfg$global[[key]])
gvec <- function(cfg, key) as.numeric(strsplit(cfg$global[[key]], ",")[[1]])

# Evaluate a basis expression on a time grid. `d` is used only by mean functions.
eval_basis <- function(expr, t, d = 0) eval(parse(text = expr), list(t = t, d = d))

# Score dispersion for one class at one M_true, as a diagonal matrix.
Z_matrix <- function(cfg, class, M_true) {
  row <- cfg$Z[cfg$Z$class == class & cfg$Z$M_true == M_true, ]
  stopifnot(nrow(row) == 1L)
  v <- as.numeric(row[, paste0("m", seq_len(M_true))])
  diag(v, nrow = M_true)
}

# The eigenfunction role name for a given M_true, since the 4-component
# extension lives in its own rows rather than overwriting the published pair.
eigen_role <- function(M_true) if (M_true == 2L) "eigen" else paste0("eigen", M_true)

# Eigenfunctions for one class, evaluated on a grid: an n_grid x M matrix.
eval_eigen <- function(cfg, class, M_true, t, contrast = "parallel") {
  rows <- cfg$basis[cfg$basis$contrast == contrast &
                      cfg$basis$class == class &
                      cfg$basis$role == eigen_role(M_true), ]
  rows <- rows[order(rows$component), ]
  stopifnot(nrow(rows) == M_true)
  vapply(rows$expr, eval_basis, numeric(length(t)), t = t, USE.NAMES = FALSE)
}

# Mean function for one class on a grid.
eval_mean <- function(cfg, class, t, d, contrast = "parallel") {
  row <- cfg$basis[cfg$basis$contrast == contrast &
                     cfg$basis$class == class &
                     cfg$basis$role == "mean", ]
  stopifnot(nrow(row) == 1L)
  eval_basis(row$expr, t, d = d)
}

# ---------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  cfg  <- load_config()
  fail <- 0L
  ok <- function(lab, cond, extra = "") {
    good <- isTRUE(cond)
    if (!good) fail <<- fail + 1L
    cat(sprintf("  [%s] %-58s %s\n", if (good) "ok" else "FAIL", lab, extra))
  }

  cat("\n== global parameters ==\n")
  ok("n is a positive integer", gint(cfg, "n") > 0)
  ok("K_true == 2", gint(cfg, "K_true") == 2L)
  ok("class proportions sum to 1", abs(sum(gvec(cfg, "pi")) - 1) < 1e-12)
  ok("obs_min >= 2 (predicted_yhat cannot take a singleton)", gint(cfg, "obs_min") >= 2L)
  ok("LOO is FALSE", identical(glog(cfg, "LOO"), FALSE))
  ok("bw is numeric (KFPCA has NO GCV path for it)", is.finite(gnum(cfg, "bw")))
  ok("recon == 'own' (no common yardstick imposed)", cfg$global$recon == "own")
  ok("separations are set and ordered",
     gnum(cfg, "d_hard") < gnum(cfg, "d_easy"),
     sprintf("hard %.1f, easy %.1f - set on the UNIVARIATE Delta sweep, 2026-09-05",
             gnum(cfg, "d_hard"), gnum(cfg, "d_easy")))

  cat("\n== arm table, and the mapping that is easy to get wrong ==\n")
  a <- cfg$arms
  ok("6 arms", nrow(a) == 6L)
  ok("design arm 4 (t nu=1.5) maps to gen_scores arm 6",
     a$score_arm[a$arm == 4L] == 6L, "NOT arm 4, which is the contaminated normal")
  ok("arms 5 and 6 use Gaussian scores",
     all(a$score_arm[a$arm %in% c(5L, 6L)] == 1L),
     "their departure is in the errors / curves")
  ok("exactly one arm has t3 errors", sum(a$error_dist == "t3") == 1L)
  ok("exactly one arm contaminates curves", sum(a$contaminate) == 1L)

  cat("\n== basis, evaluated ON THE GRID and not merely analytically ==\n")
  ng <- gint(cfg, "grid_n")
  t  <- seq(gnum(cfg, "t_min"), gnum(cfg, "t_max"), length.out = ng)
  # TRAPEZOID, not mean(). A grid from seq(0, 1, length.out = n) carries both
  # endpoints, so mean() double-counts them and understates every norm by
  # exactly 1/n -- 0.0099 at n = 101, which looks like a basis defect and is not.
  ip <- function(x, y) { z <- x * y; sum(diff(t) * (utils::head(z, -1) + utils::tail(z, -1)) / 2) }
  gram <- function(P) { M <- ncol(P); G <- matrix(0, M, M)
    for (i in seq_len(M)) for (j in seq_len(M)) G[i, j] <- ip(P[, i], P[, j]); G }
  for (M in c(2L, 4L)) {
    for (k in 1:2) {
      P <- eval_eigen(cfg, k, M, t)
      G <- gram(P)
      ok(sprintf("class %d, M_true = %d: orthonormal on the grid", k, M),
         max(abs(G - diag(M))) < 1e-6,
         sprintf("max |G - I| = %.2e", max(abs(G - diag(M)))))
    }
  }
  # Class 3 exists at M_true = 2 only, for the K-recovery block. It is checked
  # here rather than left to fail inside a runner at 3 a.m. on the cluster.
  P3 <- eval_eigen(cfg, 3L, 2L, t); G3 <- gram(P3)
  ok("class 3, M_true = 2: orthonormal on the grid",
     max(abs(G3 - diag(2))) < 1e-6,
     sprintf("max |G - I| = %.2e -- the K-recovery third class", max(abs(G3 - diag(2)))))
  ok("class 3 has its own Z row", nrow(cfg$Z[cfg$Z$class == 3L, ]) == 1L)
  # The three means must sit at -d, 0, +d, so the two ADJACENT gaps both equal d
  # and the hardest pair at K = 3 is exactly as hard as the only pair at K = 2.
  # Without that, a K-recovery result would confound "can it find three classes"
  # with "are three classes easier to find".
  m <- vapply(1:3, function(k) mean(eval_mean(cfg, k, t, d = 2)), numeric(1))
  ok("the three class means are equally spaced at d",
     abs(abs(m[2] - m[1]) - 2) < 1e-6 && abs(abs(m[1] - m[3]) - 2) < 1e-6,
     sprintf("levels %.2f, %.2f, %.2f at d = 2", m[1], m[2], m[3]))

  cat("\n== mean vs eigenfunctions: a property of the published design, reported ==\n")
  for (k in 1:2) {
    P  <- eval_eigen(cfg, k, 2L, t)
    mu <- eval_mean(cfg, k, t, d = 0)
    ov <- vapply(seq_len(ncol(P)), function(j) ip(P[, j], mu), numeric(1))
    cat(sprintf("  [--] class %d mean vs its own eigenfunctions            |<mu, phi_m>| = %s\n",
                k, paste(sprintf("%.4f", abs(ov)), collapse = ", ")))
  }
  cat("       Class 1 is orthogonal; class 2 is NOT, because\n")
  cat("       <sqrt(2)sin(4pi t), sqrt(2)cos(3pi t)> = 2/(7pi) + 2/pi = 0.7276.\n")
  cat("       This is not an error: a Karhunen-Loeve model does not require mu\n")
  cat("       orthogonal to phi, and the mean is removed before the FPCA step.\n")
  cat("       But the two classes are NOT treated symmetrically by the published\n")
  cat("       design, and for class 2 the mean and the second component direction\n")
  cat("       are heavily aligned, so error in one leaks into the other. This\n")
  cat("       matters for the estimation study, whose mean-error and eigenfunction\n")
  cat("       endpoints are therefore correlated for class 2 and not for class 1.\n")

  cat("\n== eigen-gap: individually identifiable ==\n")
  for (k in 1:2) {
    v <- diag(Z_matrix(cfg, k, 2L))
    ok(sprintf("class %d relative eigen-gap > 40%%", k),
       (v[1] - v[2]) / v[1] > 0.4,
       sprintf("%.1f%%",
               100 * (v[1] - v[2]) / v[1]))
  }

  cat("\n== crossing means really do have equal average level ==\n")
  m1 <- eval_mean(cfg, 1, t, d = 2, contrast = "crossing")
  m2 <- eval_mean(cfg, 2, t, d = 2, contrast = "crossing")
  ok("average levels equal", abs(mean(m1) - mean(m2)) < 5e-3,
     sprintf("difference %.2e", abs(mean(m1) - mean(m2))))
  ok("the curves actually cross", any(diff(sign(m1 - m2)) != 0))

  cat("\n== mean-identical contrast is genuinely identical ==\n")
  i1 <- eval_mean(cfg, 1, t, d = 5, contrast = "identical")
  i2 <- eval_mean(cfg, 2, t, d = 5, contrast = "identical")
  ok("means coincide even when d is nonzero", max(abs(i1 - i2)) < 1e-12,
     "so only the eigenfunctions differ")

  cat("\n== methods and scenarios ==\n")
  m <- cfg$methods
  ok("7 methods", nrow(m) == 7L)
  ok("3 methods reconstruct (M1, M2, M3)",
     identical(sort(m$method_id[m$reconstructs]), c("M1", "M2", "M3")),
     "M4, M5, M7 return labels only; M6 fits a different model family")
  ok("pilot uses only M1 and M2", sum(m$include_pilot) == 2L)
  ok("pilot has 6 scenarios", nrow(cfg$pilot) == 6L)
  ok("pilot includes d = 0 on both arms", sum(cfg$pilot$d == 0) == 2L,
     "the shape-only floor")
  # The pilot hardcodes its d values while the core grid resolves "hard"/"easy"
  # through d_of(). They therefore drift apart silently: when d_easy moved 3 -> 2
  # on 2026-09-05 the pilot was left confirming a value the core grid no longer
  # used, and never touching the one it did. Caught by hand in the design audit;
  # this makes it impossible to miss again.
  ok("every pilot d is 0 or a configured design level",
     all(cfg$pilot$d %in% c(0, gnum(cfg, "d_hard"), gnum(cfg, "d_easy"))),
     sprintf("pilot has %s; design levels are %g and %g",
             paste(sort(unique(cfg$pilot$d)), collapse = ", "),
             gnum(cfg, "d_hard"), gnum(cfg, "d_easy")))
  ok("core grid has 12 scenarios", nrow(cfg$core) == 12L)
  ok("core covers 6 arms x 2 separations",
     nrow(unique(cfg$core[, c("arm", "d_level")])) == 12L)
  ok("K recovery has 3 scenarios", nrow(cfg$krecovery) == 3L)
  ok("K recovery spans true K = 1, 2, 3",
     identical(sort(cfg$krecovery$K_true), 1:3))
  ok("variations have a baseline row", any(cfg$variations$factor == "baseline"),
     "needed so each panel has a reference at matched methods and replicates")
  ok("mean-identical row uses the identical contrast",
     cfg$variations$contrast[cfg$variations$level == "mean-identical"] == "identical")

  cat("\n== fit budget ==\n")
  nb <- nrow(cfg$variations) * 2L
  cat(sprintf("  scenarios: pilot %d + core %d + variations %d + K %d = %d\n",
              nrow(cfg$pilot), nrow(cfg$core), nb, nrow(cfg$krecovery),
              nrow(cfg$pilot) + nrow(cfg$core) + nb + nrow(cfg$krecovery)))
  cat(sprintf("  fits:      pilot %d + core %d + variations %d + K %d\n",
              nrow(cfg$pilot) * 2L * gint(cfg, "n_reps_pilot"),
              nrow(cfg$core) * nrow(cfg$methods) * gint(cfg, "n_reps_core"),
              nb * 4L * gint(cfg, "n_reps"),
              nrow(cfg$krecovery) * 5L * 2L * gint(cfg, "n_reps")))

  cat(sprintf("\n%s  %d check%s failed\n\n",
              if (fail == 0L) "ALL CHECKS PASSED." else "CONFIGURATION IS INCONSISTENT.",
              fail, if (fail == 1L) "" else "s"))
  quit(status = if (fail == 0L) 0L else 1L)
}
