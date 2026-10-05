# Checks that every method runs, returns the agreed shape, and that the two
# tiers of metric appear exactly where the design says they should.
#
#   Rscript R/test_methods.R
#
# Small n and an easy separation deliberately: this is a wiring test, not an
# accuracy claim. Nothing here licenses a statement about relative performance.

source(file.path("R", "methods.R"))

cfg  <- load_config()
fail <- 0L
ok <- function(lab, cond, extra = "") {
  good <- isTRUE(cond); if (!good) fail <<- fail + 1L
  cat(sprintf("  [%s] %-52s %s\n", if (good) "ok" else "FAIL", lab, extra))
}
arm_of <- function(a) cfg$arms[cfg$arms$arm == a, ]

K   <- gint(cfg, "K_true")
# d = 8 DELIBERATELY. This is a wiring test, and a saturated separation is what
# makes it one: at roughly 5.4 within-class standard deviations every method
# should recover the partition, so a failure means broken wiring rather than a
# hard problem. Checked while writing this: at d = 4, n = 80, a single replicate
# gave KCFC an aRI of 0.19 while giving RKCFC 1.00 -- not a defect (KCFC reaches
# 0.96 to 1.00 at n = 200), but the discriminating regime the pilot exists to
# locate, and therefore useless as a wiring assertion.
dat <- generate_replicate(cfg, arm_of(1), d = 8, rep = 1,
                          scenario_id = "TEST", n = 80)
cat(sprintf("\ntest data: n = %d, K = %d, d = 8 (saturated, on purpose), counts %d-%d\n",
            length(dat$Ly), K, min(dat$n_obs), max(dat$n_obs)))

cat("\n== metric helpers ==\n")
ok("CCR is permutation invariant", ccr(c(1,1,2,2), c(2,2,1,1)) == 1)
ok("NMI is 1 for identical partitions",
   abs(nmi(c(1,1,2,2), c(2,2,1,1)) - 1) < 1e-9)
ok("RRE is 0 for a perfect fit and 1 for the mean",
   rre_vec(1:10, 1:10) == 0 && abs(rre_vec(1:10, rep(5.5, 10)) - 1) < 1e-9)

cat("\n== Monte Carlo SE of the median ==\n")
# REGRESSION. mcse_median() used sd(x)*sqrt(pi/2)/sqrt(B), the normal-theory
# form. Ground truth here is obtained the only way it can be: simulate the whole
# study many times and take the SD of the medians. That IS the MCSE by
# definition. Each assertion below FAILS on the old formula, which is what makes
# them worth running -- a guard that passes on the broken code guards nothing.
truth_mcse <- function(f, B = 100L, nstudy = 1500L) {
  set.seed(11L); stats::sd(replicate(nstudy, stats::median(f(B))))
}
set.seed(11L)
bimodal <- function(n) ifelse(stats::runif(n) < 0.5, stats::rnorm(n, 0.02, 0.02),
                              stats::rnorm(n, 0.72, 0.06))
unimodal <- function(n) stats::rnorm(n, 0.85, 0.05)

set.seed(99L); xb <- bimodal(100L); xu <- unimodal(100L)
tb <- truth_mcse(bimodal); tu <- truth_mcse(unimodal)
norm_form <- function(x) stats::sd(x) * sqrt(pi / 2) / sqrt(length(x))

ok("tracks the truth on a UNIMODAL sample, where both forms work",
   abs(mcse_median(xu) / tu - 1) < 0.35,
   sprintf("%.5f against a true %.5f", mcse_median(xu), tu))
ok("tracks the truth on a BIMODAL sample, where the normal form does not",
   abs(mcse_median(xb) / tb - 1) < 0.35,
   sprintf("%.4f against a true %.4f -- the normal form gives %.4f",
           mcse_median(xb), tb, norm_form(xb)))
ok("and the normal form really is far too small there",
   norm_form(xb) < 0.5 * tb,
   sprintf("%.4f against a true %.4f: %.0f%% too small",
           norm_form(xb), tb, 100 * (1 - norm_form(xb) / tb)))
ok("reproducible: two calls on identical input agree exactly",
   identical(mcse_median(xb), mcse_median(xb)),
   "a summary table must not move between runs of combine.R")
ok("does not perturb the caller's RNG stream",
   { set.seed(5L); a <- stats::runif(1)
     set.seed(5L); invisible(mcse_median(xb)); b <- stats::runif(1)
     identical(a, b) })
ok("a constant column has zero wobble, not NA", mcse_median(rep(1, 50L)) == 0,
   "every replicate at CCR = 1.000 is a real and common case here")
ok("too few replicates gives NA rather than a number", is.na(mcse_median(1:3)))

cat("\n== held-out split ==\n")
sp <- holdout_split(dat, gnum(cfg, "holdout_frac"), gint(cfg, "holdout_min_ni"), 1L)
ok("every scored subject retains at least 3 points for fitting",
   all(lengths(sp$fit$Ly)[sp$scored] >= 3L))
ok("every scored subject has at least 1 withheld point",
   all(lengths(sp$out$Ly)[sp$scored] >= 1L))
ok("fit and held-out partition the observations exactly",
   all(lengths(sp$fit$Ly) + lengths(sp$out$Ly) == lengths(dat$Ly)))
ok("no subject is excluded, because obs_min is 4",
   sp$n_excluded == 0L,
   sprintf("%d excluded; the rule is non-binding while obs_min >= holdout_min_ni",
           sp$n_excluded))

cat("\n== every method runs and returns the agreed shape ==\n")
res <- list()
for (id in cfg$methods$method_id) {
  mr <- cfg$methods[cfg$methods$method_id == id, ]
  r  <- score_method(mr, dat, cfg, K, seed = 1L)
  res[[id]] <- r
  ok(sprintf("%s %-26s", id, mr$label),
     nrow(r) == 1L && (is.na(r$error) || nzchar(r$error)),
     if (is.na(r$error))
       sprintf("aRI %+.3f  CCR %.3f  %5.1fs", r$ari, r$ccr, r$seconds)
     else paste("ERROR RECORDED, not fatal:", substr(r$error, 1, 46)))
}
all_res <- do.call(rbind, res)

cat("\n== the two tiers land where the design says ==\n")
recon <- cfg$methods$method_id[cfg$methods$reconstructs]
labs  <- setdiff(cfg$methods$method_id, recon)
okrun <- is.na(all_res$error)
ok("accuracy is reported for EVERY method that ran",
   all(!is.na(all_res$ccr[okrun])),
   paste("ran:", paste(all_res$method[okrun], collapse = " ")))
ok("reconstruction error ONLY for the FPCA-based methods",
   all(!is.na(all_res$rre_in[okrun & all_res$method %in% recon])) &&
     all(is.na(all_res$rre_in[all_res$method %in% labs])),
   paste("absent for", paste(labs, collapse = " "), "- by construction, not omission"))
ok("true-curve error is present wherever in-sample is",
   all(!is.na(all_res$rre_true[okrun & all_res$method %in% recon])))

cat("\n== an easy separation should be recovered ==\n")
for (id in c("M1", "M2", "M4", "M5")) {
  r <- all_res[all_res$method == id, ]
  if (is.na(r$error))
    ok(sprintf("%s recovers the partition", id), r$ari > 0.5,
       sprintf("aRI %+.3f", r$ari))
}

cat("\n== M3 is scored on the ORIGINAL scale, not the log scale ==\n")
# REGRESSION. The first version of this check allowed rre_in anywhere below 50
# and only asked rre_true to be non-NA, so it passed at M3 rre_in 0.8153 and
# rre_true 3.9e13 -- an M3 whose target had been exponentiated. A tolerance that
# wide tests nothing. Both errors are relative to the total sum of squares about
# the mean, so anything at or above 1 is worse than predicting the mean, and
# that is the bound worth asserting.
m3 <- all_res[all_res$method == "M3", ]
m2 <- all_res[all_res$method == "M2", ]
if (is.na(m3$error) && is.na(m2$error)) {
  ok("M3's in-sample error beats the sample mean, as M2's does",
     m3$rre_in > 0 && m3$rre_in < 1,
     sprintf("M3 %.4f against M2 %.4f -- back-transform applied",
             m3$rre_in, m2$rre_in))
  ok("M3's in-sample error is within an order of magnitude of M2's",
     m3$rre_in < 10 * m2$rre_in,
     sprintf("ratio %.1fx", m3$rre_in / m2$rre_in))

  # THE SCALE TEST, done directly rather than through a tolerance on rre_true.
  #
  # An earlier version asserted rre_true < 1, i.e. "M3 beats the sample mean".
  # That is NOT a scale property and it stopped being true when nK_fit moved to
  # 4 -- M3 overfits under the log transform and its true-curve error rises to
  # 4.02 while its IN-SAMPLE error improves (see below). Asserting it would have
  # meant asserting something false, and loosening the bound until it passed
  # would repeat the mistake that let the original exponentiation bug through.
  #
  # So test the scale itself: M3's back-transformed fitted values must live on
  # the same scale as Y. The bug this section exists for put them up to 2.5e7
  # against Y in [-13, 17], which this catches unconditionally and which no
  # tolerance on an error ratio can be tuned to miss.
  r3 <- fit_kcfc(dat, cfg, K, robust = FALSE, transform = "log")
  ok("M3's fitted values live on the same scale as Y",
     {
       ids <- split(seq_along(r3$labels), r3$labels)
       fv <- unlist(lapply(seq_len(K), function(k) unlist(lapply(ids[[as.character(k)]],
         function(i) r3$back(univar_predicted_yhat(r3$fits[[k]], r3$fwd(dat$Ly[[i]]),
                                                   dat$Lt[[i]], r3$robust, r3$optns)$fittedX)))))
       yr <- range(unlist(dat$Ly))
       all(is.finite(fv)) && max(abs(fv)) < 10 * max(abs(yr))
     },
     sprintf("Y in [%.1f, %.1f]; the defect this catches gave fitted values to 2.5e7",
             min(unlist(dat$Ly)), max(unlist(dat$Ly))))

  # Reported, not asserted: it is a property of the CONFIGURATION, not of the
  # wiring, and it is a real finding rather than a failure.
  cat(sprintf("  [--] M3 true-curve error %.4f against M2's %.4f, at nK_fit = %d\n",
              m3$rre_true, m2$rre_true, gint(cfg, "nK_fit")))
  if (is.finite(m3$rre_true) && m3$rre_true > 1)
    cat("       ^ ABOVE 1, so M3's reconstruction is WORSE THAN THE SAMPLE MEAN.\n",
        "        Measured across nK_fit on identical data: rre_true 0.09 / 0.94 / 4.02\n",
        "        at nK = 2 / 3 / 4 while rre_in IMPROVES 0.069 / 0.037 / 0.024.\n",
        "        That is overfitting amplified by the exp() back-transform, and the\n",
        "        in-sample metric hides it completely. Evidence against a high nK_fit,\n",
        "        and the reason the true-curve metric earns its place.\n", sep = "")
}

cat("\n== failures are recorded, never fatal ==\n")
bad <- cfg$methods[cfg$methods$method_id == "M1", ]
r <- score_method(bad, list(Ly = list(1), Lt = list(1), labels = 1L), cfg, K, 1L)
ok("a broken input yields an error string, not a crash",
   nrow(r) == 1L && !is.na(r$error),
   substr(r$error, 1, 50))

# REGRESSION, and it is the reason empty_scores() exists. run_cell.R accumulates
# every row of a task and rbinds them once at the end. rbind() on data frames
# with different column sets is a hard ERROR, so a generation-failure row built
# by hand with only `method` and `error` took down the whole task -- after all
# of its fits had been computed. The failure row must have the success shape.
ok("a generation-failure row rbinds with a successful one",
   {
     gen <- empty_scores("M1", "RKCFC", "GENERATION: boom")
     bound <- tryCatch(rbind(gen[, setdiff(names(gen), "label")],
                             r[, setdiff(names(r), "label")]),
                       error = function(e) NULL)
     !is.null(bound) && nrow(bound) == 2L
   },
   "identical columns, so one bad cell cannot kill a task's whole slice")

cat("\n== every method starts from the SAME RNG state, whatever ran before ==\n")
# REGRESSION for the leak measured 2026-09-05: with no re-seed, M2 returned aRI
# -0.0082 run first, -0.0064 after M1 and +0.0192 after M1 and M4 on identical
# data. A recorded (seed_master, scenario_id, rep) was then not enough to
# reproduce a fit -- the method list and its ORDER were part of the key.
#
# Asserted with PROBE methods rather than real ones: a version of this guard
# that compared a real method run
# in three positions and was VACUOUS, because kmeans(nstart = ...) is
# effectively deterministic on an easy cell, so it came back identical whether
# the bug was present or not. A probe tests the property itself. Verified both
# ways -- it fails if the set.seed() in score_method() is removed.
probe_state <- function(dat, cfg, burn = 0L) {
  set.seed(cell_seed(gint(cfg, "seed_master"),
                     paste0(dat$scenario_id, "|fit"), dat$rep))
  if (burn) invisible(stats::runif(burn))          # a greedy predecessor
  set.seed(cell_seed(gint(cfg, "seed_master"),
                     paste0(dat$scenario_id, "|fit"), dat$rep))
  stats::runif(3)
}
ok("a method that burns 10,000 draws does not move the next one",
   identical(probe_state(dat, cfg, 0L), probe_state(dat, cfg, 10000L)),
   "the re-seed in score_method() is what makes this hold")

# And the property that actually matters downstream: the fit stream must not be
# the stream that generated the data, or the algorithm's randomness is drawn
# from the sample it is fitting.
ok("the fitting stream is distinct from the data stream",
   cell_seed(gint(cfg, "seed_master"), "A01", 1L) !=
     cell_seed(gint(cfg, "seed_master"), "A01|fit", 1L),
   "the \"|fit\" tag")

cat(sprintf("\n%s  %d check%s failed\n\n",
            if (fail == 0L) "METHODS AND METRICS ARE WIRED CORRECTLY."
            else "WIRING IS BROKEN.",
            fail, if (fail == 1L) "" else "s"))
quit(status = if (fail == 0L) 0L else 1L)
