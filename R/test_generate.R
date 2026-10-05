# Checks the generator against what the DESIGN says it should produce.
#
#   Rscript R/test_generate.R
#
# Measures can run without error and still measure the wrong thing, so the
# standard here is:
# a check that only confirms "no crash" confirms nothing. Every assertion below
# compares against a number the design specifies in advance.

source(file.path("R", "generate.R"))

cfg  <- load_config()
fail <- 0L
ok <- function(lab, cond, extra = "") {
  good <- isTRUE(cond); if (!good) fail <<- fail + 1L
  cat(sprintf("  [%s] %-56s %s\n", if (good) "ok" else "FAIL", lab, extra))
}
arm_of <- function(a) cfg$arms[cfg$arms$arm == a, ]

cat("\n== reproducibility: the property the HPC array depends on ==\n")
a <- generate_replicate(cfg, arm_of(1), d = 2, rep = 1, scenario_id = "A01", n = 60)
b <- generate_replicate(cfg, arm_of(1), d = 2, rep = 1, scenario_id = "A01", n = 60)
c3 <- generate_replicate(cfg, arm_of(1), d = 2, rep = 2, scenario_id = "A01", n = 60)
d4 <- generate_replicate(cfg, arm_of(1), d = 2, rep = 1, scenario_id = "A02", n = 60)
ok("same scenario + rep reproduces exactly", identical(a$Ly, b$Ly))
ok("a different rep gives different data", !identical(a$Ly, c3$Ly))
ok("a different scenario gives different data", !identical(a$Ly, d4$Ly))

cat("\n== observation counts: truncated, not clamped ==\n")
g <- generate_replicate(cfg, arm_of(1), d = 2, rep = 1, scenario_id = "A01", n = 4000)
lam <- gnum(cfg, "obs_poisson_lambda"); flo <- gint(cfg, "obs_min")
ok("no subject below obs_min", min(g$n_obs) >= flo)
# The right test is not "few at the floor" but "the floor has EXACTLY the mass a
# truncated Poisson gives it". Redrawing yields Poisson(lambda) conditioned on
# X >= floor; clamping would instead pile all the removed mass onto the floor.
p_floor <- stats::dpois(flo, lam) / (1 - stats::ppois(flo - 1, lam))
ok("floor carries exactly the truncated-Poisson mass, not a clamp's pile-up",
   abs(mean(g$n_obs == flo) - p_floor) < 0.015,
   sprintf("observed %.2f%% against %.2f%% expected; clamping would give %.2f%%",
           100 * mean(g$n_obs == flo), 100 * p_floor,
           100 * stats::ppois(flo, lam)))
x <- flo:60; mu_trunc <- sum(x * stats::dpois(x, lam)) / (1 - stats::ppois(flo - 1, lam))
ok("mean count matches the truncated Poisson",
   abs(mean(g$n_obs) - mu_trunc) < 0.15,
   sprintf("observed %.2f against %.2f expected (nominal lambda = %.0f)",
           mean(g$n_obs), mu_trunc, lam))
ok("floor is honoured and the ceiling is gone",
   min(g$n_obs) == gint(cfg, "obs_min") && max(g$n_obs) > 12,
   sprintf("observed %d-%d, mean %.2f. The ORIGINAL filtered to 4-12; the revision keeps the floor and drops the ceiling",
           min(g$n_obs), max(g$n_obs), mean(g$n_obs)))

cat("\n== score arms hit their target moments ==\n")
mom <- function(arm) {
  gg <- generate_replicate(cfg, arm_of(arm), d = 0, rep = 7,
                           scenario_id = paste0("M", arm), n = 20000)
  x <- gg$xi[[1]][, 1]; x <- (x - mean(x)) / stats::sd(x)
  c(skew = mean(x^3), exkurt = mean(x^4) - 3)
}
for (arm in 1:3) {
  r  <- arm_of(arm); m <- mom(arm)
  tgt <- c(as.numeric(r$skew), as.numeric(r$ex_kurt))
  ok(sprintf("arm %d (%s)", arm, r$label),
     abs(m[1] - tgt[1]) < 0.15 && abs(m[2] - tgt[2]) < 0.60,
     sprintf("skew %+.2f (target %+.1f), ex.kurt %+.2f (target %+.1f)",
             m[1], tgt[1], m[2], tgt[2]))
}
m4 <- mom(4)
ok("arm 4 (t nu=1.5) has an explosive sample kurtosis",
   is.finite(m4[2]) && m4[2] > 50,
   sprintf("ex.kurt %.0f -- the population value is infinite, so a large and unstable sample value is CORRECT",
           m4[2]))

cat("\n== IQR matching: tails differ, spread does not ==\n")
iqr_of <- function(arm) {
  gg <- generate_replicate(cfg, arm_of(arm), d = 0, rep = 9,
                           scenario_id = paste0("Q", arm), n = 20000)
  stats::IQR(gg$xi[[1]][, 1])
}
ref <- iqr_of(1)
for (arm in 2:4) {
  ok(sprintf("arm %d score IQR matches the Gaussian arm", arm),
     abs(iqr_of(arm) / ref - 1) < 0.06,
     sprintf("ratio %.3f", iqr_of(arm) / ref))
}

cat("\n== the t3 error arm ==\n")
ge <- generate_replicate(cfg, arm_of(5), d = 2, rep = 3, scenario_id = "A05", n = 3000)
gg <- generate_replicate(cfg, arm_of(1), d = 2, rep = 3, scenario_id = "A05", n = 3000)
res_e <- unlist(ge$Ly) - ge$X_true[cbind(rep(seq_along(ge$Lt), ge$n_obs),
                 unlist(lapply(ge$Lt, function(t) match(t, ge$tgrid))))]
res_g <- unlist(gg$Ly) - gg$X_true[cbind(rep(seq_along(gg$Lt), gg$n_obs),
                 unlist(lapply(gg$Lt, function(t) match(t, gg$tgrid))))]
ok("error IQR matches the Gaussian arm",
   abs(stats::IQR(res_e) / stats::IQR(res_g) - 1) < 0.08,
   sprintf("ratio %.3f", stats::IQR(res_e) / stats::IQR(res_g)))
ok("but the tails are much heavier",
   mean(abs(res_e) > 3 * stats::IQR(res_g)) > 10 * mean(abs(res_g) > 3 * stats::IQR(res_g)),
   sprintf("P(|e| > 3 IQR): t3 %.4f vs Gaussian %.4f",
           mean(abs(res_e) > 3 * stats::IQR(res_g)),
           mean(abs(res_g) > 3 * stats::IQR(res_g))))

cat("\n== curve contamination ==\n")
gc <- generate_replicate(cfg, arm_of(6), d = 2, rep = 4, scenario_id = "A06", n = 3000)
ok("about 5 percent of subjects are contaminated",
   abs(mean(gc$contaminated) - gnum(cfg, "contam_frac")) < 0.012,
   sprintf("%.1f%%", 100 * mean(gc$contaminated)))
ok("X_true EXCLUDES the contaminant, by design",
   {
     Phi <- lapply(1:2, function(k) eval_eigen(cfg, k, 2L, gc$tgrid))
     mu  <- lapply(1:2, function(k) eval_mean(cfg, k, gc$tgrid, 2))
     i <- which(gc$contaminated)[1]; k <- gc$labels[i]
     max(abs(gc$X_true[i, ] - (mu[[k]] + as.vector(Phi[[k]] %*% gc$xi[[k]][i, ])))) < 1e-9
   },
   "so a robust method is rewarded for resisting the outlier, not punished")
ok("contaminated subjects sit far from their class mean",
   {
     dc <- mean(abs(gc$X_true[gc$contaminated, ] -
                    gc$X_true[!gc$contaminated, ][seq_len(sum(gc$contaminated)), ]))
     is.finite(dc)
   }, "")

cat("\n== class-dependent windows ==\n")
gw <- generate_replicate(cfg, arm_of(1), d = 2, rep = 5, scenario_id = "B07",
                         n = 3000, window = "class")
last <- vapply(gw$Lt, max, numeric(1))
ok("class 2 is followed for a shorter time",
   mean(last[gw$labels == 2]) < mean(last[gw$labels == 1]) - 0.15,
   sprintf("mean last observation: class 1 %.3f, class 2 %.3f",
           mean(last[gw$labels == 1]), mean(last[gw$labels == 2])))
ok("no observation exceeds its subject's window",
   max(last) <= gnum(cfg, "t_max") + 1e-12)
gwc <- generate_replicate(cfg, arm_of(1), d = 2, rep = 5, scenario_id = "B07", n = 500)
ok("the common window spans the whole domain",
   max(vapply(gwc$Lt, max, numeric(1))) > 0.98)

# The non-informative level is the CONTROL for the class-dependent one: windows
# vary, but carry no information about class. If accuracy holds here and drops
# under "class", the culprit is informativeness rather than uneven follow-up.
gn <- generate_replicate(cfg, arm_of(1), d = 2, rep = 5, scenario_id = "B10",
                         n = 4000, window = "noninf")
lastn <- vapply(gn$Lt, max, numeric(1))
ok("non-informative window: follow-up varies",
   stats::sd(lastn) > 0.1,
   sprintf("sd of last observation %.3f", stats::sd(lastn)))
ok("non-informative window: but carries NO information about class",
   abs(mean(lastn[gn$labels == 1]) - mean(lastn[gn$labels == 2])) < 0.03,
   sprintf("class means %.3f vs %.3f -- against %.3f vs %.3f under the class-dependent level",
           mean(lastn[gn$labels == 1]), mean(lastn[gn$labels == 2]),
           mean(last[gw$labels == 1]), mean(last[gw$labels == 2])))
# Fixed-count windows must NOT change dispersion. An earlier version of this
# test asserted the opposite and passed by luck at 0.77 against 0.76; the claim
# was wrong. Truncation alone pushes var/mean BELOW 1, and only tying the count
# to the window can push it above.
ok("fixed-count windows leave dispersion unchanged",
   abs(stats::var(gn$n_obs) / mean(gn$n_obs) -
       stats::var(g$n_obs) / mean(g$n_obs)) < 0.06,
   sprintf("var/mean %.2f under windows against %.2f under a common window",
           stats::var(gn$n_obs) / mean(gn$n_obs),
           stats::var(g$n_obs) / mean(g$n_obs)))

# The three-row ladder: each step adds exactly one thing.
gr <- generate_replicate(cfg, arm_of(1), d = 2, rep = 5, scenario_id = "B11",
                         n = 4000, window = "noninf", count_mode = "rate",
                         obs_min_override = 2)
Tr <- vapply(gr$Lt, max, numeric(1))
ok("rate mode ties count to follow-up length",
   stats::cor(gr$n_obs, Tr) > 0.30,
   sprintf("corr(count, last obs) %+.3f against %+.3f under fixed count",
           stats::cor(gr$n_obs, Tr), stats::cor(gn$n_obs, lastn)))
ok("rate mode gives OVER-dispersed counts, as real cohorts do",
   stats::var(gr$n_obs) / mean(gr$n_obs) > 1,
   sprintf("var/mean %.2f under rate against %.2f under fixed",
           stats::var(gr$n_obs) / mean(gr$n_obs),
           stats::var(gn$n_obs) / mean(gn$n_obs)))
ok("the relaxed floor lets short-follow-up subjects stay sparse",
   min(gr$n_obs) == 2L,
   "a floor of 4 would force them back up, reproducing the selection of the cohort")

cat("\n== contrasts ==\n")
gx <- generate_replicate(cfg, arm_of(1), d = 3, rep = 6, scenario_id = "B02",
                         n = 4000, contrast = "crossing")
lv <- tapply(rowMeans(gx$X_true), gx$labels, mean)
ok("crossing means: class average levels are equal",
   abs(diff(lv)) < 0.25, sprintf("difference %.3f", abs(diff(lv))))
gi <- generate_replicate(cfg, arm_of(1), d = 5, rep = 6, scenario_id = "B03",
                         n = 4000, contrast = "identical")
lvi <- tapply(rowMeans(gi$X_true), gi$labels, mean)
ok("mean-identical: level difference vanishes even at d = 5",
   abs(diff(lvi)) < 0.25, sprintf("difference %.3f", abs(diff(lvi))))
gp <- generate_replicate(cfg, arm_of(1), d = 5, rep = 6, scenario_id = "B01", n = 4000)
lvp <- tapply(rowMeans(gp$X_true), gp$labels, mean)
ok("parallel: level difference is about d",
   abs(abs(diff(lvp)) - 5) < 0.5, sprintf("difference %.3f against d = 5", abs(diff(lvp))))

cat("\n== M_true = 4 ==\n")
g4 <- generate_replicate(cfg, arm_of(1), d = 2, rep = 8, scenario_id = "B06",
                         n = 500, M_true = 4L)
ok("four score columns per class", ncol(g4$xi[[1]]) == 4L)
ok("score variances follow the configured spectrum",
   {
     v <- apply(g4$xi[[1]], 2, stats::var)
     all(diff(v) < 0) && abs(v[1] / v[4] - 9) < 4
   },
   sprintf("variances %s against Z = 9, 4, 2, 1",
           paste(sprintf("%.1f", apply(g4$xi[[1]], 2, stats::var)), collapse = ", ")))

cat(sprintf("\n%s  %d check%s failed\n\n",
            if (fail == 0L) "GENERATOR MATCHES THE DESIGN." else "GENERATOR DOES NOT MATCH THE DESIGN.",
            fail, if (fail == 1L) "" else "s"))
quit(status = if (fail == 0L) 0L else 1L)
