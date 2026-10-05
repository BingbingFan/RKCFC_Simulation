# STUDY F3 -- COX MODELS.
# 36-MONTH VERSION (2026-10-02) of run_realdata_f3.R: the same models, with the
# window, landmark and results folder from realdata36_settings.R.
# SINCE 2026-10-03 FITTED UNDER BOTH TIME ORIGINS: time from
# SURGERY on all patients, the submitted design, and the month-36 LANDMARK. Which
# one the tables lead with is REAL_COX_ORIGIN in realdata36_settings.R ("surgery");
# the other is the sensitivity analysis. Both are printed and saved in full.
#
#   Rscript R/run_realdata36_f3.R [f1dir] [outdir]
#
# THE DEFECT. The submitted text took the hazard ratio from Model 3 and the P
# value from Model 1. That is not a reporting slip: it makes the adjusted effect
# look as though it carried the unadjusted model's significance. Every quantity
# below is therefore printed FOR EACH MODEL AS A COMPLETE SET -- HR, confidence
# interval and P together, from the same fit -- so no combination across models
# is possible without deliberately reaching into two rows.
#
# ON THE TWO TIME ORIGINS. Group membership is built from CEA over months 0-36.
# From surgery, the exposure is not known at time zero and a patient had to live
# long enough to be classified (guarantee-time bias); the landmark at
# month 36 removes that by starting the clock when membership is known, at the
# cost of the 146 patients who died or were censored before it and of the deaths
# inside the window. The from-surgery estimate is reported, with the landmark
# one beside it.
#
# ON CEA AS A COVARIATE.
#
# A "postoperative CEA" covariate may be over-adjustment. The
# extract carries no postoperative CEA column, so that exact model cannot be
# reproduced. But **CEA at time 0 IS the
# pre-surgical measurement**, and checking it against the baseline table settles
# the matter: `Preoperative_CEA` is IDENTICAL to the t = 0 point of the
# trajectory series in 1403 of 1405 subjects, correlation 0.999997.
#
# So `Preoperative_CEA` is not a safe baseline covariate. It is a LITERAL
# COMPONENT OF THE EXPOSURE -- one of the measurements the clustering used to
# build the groups, present for every subject and making up 13 percent of all
# observations. Adjusting the group effect for it conditions on part of the
# thing being estimated, which attenuates toward the null for reasons that have
# nothing to do with biology.
#
# THEREFORE: the primary adjusted model EXCLUDES it, and a separate model
# includes it, labelled as the over-adjustment sensitivity it is. That is the
# treatment over-adjustment calls for, applied to the variable that actually has
# the problem.
#
# MISSING DATA. Stated rather than assumed: the counts are printed per covariate
# and the analysis is complete-case, which is defensible only because the
# covariates are in fact complete. If a future extract is not, this prints it.

source(file.path("R", "methods.R"))
source(file.path("R", "realdata36_settings.R"))

args   <- commandArgs(trailingOnly = TRUE)
args   <- args[args != "--args"]
f1dir  <- if (length(args) >= 1) args[1] else file.path(REAL_RESULTS, "F1")
outdir <- if (length(args) >= 2) args[2] else file.path(REAL_RESULTS, "F3")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("survival", quietly = TRUE))
  stop("package `survival` is required for F3.", call. = FALSE)

LANDMARK <- REAL_LANDMARK
ORIGINS  <- c(REAL_COX_ORIGIN, setdiff(c("surgery", "landmark"), REAL_COX_ORIGIN))   # primary first
f1 <- readRDS(file.path(f1dir, "F1.rds"))
# the clustering window stays [0, Tmax] whatever the time origin
d  <- truncate_cea36(load_cea36(), REAL_TMAX); b <- d$base
pk <- peak_cea(d)
prim <- if (exists("REAL_PRIMARY")) REAL_PRIMARY else "raw_M1"
g <- f1$fits[[prim]]
if (is.null(g)) stop("F1 has no `", prim, "` partition.", call. = FALSE)
if (length(g) != nrow(b)) stop("partition length does not match the cohort", call. = FALSE)

b$group <- factor(g)
# PEAK CEA, the comparator F0C made necessary: a single summary reproduces a
# level-split partition at AUC 0.95, so a model set that never puts the two side
# by side leaves the obvious question unanswered. M5 below is that model.
b$peak_cea <- pk

# THE ELEVEN COVARIATES OF THE SUBMITTED ANALYSIS PLUS BMI.
# Tumor_differentiation
# carries an "Unknown" level, a missing-indicator category rather than a missing
# value; every other covariate is complete, and the complete-case count printed
# below shows it.
COV12 <- c("Age", "Sex", "BMI", "Surgical_approach", "Primary_site", "Tumor_differentiation",
           "AJCC_8th_ed_Stage", "Lymph_node_yield", "Mucinous_colloid_type",
           "Lymphovascular_invasion", "Perineural_invasion", "Adjuvant_chemotherapy")
COVS <- list(
  M1 = character(0),
  M2 = c("Age", "Sex"),
  # PRIMARY. No CEA term: the t = 0 CEA is part of the exposure.
  M3 = COV12,
  # SENSITIVITY, and reported as one. Adds the preoperative CEA, which is the
  # trajectory's own t = 0 point. Its purpose is to SHOW the attenuation, not to
  # provide an estimate anyone should quote.
  M4 = c(COV12, "Preoperative_CEA"),
  # THE COMPARATOR, not a confounder adjustment. Peak CEA over the trajectory
  # window reproduces much of the partition, so this model asks the question a
  # reader will: once peak CEA is in the model, does the GROUP term still carry
  # anything? Unlike M4's Preoperative_CEA this is not part of the exposure in the
  # same conditioning sense -- it is a rival summary of the same trajectory, and
  # the comparison is the point rather than a hazard to avoid.
  M5 = c(COV12, "peak_cea"))
PRIMARY_MODEL <- "M3"

# the analysis set under one time origin
cohort_of <- function(origin) {
  if (origin == "landmark") {
    bb <- b[b$Survival_Time >= LANDMARK, , drop = FALSE]
    bb$.t <- bb$Survival_Time - LANDMARK
  } else {
    bb <- b
    bb$.t <- bb$Survival_Time
  }
  bb$.e <- bb$Death
  bb
}
fit_models <- function(bb, origin) {
  fits <- list(); rows <- list()
  for (mn in names(COVS)) {
    rhs <- c("group", COVS[[mn]])
    fml <- stats::as.formula(paste("survival::Surv(.t, .e) ~", paste(rhs, collapse = " + ")))
    fit <- survival::coxph(fml, data = bb)
    s   <- summary(fit)
    fits[[mn]] <- fit
    # THE GROUP EFFECT, AS ONE COMPLETE SET FROM THIS MODEL. HR, interval and P
    # come from the same `summary()` object and are written on one row.
    i <- grep("^group", rownames(s$coefficients))[1]
    rows[[mn]] <- data.frame(
      origin = origin, model = mn,
      adjusted_for = if (length(COVS[[mn]])) paste(COVS[[mn]], collapse = ", ") else "nothing",
      n = s$n, events = s$nevent,
      hr = unname(s$conf.int[i, 1]), lo = unname(s$conf.int[i, 3]),
      hi = unname(s$conf.int[i, 4]), p = unname(s$coefficients[i, 5]),
      stringsAsFactors = FALSE)
  }
  list(fits = fits, table = do.call(rbind, rows))
}

cat(sprintf("\n== F3: Cox models, partition %s; reported origin: %s ==\n", prim, REAL_COX_ORIGIN))
allc <- unique(unlist(COVS))
cat("\nMISSING DATA, per covariate (analysis is complete-case; all patients)\n")
for (v in allc)
  cat(sprintf("  %-24s %d missing of %d\n", v, sum(is.na(b[[v]])), nrow(b)))
cat(sprintf("  complete cases across all covariates: %d of %d\n", sum(stats::complete.cases(b[, allc, drop = FALSE])), nrow(b)))

res <- list(); tab <- list(); ph_tab <- list()
for (origin in ORIGINS) {
  bb <- cohort_of(origin)
  r  <- fit_models(bb, origin)
  res[[origin]] <- r; tab[[origin]] <- r$table
  cat(sprintf("\n%s TIME FROM %s: %d patients, %d events%s\n",
              if (origin == REAL_COX_ORIGIN) "REPORTED --" else "SENSITIVITY --",
              if (origin == "surgery") "SURGERY (all patients)" else sprintf("THE MONTH-%d LANDMARK", LANDMARK),
              nrow(bb), sum(bb$.e),
              if (origin == "landmark") sprintf(" (%d of %d died or were censored before the landmark)", nrow(b) - nrow(bb), nrow(b)) else ""))
  cat(sprintf("  %-6s %-52s %5s %6s %20s %10s\n", "model", "adjusted for", "n", "events", "HR (95% CI)", "p"))
  t1 <- r$table
  for (i in seq_len(nrow(t1)))
    cat(sprintf("  %-6s %-52s %5d %6d %20s %10.4g\n", t1$model[i], substr(t1$adjusted_for[i], 1, 52), t1$n[i], t1$events[i],
                sprintf("%.2f (%.2f-%.2f)", t1$hr[i], t1$lo[i], t1$hi[i]), t1$p[i]))
  # the group's own contribution once peak CEA is in the model (M5 against M5 without the group)
  f5r <- stats::as.formula(paste("survival::Surv(.t, .e) ~", paste(COVS[["M5"]], collapse = " + ")))
  lr5 <- 2 * (as.numeric(stats::logLik(r$fits[["M5"]])) - as.numeric(stats::logLik(survival::coxph(f5r, data = bb))))
  i3 <- which(t1$model == "M3"); i5 <- which(t1$model == "M5"); i4 <- which(t1$model == "M4")
  cat(sprintf("  does the group survive peak CEA?  HR %.2f -> %.2f, p %.4g -> %.4g; group's LR chisq in M5 %.2f on 1 df, p = %.4g\n",
              t1$hr[i3], t1$hr[i5], t1$p[i3], t1$p[i5], lr5, stats::pchisq(lr5, df = 1, lower.tail = FALSE)))
  l3 <- log(t1$hr[i3]); l4 <- log(t1$hr[i4])
  if (is.finite(l3) && abs(l3) > 0.01)
    cat(sprintf("  attenuation from conditioning on the preoperative CEA (M4): HR %.2f -> %.2f (%.0f%% of the log hazard ratio removed)\n",
                t1$hr[i3], t1$hr[i4], 100 * (1 - l4 / l3)))
  # proportional hazards, scaled Schoenfeld residuals, primary model
  zp <- survival::cox.zph(r$fits[[PRIMARY_MODEL]]); tz <- zp$table; ph_tab[[origin]] <- tz
  cat(sprintf("  proportional hazards (%s): group p %.3f; global p %.3f%s\n", PRIMARY_MODEL,
              tz["group", "p"], tz["GLOBAL", "p"],
              if (any(tz[rownames(tz) != "GLOBAL", "p"] < 0.05)) sprintf(" (violations: %s)", paste(rownames(tz)[rownames(tz) != "GLOBAL" & tz[, "p"] < 0.05], collapse = ", ")) else ""))
  res[[origin]]$lr5 <- lr5
}
tab <- do.call(rbind, tab); rownames(tab) <- NULL

cat("\n  PRIMARY IS M3 under the reported origin. M1 and M2 make the adjustment path\n")
cat("  visible; they are not a menu. M4 is the OVER-ADJUSTMENT SENSITIVITY and must\n")
cat("  never be quoted as an estimate: it adds Preoperative_CEA, which IS the\n")
cat("  trajectory's t = 0 point. M5 IS THE COMPARATOR: if the group term survives\n")
cat("  it, the clustering carries survival information a threshold on peak CEA\n")
cat("  does not; if not, say so.\n")
cat("\n  ONE THING FOR THE METHODS SECTION. Every subject's trajectory begins with a\n")
cat("  PRE-treatment value, so the groups partly encode baseline tumour burden\n")
cat("  rather than post-operative dynamics alone. That is a property of the design\n")
cat("  as submitted, not a defect introduced here, and it should be stated.\n")

saveRDS(list(landmark = LANDMARK, partition = prim, origin = REAL_COX_ORIGIN, origins = ORIGINS,
             table = tab, ph = ph_tab[[REAL_COX_ORIGIN]], ph_by_origin = ph_tab,
             lr5_by_origin = lapply(res, `[[`, "lr5"),
             n_at_risk = nrow(cohort_of(REAL_COX_ORIGIN)), n_events = sum(cohort_of(REAL_COX_ORIGIN)$.e),
             primary_model = PRIMARY_MODEL,
             code_hash = source_hash(c(file.path("R", "run_realdata36_f3.R"), REAL36_SETTINGS_FILES))),
        file.path(outdir, "F3.rds"))
cat(sprintf("\nwrote %s\n", file.path(outdir, "F3.rds")))
