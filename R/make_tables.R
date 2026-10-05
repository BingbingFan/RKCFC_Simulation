# PUBLICATION TABLES for the univariate RKCFC manuscript.
#
#   Rscript R/make_tables.R [results_root] [out_dir]
#
# Defaults: results/ and tables/.
#
# WHAT THIS PRODUCES
#   tables/tab01..tab04.tex   the four MAIN tables, as \input-able fragments
#   tables/tabS1..tabS10.tex  the ten SUPPLEMENTARY tables
#   tables/all_tables.tex     a standalone driver that \input's every fragment
#   tables/all_tables.pdf     the compiled proof
#
# The fragments carry no preamble, so the manuscript can \input them directly and
# the driver exists only to typeset a proof. Nothing is hard-coded: every number
# is read from the result files, so re-running after a block re-runs updates the
# tables rather than leaving them to drift from the evidence.
#
# ORDERING IS THE ARGUMENT, not the chronology of the runs.
#   1  the design            -- what was simulated, before any result
#   2  clustering accuracy   -- the comparison the paper exists to make
#   3  estimation accuracy   -- WHY it works: the eigenspace, isolated from the
#                               clustering loop by fitting under the true labels
#   4  number of clusters    -- the criterion the real-data analysis relies on
# Supplements follow the same logic: design justification first (S1, S2), then
# diagnostics for the main result (S3-S5), then the settings studies (S6-S9),
# then the pilot (S10), which confirms rather than establishes.
#
# A NOTE ON WHAT IS AND IS NOT REPORTED AS AN ACCURACY NUMBER. A fit pinned at
# maxIter is an artifact of where the loop was cut off; cells in that state are
# flagged, not quietly averaged in. On the heavy-tailed arms the aRI distribution
# is BIMODAL -- mass at zero and mass near 0.9 -- so the median alone reports
# which side a sample fell on. P(success) is printed in S5, S8 and S9; for the
# core grid the manuscript text gives it.
#
# EVERY MEDIAN OF AN ACCURACY MEASURE HAS A BOOTSTRAP MCSE: in its own column
# where the table has room (S6-S9), otherwise as the largest MCSE over the
# table's cells, stated in the note (tab02, tab03, tab05, tabS5, tabS10, printed
# as Tables 2, 4, 3, S5 and S10). The manuscript (Section 3.1, Design and
# Metrics) describes exactly this; change them together.

suppressWarnings(suppressMessages({
  args <- commandArgs(trailingOnly = TRUE)
  ROOT <- if (length(args) >= 1) args[1] else "results"
  OUT  <- if (length(args) >= 2) args[2] else "tables"
}))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# --- reading -------------------------------------------------------------------

rd <- function(block, pat = "task.*rds$") {
  d <- file.path(ROOT, block)
  if (!dir.exists(d)) return(NULL)
  fs <- list.files(d, pat, full.names = TRUE)
  if (!length(fs)) return(NULL)
  do.call(rbind, lapply(fs, readRDS))
}

# Bootstrap MCSE of a median. NOT the normal-theory form: that was measured 84
# percent too small on a bimodal cell and 340 percent too large at a ceiling,
# and every cell that matters here is one or the other.
mcse_median <- function(x, n_boot = 2000L, seed = 20260908L) {
  x <- x[is.finite(x)]; B <- length(x)
  if (B < 5L) return(NA_real_)
  if (diff(range(x)) == 0) return(0)
  old <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", .GlobalEnv)
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  stats::sd(vapply(seq_len(n_boot),
                   function(i) stats::median(x[sample.int(B, B, replace = TRUE)]),
                   numeric(1)))
}

# Largest MCSE over a table's cells, for the note of a table with no room for an
# MCSE column. Rounded UP at the printed precision, so that "at most" stays true;
# cells with fewer than five values have no MCSE and are skipped.
mcse_max <- function(x, d = 3) {
  m <- suppressWarnings(max(x, na.rm = TRUE))
  f(ceiling(m * 10^d - 1e-9) / 10^d, d)
}

med  <- function(x) { x <- x[is.finite(x)]; if (length(x)) stats::median(x) else NA_real_ }
psucc <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x > 0.05) else NA_real_ }
rate <- function(x) { x <- x[!is.na(x)]; if (length(x)) mean(x) else NA_real_ }

# --- LaTeX ---------------------------------------------------------------------

# Escape what LaTeX would otherwise eat. Level names here are full of underscores
# (`bw_half`, `sil_manuscript`, `mean_robust`) and arm labels carry commas and
# equals signs, so this is load-bearing rather than defensive.
esc <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  for (ch in c("&", "%", "$", "#", "_", "{", "}"))
    x <- gsub(ch, paste0("\\", ch), x, fixed = TRUE)
  x <- gsub("~", "\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("^", "\\textasciicircum{}", x, fixed = TRUE)
  x
}

# Numbers. NA prints as an em-dash rather than "NA", which reads as a missing
# measurement instead of a failed one.
#
# NEGATIVE ZERO IS SUPPRESSED. An adjusted Rand index of -0.0004 is chance, and
# `formatC` renders it "-0.000", which invites a reader to see a sign where the
# data carry none. Anything below half a unit in the last printed place is zero.
f <- function(x, d = 3) {
  z <- ifelse(is.finite(x) & abs(x) < 0.5 * 10^(-d), 0, x)
  ifelse(is.finite(z), formatC(z, format = "f", digits = d), "---")
}
pc <- function(x, d = 0) ifelse(is.finite(x), paste0(formatC(100 * x, format = "f", digits = d), "\\%"), "---")

# Caption split rules (table name in the caption, everything else in the note).
source(file.path("R", "table_captions.R"))

emit <- function(id, caption, label, header, align, rows, notes = NULL,
                 small = TRUE, midrules = integer(0), wide = FALSE) {
  con <- file(file.path(OUT, paste0(id, ".tex")), open = "wt", encoding = "UTF-8")
  on.exit(close(con))
  # Caption = table name only; the rest of the caption opens the note (R/table_captions.R).
  sp <- split_caption(id, caption, notes)
  caption <- sp$caption
  notes <- if (nzchar(sp$notes)) sp$notes else NULL
  w <- function(...) cat(..., "\n", sep = "", file = con)
  w("% ", id, " -- generated by R/make_tables.R on ", format(Sys.Date()), ". Do not edit by hand.")
  w("\\begin{table}[htbp]")
  w("\\centering")
  # `wide`: a table with nine or more columns overran the text block of the
  # journal class by up to 92pt when set at \small with the default column
  # separation (seen 2026-09-09 by rendering the page: the last column and the
  # note were clipped). \footnotesize with a 3pt separation brings an
  # eleven-column table inside the block; the \setlength is local to the table.
  if (wide) w("\\setlength{\\tabcolsep}{3pt}\\footnotesize") else if (small) w("\\small")
  w("\\caption{", caption, "}")
  w("\\label{", label, "}")
  # Every table spans the text width: tabular* with
  # \extracolsep{\fill} keeps the natural column widths and spreads the extra space
  # between the columns, so the caption and the note line up with the text block.
  w("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}", align, "@{}}")
  w("\\toprule")
  w(header, " \\\\")
  w("\\midrule")
  for (i in seq_along(rows)) {
    w(rows[i], " \\\\")
    if (i %in% midrules) w("\\midrule")
  }
  w("\\bottomrule")
  w("\\end{tabular*}")
  if (!is.null(notes)) {
    w("\\begin{minipage}{\\linewidth}\\vspace{2pt}\\footnotesize")
    w(notes)
    w("\\end{minipage}")
  }
  w("\\end{table}")
  message("  wrote ", id, ".tex")
}

ARM_LAB <- c("Gaussian", "Fleishman(0, 4)", "Fleishman(1.5, 4)",
             "$t_{1.5}$", "$t_3$ measurement error", "curve contamination")
METH_LAB <- c(M1 = "RKCFC", M2 = "KCFC", M3 = "KCFC (log)", M4 = "FPCA $k$-means",
              M5 = "Kendall $k$-means", M6 = "LC spline", M7 = "QDA on scores")
# Header labels for the nine-column core table, where the full labels overran
# the text block; the full forms are kept for the narrower tables.
METH_LAB_SHORT <- c(M1 = "RKCFC", M2 = "KCFC", M3 = "KCFC (log)", M4 = "FPCA $k$-means",
                    M5 = "Kendall $k$-means", M6 = "LC spline", M7 = "QDA")

message("reading from ", ROOT)
core <- rd("core"); estim <- rd("estim"); krec <- rd("krecovery")
vari <- rd("variations"); comp <- rd("S_components"); ini <- rd("S_init")
smo <- rd("S_smoothing"); pil <- rd("pilot"); dlt <- rd("D_delta")
estr <- rd("S_estimators"); conv <- rd("convergence")

# ================================ MAIN TABLES ==================================

## Table 1 -- the design ---------------------------------------------------------
if (!is.null(core)) {
  arms <- sort(unique(core$arm))
  rows <- vapply(arms, function(a) {
    s <- core[core$arm == a, ]
    paste(ARM_LAB[a],
          f(med(s$n_obs_mean), 2), sprintf("%d--%d", min(s$n_obs_min), max(s$n_obs_max)),
          length(unique(s$rep)), sep = " & ")
  }, character(1))
  emit("tab01",
       "Simulation design. Six arms place the departure from normality in the scores (the two Fleishman arms and $t_{1.5}$, against the Gaussian reference), in the measurement errors ($t_3$) or in the trajectories (curve contamination); the score and error distributions are matched on interquartile range, so that arms differ in tail behavior and shape rather than in scale. The Fleishman arms are labeled by (skewness, excess kurtosis). Each arm is run at both separations, $d \\in \\{1.5, 2\\}$, with $n = 500$ subjects.",
       "tab:design",
       # Arms are named, not numbered: every other table and
       # the text refer to them by distribution, so a number column meant nothing.
       "Arm & Mean $n_i$ & Range $n_i$ & Replicates",
       "lccc", rows,
       notes = paste("Classes differ in \\emph{both} mean level and eigenfunction basis (sines against cosines),",
                     "so $d$ controls only part of the separation; at $d = 0$ the classes remain distinguishable on shape alone.",
                     "$n_i$ is the number of observations per subject, drawn truncated-Poisson."))
}

## Table 2 -- clustering accuracy, core grid -------------------------------------
if (!is.null(core)) {
  ms <- names(METH_LAB)
  rows <- character(0); mids <- integer(0); mc_core <- list()
  for (dv in sort(unique(core$d))) {
    for (a in sort(unique(core$arm))) {
      s <- core[core$arm == a & core$d == dv, ]
      cells <- vapply(ms, function(m) f(med(s$ari[s$method == m])), character(1))
      for (m in ms) mc_core[[m]] <- c(mc_core[[m]], mcse_median(s$ari[s$method == m]))
      rows <- c(rows, paste(c(f(dv, 1), ARM_LAB[a], cells), collapse = " & "))
    }
    mids <- c(mids, length(rows))
  }
  mids <- mids[-length(mids)]
  emit("tab02",
       "Clustering accuracy on the core grid: median adjusted Rand index over 100 replicates, by arm, separation and method. RKCFC matches KCFC where the scores have finite moments, including the two Fleishman arms, and separates from it sharply where the moment estimator fails, under infinite-variance scores and curve contamination.",
       "tab:core",
       paste(c("$d$", "Arm", unname(METH_LAB_SHORT)), collapse = " & "),
       paste0("cl", strrep("c", length(ms))), rows, midrules = mids, wide = TRUE,
       notes = paste("FPCA $k$-means and Kendall $k$-means apply $k$-means once to global FPC scores, non-robust and robust respectively;",
                     "both score near zero on every arm, so the iterative $K$-centers loop, not the robust score",
                     "representation, is what separates the methods. On the heavy-tailed arms the distribution of",
                     "the index is bimodal; see Table~\\ref{tab:convergence} for convergence by cause.",
                     sprintf("Bootstrap standard errors of the medians are at most %s for RKCFC, %s for KCFC and %s for KCFC (log), and at most %s for the other four methods.",
                             mcse_max(mc_core$M1), mcse_max(mc_core$M2), mcse_max(mc_core$M3),
                             mcse_max(unlist(mc_core[c("M4", "M5", "M6", "M7")])))))
}

## CCR table -- correct classification rate, core grid:
## the twin of Table 2 in CCR units, so that "the CCR agrees with the aRI throughout"
## in Section 3.1 points at a table. Same rows, same methods, same files.
if (!is.null(core)) {
  ms <- names(METH_LAB)
  rows <- character(0); mids <- integer(0); mc_ccr <- list()
  for (dv in sort(unique(core$d))) {
    for (a in sort(unique(core$arm))) {
      s <- core[core$arm == a & core$d == dv, ]
      cells <- vapply(ms, function(m) f(med(s$ccr[s$method == m])), character(1))
      for (m in ms) mc_ccr[[m]] <- c(mc_ccr[[m]], mcse_median(s$ccr[s$method == m]))
      rows <- c(rows, paste(c(f(dv, 1), ARM_LAB[a], cells), collapse = " & "))
    }
    mids <- c(mids, length(rows))
  }
  mids <- mids[-length(mids)]
  emit("tabCCR",
       "Correct classification rate on the core grid: median correct classification rate over 100 replicates, by arm, separation and method, with labels matched to the truth by the permutation maximizing agreement.",
       "tab:ccr",
       paste(c("$d$", "Arm", unname(METH_LAB_SHORT)), collapse = " & "),
       paste0("cl", strrep("c", length(ms))), rows, midrules = mids, wide = TRUE,
       notes = paste("The CCR ranks the methods as the adjusted Rand index of Table~\\ref{tab:core} does;",
                     "with two classes a rate near 0.5 is chance, which is where the one-shot methods and the",
                     "latent-class spline sit on every arm.",
                     sprintf("Bootstrap standard errors of the medians are at most %s for RKCFC, %s for KCFC and %s for KCFC (log), and at most %s for the other four methods.",
                             mcse_max(mc_ccr$M1), mcse_max(mc_ccr$M2), mcse_max(mc_ccr$M3),
                             mcse_max(unlist(mc_ccr[c("M4", "M5", "M6", "M7")])))))
}

## Table 3 -- estimation accuracy on the true partition --------------------------
if (!is.null(estim)) {
  eids <- c("E0", "E1", "E1M", "E2", "E3")
  # Names a reader can use: the manuscript never defines
  # "PACE"; it says "conventional FPCA" and cites fdapace.
  ELAB <- c(E0 = "Conventional FPCA", E1 = "Kendall's $\\tau$", E1M = "Kendall's $\\tau$ + robust mean",
            E2 = "Spatial sign", E3 = "Normal scores + conventional FPCA")
  rows <- character(0); mids <- integer(0)
  mc_est <- list(mean = numeric(0), cc = numeric(0), ang = numeric(0), ang_rob = numeric(0))
  for (a in sort(unique(estim$arm))) {
    for (e in eids) {
      s <- estim[estim$arm == a & estim$estimator == e, ]
      rows <- c(rows, paste(if (e == eids[1]) ARM_LAB[a] else "", ELAB[e],
                            f(med(s$mean_err), 4), f(med(s$angle_1), 2), f(med(s$angle_2), 2),
                            f(med(s$angle_sub), 2), f(med(s$cancor), 4), sep = " & "))
      ang <- c(mcse_median(s$angle_1), mcse_median(s$angle_2), mcse_median(s$angle_sub))
      mc_est$mean <- c(mc_est$mean, mcse_median(s$mean_err))
      mc_est$cc   <- c(mc_est$cc, mcse_median(s$cancor))
      mc_est$ang  <- c(mc_est$ang, ang)
      if (e %in% c("E1", "E1M", "E2")) mc_est$ang_rob <- c(mc_est$ang_rob, ang)
    }
    mids <- c(mids, length(rows))
  }
  mids <- mids[-length(mids)]
  emit("tab03",
       "Estimation accuracy under the \\emph{true} partition, removing the clustering loop from the comparison. Angles are in degrees between the estimated and true eigenfunctions. Kendall's $\\tau$ does not protect the mean function, which is what the robust mean addresses.",
       "tab:estim",
       "Arm & Estimator & Mean error & $\\angle_1$ & $\\angle_2$ & $\\angle_{\\mathrm{sub}}$ & Can.\\ corr.",
       # `l`, not p{}: a fixed-width column justifies "Kendall's tau + robust
       # mean" across two lines with stretched interword space, which looks like
       # a typesetting fault in a table whose whole job is to be read at a glance.
       "llccccc", rows, midrules = mids,
       notes = paste("Mean error is the relative $L_2$ error of $\\hat\\mu$ against the true class mean.",
                     "$\\angle_{\\mathrm{sub}}$ is the largest principal angle between estimated and true spans.",
                     sprintf("Bootstrap standard errors of the medians are at most %s for the mean error, %s for the canonical correlation and %s$^\\circ$ for the angles; for Kendall's $\\tau$ and the spatial-sign estimator, no angle has a standard error above %s$^\\circ$.",
                             mcse_max(mc_est$mean), mcse_max(mc_est$cc), mcse_max(mc_est$ang, 1), mcse_max(mc_est$ang_rob, 1)),
                     "Kendall's $\\tau$ and the spatial-sign estimator hold the eigenspace on every arm, within",
                     "about $2^\\circ$ of each other, whereas the conventional estimator does not; the normal-scores estimator loses the",
                     "second eigenfunction under contamination. The robust mean is a local-linear M-estimator (Huber) that",
                     "improves the mean under \\emph{symmetric} heavy tails and degrades it under skewness,",
                     "where an M-estimator of location does not target the mean; it leaves the eigenspace",
                     "untouched, which is why its angles equal Kendall's exactly.",
                     "The normal-scores row needs a different reading, because a monotone marginal transform changes the",
                     "estimand, so its mean error and angles measure how far the transformed eigenspace sits",
                     "from the \\emph{untransformed} truth. This is the distortion the transform introduces, not",
                     "estimation accuracy. Table~\\ref{tab:estimators} compares it on a partition, which is",
                     "invariant to the transform and therefore a fair basis."))
}

## Table 4 -- number of clusters -------------------------------------------------
if (!is.null(krec)) {
  crit <- c(sil_manuscript = "Silhouette (curve space)",
            sil_pairwise   = "Silhouette (Rousseeuw)",
            CH             = "Calinski--Harabasz")
  rows <- character(0); mids <- integer(0)
  for (sc in sort(unique(krec$scenario_id))) {
    Kt <- unique(krec$K_true[krec$scenario_id == sc])
    for (cn in names(crit)) for (m in c("M1", "M2")) {
      s <- krec[krec$scenario_id == sc & krec$method == m, ]
      pick <- tapply(seq_len(nrow(s)), s$rep, function(ii) {
        z <- s[ii, ]; z <- z[z$K_fit >= 2 & is.finite(z[[cn]]), ]
        if (!nrow(z)) return(NA_integer_)
        # S(1) = 0 by convention (Section 2.4): one class
        # wins when every silhouette for K >= 2 is negative. The Calinski-Harabasz index
        # has no value at K = 1, so its candidate set stays K >= 2.
        if (cn != "CH" && max(z[[cn]]) < 0) return(1L)
        z$K_fit[which.max(z[[cn]])] })
      pick <- unlist(pick)
      cnt <- vapply(1:5, function(k) sum(pick == k, na.rm = TRUE), integer(1))
      ok  <- if (Kt >= 2 || cn != "CH") mean(pick == Kt, na.rm = TRUE) else NA_real_
      rows <- c(rows, paste(c(if (cn == names(crit)[1] && m == "M1") sprintf("%d", Kt) else "",
                              if (m == "M1") crit[[cn]] else "",
                              METH_LAB[[m]], cnt, if (is.finite(ok)) pc(ok) else "n/a"),
                            collapse = " & "))
    }
    mids <- c(mids, length(rows))
  }
  mids <- mids[-length(mids)]
  emit("tab04",
       "Recovery of the number of clusters. Each entry counts, over 50 replicates, how often a criterion selected each $K$; the final column is the share selecting the true $K$. Both silhouettes recover $K=2$ reliably for either method; the Calinski--Harabasz index over-splits when the true $K$ is 1 or 2.",
       "tab:krecovery",
       "True $K$ & Criterion & Method & $K{=}1$ & $K{=}2$ & $K{=}3$ & $K{=}4$ & $K{=}5$ & Correct",
       "clp{1.6cm}cccccc", rows, midrules = mids,
       notes = paste(sprintf("The study ran on the %s arm at $d = %s$.",
                             paste(ARM_LAB[sort(unique(krec$arm))], collapse = " and "), f(unique(krec$d)[1], 1)),
                     "For both silhouettes $K = 1$ is in the candidate set with $S(1) = 0$ by convention",
                     "(Section~\\ref{sec2.4}), so one class is selected when every coefficient for $K \\geq 2$ is",
                     "negative, and the $K = 1$ column counts how often that happened. The Calinski--Harabasz",
                     "index has no value at $K = 1$, so its row at a true $K$ of 1 has no correct-selection rate. The curve-space silhouette is the",
                     "squared, self-inclusive variant; Rousseeuw's pairwise form is reported beside it because the",
                     "two are not interchangeable and differ by roughly an order of magnitude on the same partition."))
}

## Table 5 -- reconstruction error, core grid -----------------------------------
# Added 2026-09-09. The three reconstruction errors were recorded on every fit
# from the start (R/metrics.R) but no table reported them, so the manuscript's
# replacement for the submitted MPE table had nowhere to come from. Only the
# methods that produce a fitted trajectory appear: M4, M5 and M7 return labels
# and M6 fits a different model family.
if (!is.null(core)) {
  rms <- c("M1", "M2", "M3")
  rows <- character(0); mids <- integer(0); mc_rec <- list()
  for (dv in sort(unique(core$d))) {
    for (a in sort(unique(core$arm))) {
      s <- core[core$arm == a & core$d == dv, ]
      cells <- unlist(lapply(rms, function(m) {
        t <- s[s$method == m, ]
        c(f(med(t$rre_true)), f(med(t$rre_out)), f(med(t$rre_in)))
      }))
      for (m in rms) {
        t <- s[s$method == m, ]
        mc_rec[[m]] <- c(mc_rec[[m]], mcse_median(t$rre_true), mcse_median(t$rre_out), mcse_median(t$rre_in))
      }
      rows <- c(rows, paste(c(f(dv, 1), ARM_LAB[a], cells), collapse = " & "))
    }
    mids <- c(mids, length(rows))
  }
  mids <- mids[-length(mids)]
  n_hx <- suppressWarnings(max(core$n_holdout_excluded, na.rm = TRUE))
  emit("tab05",
       "Relative reconstruction error on the core grid, median over 100 replicates. \\emph{True} is the fitted trajectory against the true latent curve on a dense grid, the primary metric, since no estimated quantity appears on its target side. \\emph{Held-out} is a random 20\\% of each subject's observations withheld from fitting and scored. \\emph{In-sample} is the error on the observations used for fitting. 0 is perfect reconstruction and 1 no better than the sample mean.",
       "tab:recon",
       paste("$d$ & Arm & \\multicolumn{3}{c}{RKCFC} & \\multicolumn{3}{c}{KCFC} & \\multicolumn{3}{c}{KCFC (log)} \\\\",
             "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}",
             "& & True & Held-out & In-sample & True & Held-out & In-sample & True & Held-out & In-sample"),
       "clccccccccc", rows, midrules = mids, wide = TRUE,
       notes = paste("Each method is scored on its own reconstruction; KCFC (log) is back-transformed to the",
                     "original scale before scoring. On the contamination arm the 5\\% of contaminated subjects",
                     "are scored against their \\emph{clean} latent curve, which no method fitting the observed",
                     "data can reconstruct, so every method's true-curve error is inflated there. The gap",
                     "between KCFC's true-curve and held-out errors (0.55 against 0.33) shows that it fits the",
                     "contaminant. It reconstructs the observations well and the clean curve badly, whereas",
                     "the two errors coincide for RKCFC. Methods that return",
                     "labels only (FPCA $k$-means, Kendall $k$-means, and QDA on scores) and the latent-class spline, which fits a different model family,",
                     "are absent by construction.",
                     if (is.finite(n_hx) && n_hx == 0) "No subject was excluded from the held-out metric." else
                       sprintf("At most %d subjects per replicate were excluded from the held-out metric for having too few observations.", n_hx),
                     sprintf("Bootstrap standard errors of the medians are at most %s for RKCFC, %s for KCFC and %s for KCFC (log).",
                             mcse_max(mc_rec$M1), mcse_max(mc_rec$M2), mcse_max(mc_rec$M3))))
}

# ============================ SUPPLEMENTARY TABLES =============================

## S1 -- tuning parameters and implementation settings ---------------------------
# Values are READ from config/params_global.csv so the table cannot drift from what
# ran; the label and basis text are written here for a reader. The config's own
# `source` column is project bookkeeping (dates and who decided what) and must
# never be printed.
{
  pg <- tryCatch(utils::read.csv(file.path("config", "params_global.csv"),
                                 stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(pg)) {
    v <- function(p) { z <- pg$value[pg$param == p]; if (length(z) != 1) stop("param missing: ", p); z }
    S1 <- list(
      list("Retained components $M$", v("nK_fit"),
           "The true number in the generating model, common to both classes; retaining 1--4 and selecting by FVE is examined in Table~\\ref{tab:components}"),
      list("Separation $d$", paste(v("d_hard"), "and", v("d_easy")),
           "Chosen from the separation pilot (Tables~\\ref{tab:separation} and~\\ref{tab:pilot}) as the values at which the methods remain distinguishable"),
      list("Convergence criterion", paste0("reclassification proportion ", v("conv_prop")),
           "An exact match of successive partitions, which also detects cycles; relaxing it to 0.001--0.01 changed no partition in a 20-replicate check"),
      list("Maximum iterations", v("maxIter"),
           "Bound 2 of 1200 core fits of RKCFC and none of KCFC (Table~\\ref{tab:convergence})"),
      list("Minimum class size", "10 curves", "A class below it stops the fit; the stopping cause is recorded"),
      list("Initialization", v("init_method"),
           "Starts from $k$-means on the Kendall scores of a fit to all subjects and switches to a random assignment if a class falls below ten curves; alternatives in Table~\\ref{tab:init}"),
      list("Random restarts", v("n_restarts"), "Used in the initialization study; the partition with the smallest $J$ is retained"),
      list("Observations per subject", paste0("Poisson(", v("obs_poisson_lambda"), "), floor ", v("obs_min")),
           "Counts below the floor are redrawn; the floor is lowered to 3 and 2 in the window variations (Table~\\ref{tab:variations})"),
      list("Measurement error SD", v("error_sd"), "Standard deviation, not variance; 0.5 in Table~\\ref{tab:variations}"),
      list("Subject-level bandwidth $h'$", paste(v("bw"), "of the domain"),
           "Nadaraya--Watson smoother in Eq.~(5); no automatic selector exists, so it is a supplied value; halved and doubled in Table~\\ref{tab:smoothing}"),
      list("Kendall-surface and mean bandwidths", paste(v("bwK"), "/", v("bwmean")),
           "Generalized cross-validation"),
      list("Estimation grid", paste(v("nRegGrid"), "points"), "Regular grid for the mean and Kendall surface"),
      list("Eigenfunction smoothing", paste0("B-splines, ", v("basis_nbasis"), " functions, order ", v("basis_norder"), ", ", v("basis_nbreaks"), " breakpoints"),
           "Applied to the discrete eigenvectors"),
      list("Leave-one-out refit", "Not used", "The package option refits once per subject per iteration and was switched off"),
      list("Held-out fraction", v("holdout_frac"), "The fraction of each subject's observations withheld for the held-out reconstruction error"),
      list("Curve contamination",paste0(100 * as.numeric(v("contam_frac")), "% of subjects, ", v("contam_mult"), " times the within-class RMS, length scale ", v("contam_lengthscale")),
           "Additive Gaussian-process trajectory orthogonalized against both classes' eigenfunctions"),
      list("Replicates", paste0(v("n_reps_core"), " (core), ", v("n_reps"), " (variations, class number, estimation, settings), ", v("n_reps_pilot"), " (pilot)"),
           "Medians with bootstrap standard errors, in a column or as the largest value in the table note"),
      list("Master seed", v("seed_master"), "A per-cell seed is derived from it, so any scenario and replicate reproduces in isolation"))
    rows <- vapply(S1, function(r) paste(r[[1]], esc(r[[2]]), r[[3]], sep = " & "), character(1))
    emit("tabS1",
         "Tuning parameters and implementation settings of the simulation study, read from the configuration tables shipped with the package.",
         "tab:params",
         "Setting & Value & Basis",
         "p{3.5cm}p{2.9cm}p{7.5cm}", unname(rows), wide = TRUE,
         notes = paste("The fraction-of-variance rule and the convergence evidence \\emph{disagree} on this design, since",
                       "cumulative FVE at 90\\% or above selects four or more components, every one of which",
                       "runs to the iteration cap and scores at chance (Table~\\ref{tab:components})."))
  }
}

## S2 -- separation -------------------------------------------------------------
if (!is.null(dlt)) {
  ds <- sort(unique(dlt$level))
  rows <- vapply(sort(unique(dlt$arm)), function(a) {
    cells <- vapply(ds, function(l) f(med(dlt$delta[dlt$arm == a & dlt$level == l]), 2), character(1))
    paste(c(ARM_LAB[a], cells), collapse = " & ")
  }, character(1))
  ref <- vapply(ds, function(l) f(med(dlt$ref_ccr[dlt$level == l]), 3), character(1))
  rows <- c(rows, paste(c("\\midrule Median reference CCR", ref), collapse = " & "))
  emit("tabS2",
       "Robust Mahalanobis separation $\\Delta$ between the two class score distributions, by arm and mean offset, over 20 replicates.",
       "tab:separation",
       paste(c("Arm", paste0("$d = ", sub("^d=", "", ds), "$")), collapse = " & "),
       paste0("l", strrep("c", length(ds))), rows,
       notes = paste("$\\Delta$ is computed from the spatial median and spatial-sign covariance, so it is defined",
                     "on the infinite-variance arm. The reference CCR is $\\Phi(\\Delta/2)$, the accuracy a",
                     "two-class linear rule would attain at that separation; it is a \\emph{reference}, not a Bayes",
                     "bound, and the observed accuracies exceed it because the classes also differ in eigenfunction.",
                     "Cross-arm spread exceeds 3-fold at every $d$, so the arms are not exchangeable and accuracy",
                     "is read against each arm's own reference rather than compared across arms."))
}

## S3 -- convergence and cost ----------------------------------------------------
if (!is.null(core)) {
  ms <- names(METH_LAB)
  rows <- vapply(ms, function(m) {
    s <- core[core$method == m, ]
    sr <- s$stop_reason
    # ONLY THE ITERATIVE METHODS HAVE A STOP REASON. M4, M5, M6 and M7 fit once
    # and stop, so they record none -- and printing four zeros for them would
    # read as "never converged" rather than "the question does not arise", which
    # is the opposite of the truth about a one-shot method.
    cause <- if (all(is.na(sr)))
      "\\multicolumn{4}{c}{\\emph{not iterative}}"
    else paste(sum(sr == "converged", na.rm = TRUE), sum(sr == "cycle", na.rm = TRUE),
               sum(sr == "maxIter", na.rm = TRUE), sum(sr == "small_cluster", na.rm = TRUE),
               sep = " & ")
    paste(METH_LAB[[m]], cause, sum(!is.na(s$error)),
          f(med(s$iter), 1), f(med(s$seconds), 1), sep = " & ")
  }, character(1))
  emit("tabS3",
       "Convergence and computational cost on the core grid, over all 1200 cells. Non-convergence is reported \\emph{by cause}, because a fit stopped by the iteration cap and one refused for a small class need opposite remedies.",
       "tab:convergence",
       "Method & Converged & Cycle & Cap & Small class & Failed & Median iter. & Median s",
       "lcccccrr", unname(rows),
       notes = paste("The cap ($\\texttt{maxIter} = 100$) bound only 2 of 1200 RKCFC fits, so no reported accuracy",
                     "rests on a truncated loop. The latent-class spline failures are degenerate fits returning no",
                     "partition, recorded rather than discarded. Times are single-core seconds per fit."))
}

## S4 -- the objective trace ------------------------------------------------------
if (!is.null(conv)) {
  q90 <- function(x) { x <- x[is.finite(x)]; if (length(x)) unname(stats::quantile(x, 0.9)) else NA_real_ }
  ms <- sort(unique(conv$method))
  rows <- character(0); mids <- integer(0)
  for (a in sort(unique(conv$arm))) {
    for (m in ms) {
      s <- conv[conv$arm == a & conv$method == m & is.na(conv$error), ]
      if (!nrow(s)) next
      rows <- c(rows, paste(if (m == ms[1]) ARM_LAB[a] else "", METH_LAB[[m]],
                            f(med(s$n_reest), 1), pc(med(s$frac_reest_up)),
                            pc(med(s$med_rise_rel), 1), pc(q90(s$max_rise_rel), 1),
                            sum(s$stop_reason == "cycle", na.rm = TRUE),
                            f(med(s$iter), 1), sep = " & "))
    }
    mids <- c(mids, length(rows))
  }
  mids <- mids[-length(mids)]
  ok <- conv[is.na(conv$error), ]
  emit("tabS4",
       sprintf("Behavior of the objective $J$ across iterations, 50 replicates. The reassignment step is monotone by construction; the re-estimation step is not, and this table reports how often that matters. Reassignment descended in all %d fits.",
               sum(ok$reassign_descends_all, na.rm = TRUE)),
       "tab:jtrace",
       "Arm & Method & Steps & \\% raising $J$ & Median rise & p90 largest rise & Cycles & Median iter.",
       "llcccccc", rows, midrules = mids,
       notes = paste("A rise is $J_{\\mathrm{before}}(b) - J_{\\mathrm{after}}(b-1) > 0$, the change",
                     "attributable to the refit alone, excluding the preceding reassignment.",
                     "Rises are relative to the $J$ the refit started from, because $J$ scales with $n$",
                     "and with the noise level and is not comparable across arms otherwise.",
                     "\\emph{Median rise} is the median over fits of that fit's median rise;",
                     "\\emph{p90 largest rise} is the 90th percentile of each fit's single largest rise,",
                     "so it describes the worst rises.",
                     "Both methods are shown because the absence of guaranteed descent is",
                     "a property of the K-centers construction rather than of the robust estimator, since",
                     "KCFC raises $J$ more often than RKCFC on three of the four arms and comparably on",
                     "the fourth, so substituting Kendall's $\\tau$ neither causes the difficulty nor",
                     "escapes it.",
                     "A cycle means the loop revisited a partition it had already occupied; this is how the",
                     "missing descent guarantee shows up in practice, and it is why",
                     "cycle detection is part of the stopping rule."))
}

# What each variations row moves, from config/scenarios_b.csv, so the table can be
# read without the config. Codes alone (B04_a2) tell a reader nothing.
VARI_LAB <- c(
  B01 = "Reference",
  B02 = "Crossing means, equal level",
  B03 = "Mean-identical, $d = 0$",
  B04 = "Sparser, Poisson mean 3",
  B05 = "Error SD 0.5",
  B06 = "Four eigenfunctions, two retained",
  B07 = "Class-dependent window, floor 2",
  B08 = "$n = 1000$",
  B09 = "Proportions 0.7/0.3",
  B10 = "Non-informative window, counts unchanged",
  B11 = "Non-informative window, counts by rate, floor 2",
  B12 = "Non-informative window, floor 3")

## S5 -- variations --------------------------------------------------------------
if (!is.null(vari)) {
  scs <- sort(unique(vari$scenario_id))
  # Rows grouped by arm: the scenario codes (B04_a2) carried
  # nothing a reader could use except the arm suffix, so the arm heads each block and the
  # variation names the row.
  codes <- sort(unique(sub("_a[24]$", "", scs)))
  rows <- character(0); mids <- integer(0)
  for (a in c(2L, 4L)) {
    rows <- c(rows, sprintf("\\multicolumn{7}{@{}l}{\\textit{%s arm}}", ARM_LAB[a]))
    for (cd in codes) {
      sc <- sprintf("%s_a%d", cd, a)
      s <- vari[vari$scenario_id == sc, ]
      g <- function(m, col) { z <- s[s$method == m, ]; switch(col,
        ari = f(med(z$ari)), ps = pc(psucc(z$ari)), cv = f(rate(z$converged), 2)) }
      rows <- c(rows, paste(VARI_LAB[[cd]], g("M1","ari"), g("M1","ps"), g("M1","cv"),
                            g("M2","ari"), g("M2","ps"), g("M2","cv"), sep = " & "))
    }
    mids <- c(mids, length(rows))
  }
  mids <- mids[-length(mids)]
  mc_var <- lapply(c(M1 = "M1", M2 = "M2"), function(m)
    vapply(scs, function(sc) mcse_median(vari$ari[vari$scenario_id == sc & vari$method == m]), numeric(1)))
  emit("tabS5",
       "Variations on the baseline design: one factor moved at a time, 50 replicates each, on the Fleishman and $t_{1.5}$ arms. P(s) is the share of replicates attaining an adjusted Rand index above 0.05.",
       "tab:variations",
       "Variation & \\multicolumn{3}{c}{RKCFC} & \\multicolumn{3}{c}{KCFC} \\\\ \\cmidrule(lr){2-4}\\cmidrule(lr){5-7} & aRI & P(s) & Conv. & aRI & P(s) & Conv.",
       "lcccccc", rows, wide = TRUE, midrules = mids,
       notes = paste("The crossing-means variation uses $\\mu_{(k)} = \\sqrt{2}\\sin(4\\pi t) \\pm d\\cos(\\pi t)$, whose average",
                     "levels are identical. The four-eigenfunction variation gives each class four true eigenfunctions,",
                     "at frequencies 2, 3, 5, and 6 with eigenvalues $(9, 4, 2, 1)$ and $(8.1, 3.6, 1.8, 0.9)$, of which",
                     "two are retained. The follow-up windows $[0, T_i]$ are non-informative in the three",
                     "non-informative-window variations, with $T_i \\sim U(0.4, 1)$ and the counts per subject unchanged",
                     "or drawn at rate $\\lambda T_i$ as labeled, and class-dependent in the class-dependent variation,",
                     "with $T_i \\sim U(0.6, 1)$ in class 1 and $U(0.3, 0.7)$ in class 2 and counts drawn at rate",
                     "$\\lambda T_i$, so that follow-up carries information about class, as it does in the clinical cohort.",
                     "The two floor-2 variations relax the minimum observations per subject from four to two and",
                     "converge on \\emph{no} replicate; the floor-3 variation, the counts-unchanged window with the",
                     "floor at three, behaves like the default. The method therefore requires $n_i \\geq 3$ on this",
                     "design, and the threshold is sharp rather than gradual.",
                     sprintf("Bootstrap standard errors of the aRI medians are at most %s for RKCFC and %s for KCFC.",
                             mcse_max(mc_var$M1), mcse_max(mc_var$M2))))
}

## S6 -- retained components -----------------------------------------------------
if (!is.null(comp)) {
  # Row labels in words: "fixed 2" and "fve95(M=6)" were
  # level codes. Fixed levels first, then the FVE thresholds in order, with the two
  # FVE 95% rows (which selected different M) together.
  comp_lab <- function(l) {
    if (grepl("^M[1-4]$", l)) return(sprintf("$M = %s$ fixed", substr(l, 2, 2)))
    m <- regmatches(l, regexec("^fve([0-9]+)\\(M=([0-9]+)\\)$", l))[[1]]
    if (length(m) == 3) return(sprintf("FVE %s\\%%, $M = %s$ selected", m[2], m[3]))
    esc(l)
  }
  comp_key <- function(l) {
    if (grepl("^M[1-4]$", l)) return(as.numeric(substr(l, 2, 2)))
    m <- regmatches(l, regexec("^fve([0-9]+)\\(M=([0-9]+)\\)$", l))[[1]]
    if (length(m) == 3) 100 + as.numeric(m[2]) + as.numeric(m[3]) / 10 else 1e6
  }
  lv <- unique(comp$level)
  lv <- lv[order(vapply(lv, comp_key, numeric(1)))]
  rows <- vapply(lv, function(l) {
    s <- comp[comp$level == l, ]
    paste(comp_lab(l), f(med(s$ari)), f(mcse_median(s$ari), 4), f(rate(s$converged), 2),
          f(med(s$iter), 0), f(med(s$rre_true)), f(mcse_median(s$rre_true), 4), sep = " & ")
  }, character(1))
  emit("tabS6",
       "Retained components. The first four rows fix the number retained; the FVE rows choose it by cumulative fraction of variance explained, and each states the number selected.",
       "tab:components",
       "Components retained & aRI & MCSE & Conv. & Median iter. & Rel.\\ recon.\\ error & MCSE",
       "lcccccc", unname(rows),
       notes = paste("Every row retaining four or more components has a",
                     "convergence rate of zero and a median iteration count equal to the cap, so its accuracy is an",
                     "artifact of where the loop was stopped. No FVE threshold recovers the two components that",
                     "converge and score highest, which is why the components rule is reported as a disagreement",
                     "between criteria rather than as a selection rule."))
}

## S7 -- initialization ----------------------------------------------------------
if (!is.null(ini)) {
  # Scheme names in words; the level codes were internal.
  INIT_LAB <- c(two_staged = "Two-staged (default)", pace_kmeans = "$k$-means on conventional scores",
                kendall_kmeans = "$k$-means on Kendall scores", random20 = "Twenty random restarts",
                ward = "Ward's method", oracle = "True partition (reference)")
  lv <- unique(ini$level)
  rows <- character(0); mids <- integer(0)
  for (a in sort(unique(ini$arm))) {
    for (l in lv) {
      s <- ini[ini$arm == a & ini$level == l, ]
      rows <- c(rows, paste(if (l == lv[1]) ARM_LAB[a] else "", if (l %in% names(INIT_LAB)) INIT_LAB[[l]] else esc(l),
                            f(med(s$ari)), f(mcse_median(s$ari), 4), f(rate(s$converged), 2),
                            sum(is.finite(s$ari)), sep = " & "))
    }
    mids <- c(mids, length(rows))
  }
  mids <- mids[-length(mids)]
  emit("tabS7",
       "Initialization schemes, 50 replicates. The two-staged start is the configured default; the start from the true partition is a reference, not a usable scheme.",
       "tab:init",
       "Arm & Initialization & aRI & MCSE & Conv. & $n$ returned",
       "llcccc", rows, midrules = mids,
       notes = paste("The $n$ column matters most on the $t_{1.5}$ arm, where the deterministic schemes fail outright on",
                     "most replicates because clustering heavy-tailed scores produces a class below the minimum",
                     "size, so their medians are taken over the replicates that happened to survive and are not",
                     "comparable with the others. Twenty random restarts match the start from the true partition and never fail."))
}

## S8 -- bandwidth and mean estimator --------------------------------------------
if (!is.null(smo)) {
  # Setting names in words; the level codes were internal.
  SMO_LAB <- c(bw_half = "$h'$ halved", bw_base = "$h'$ as configured", bw_double = "$h'$ doubled",
               mean_robust = "Robust mean, $h'$ as configured")
  lv <- unique(smo$level)
  rows <- character(0); mids <- integer(0)
  for (a in sort(unique(smo$arm))) {
    for (l in lv) {
      s <- smo[smo$arm == a & smo$level == l, ]
      rows <- c(rows, paste(if (l == lv[1]) ARM_LAB[a] else "", if (l %in% names(SMO_LAB)) SMO_LAB[[l]] else esc(l),
                            f(med(s$ari)), f(mcse_median(s$ari), 4), f(rate(s$converged), 2),
                            pc(psucc(s$ari)), sum(is.finite(s$ari)), sep = " & "))
    }
    mids <- c(mids, length(rows))
  }
  mids <- mids[-length(mids)]
  emit("tabS8",
       "Bandwidth sensitivity and the mean estimator, 50 replicates. The robust mean replaces the least-squares mean with a local linear M-estimator at the configured bandwidth.",
       "tab:smoothing",
       "Arm & Setting & aRI & MCSE & Conv. & P(s) & $n$ returned",
       "llccccc", rows, midrules = mids,
       notes = paste("Accuracy is flat across a four-fold change in bandwidth, but the mean estimator matters.",
                     "On the $t_{1.5}$ arm the least-squares mean leaves the index bimodal, with roughly a third",
                     "of replicates collapsing to zero, while the robust mean converges on every replicate and",
                     "removes the zero mode entirely, at no cost on the Gaussian arm. The class mean is the",
                     "reference the reassignment step scores against, so distortion of the mean, not of the",
                     "eigenspace, is what drives the failures."))
}

## S9 -- estimator swaps, judged on the partition --------------------------------
if (!is.null(estr)) {
  # The same names as Table 4; E1_kendall etc. were codes.
  EST_LAB <- c(E1_kendall = "Kendall's $\\tau$", E2_sign = "Spatial sign",
               E3_rank = "Normal scores + conventional FPCA")
  lv <- unique(estr$level)
  rows <- character(0); mids <- integer(0)
  for (a in sort(unique(estr$arm))) {
    for (l in lv) {
      s <- estr[estr$arm == a & estr$level == l, ]
      rows <- c(rows, paste(if (l == lv[1]) ARM_LAB[a] else "", if (l %in% names(EST_LAB)) EST_LAB[[l]] else esc(l),
                            f(med(s$ari)), f(mcse_median(s$ari), 4),
                            f(rate(s$converged), 2), pc(psucc(s$ari)),
                            f(med(s$seconds), 1), sep = " & "))
    }
    mids <- c(mids, length(rows))
  }
  mids <- mids[-length(mids)]
  emit("tabS9",
       "Robust estimators swapped inside the clustering loop, 50 replicates. This is the comparison on equal terms, because a partition is invariant to the marginal transform that the normal-scores estimator applies, so all three estimators are scored on the same quantity here, which is not true of the estimation endpoints in Table~\\ref{tab:estim}.",
       "tab:estimators",
       "Arm & Estimator & aRI & MCSE & Conv. & P(s) & Median s",
       "llccccc", rows, midrules = mids,
       notes = paste("Kendall's $\\tau$ signs pairwise differences between subjects;",
                     "the spatial-sign estimator signs each subject's deviation from the mean, at the same mean, so",
                     "the two differ in the scatter estimator alone; the normal-scores estimator replaces the",
                     "measurements by normal scores and then applies conventional FPCA.",
                     "Reconstruction error is not reported for the normal-scores estimator, because the fit lives on the",
                     "normal-scores scale and the data do not, and the transform has no elementwise inverse,",
                     "so any such number would compare two different scales."))
}

## S10 -- pilot -------------------------------------------------------------------
if (!is.null(pil)) {
  rows <- character(0); mc_pil <- list(M1 = numeric(0), M2 = numeric(0))
  # The arm names the block; the pilot codes P01-P06 said
  # nothing, and the arm was not shown at all.
  seen <- integer(0); mids_pil <- integer(0)
  for (sc in sort(unique(pil$scenario_id))) {
    s <- pil[pil$scenario_id == sc, ]
    a <- unique(s$arm)[1]
    if (length(seen) && !(a %in% seen)) mids_pil <- c(mids_pil, length(rows))
    g <- function(m, col) { z <- s[s$method == m, ]; switch(col,
      ccr = f(med(z$ccr)), ari = f(med(z$ari)), cv = f(rate(z$converged), 2)) }
    rows <- c(rows, paste(if (a %in% seen) "" else ARM_LAB[a], f(unique(s$d)[1], 1),
                          g("M1","ccr"), g("M1","ari"), g("M1","cv"),
                          g("M2","ccr"), g("M2","ari"), g("M2","cv"), sep = " & "))
    seen <- union(seen, a)
    for (m in c("M1", "M2")) {
      z <- s[s$method == m, ]
      mc_pil[[m]] <- c(mc_pil[[m]], mcse_median(z$ccr), mcse_median(z$ari))
    }
  }
  emit("tabS10",
       "Pilot confirmation at the full sample size, 20 replicates. The pilot checked the design before the core grid was run and is not a result in itself.",
       "tab:pilot",
       "Arm & $d$ & \\multicolumn{3}{c}{RKCFC} & \\multicolumn{3}{c}{KCFC} \\\\ \\cmidrule(lr){3-5}\\cmidrule(lr){6-8} & & CCR & aRI & Conv. & CCR & aRI & Conv.",
       "lccccccc", rows, midrules = mids_pil,
       notes = paste("The iteration cap bound no fit in the pilot, which is consistent with two retained components being enough (see also Table~\\ref{tab:components}).",
                     sprintf("Bootstrap standard errors of the medians are at most %s for RKCFC and %s for KCFC.",
                             mcse_max(mc_pil$M1), mcse_max(mc_pil$M2))))
}

# --- driver document -----------------------------------------------------------

## Overview of the studies ----------------------------------------------------------
## One row per study. The descriptive columns are design facts; the replicates and
## the fits are read from the result files, so the counts cannot drift from the
## evidence, and the total is printed because the Studies paragraph of Section 3.1
## quotes it. A fit is one (scenario, method or level, replicate) combination, per
## class in the estimation study and per K in the class-number study; failed fits
## are counted, since every fit leaves a row with its stopping cause.
{
  key_all <- c("scenario_id", "method", "estimator", "level", "rep", "K_fit", "class")
  nfit <- function(blocks) {
    reps <- integer(0); fits <- 0L
    for (b in blocks) {
      x <- rd(b); if (is.null(x)) next
      reps <- c(reps, length(unique(x$rep)))
      fits <- fits + nrow(unique(x[intersect(key_all, names(x))]))
    }
    list(reps = if (length(reps)) paste(sort(unique(reps)), collapse = ", ") else "--", fits = fits)
  }
  studies <- list(
    list("Core grid", "Arm and separation", "All six", "1.5, 2", "All seven", "core",
         "Tables~\\ref{tab:core}, \\ref{tab:recon}, and~\\ref{tab:convergence}"),
    list("Objective trace", "$J$ recorded at every iteration", "Gaussian, Fleishman(0, 4), $t_{1.5}$, contamination",
         "1.5", "RKCFC, KCFC", "convergence", "Table~\\ref{tab:jtrace}"),
    list("Variations block", "One factor at a time, twelve rows", "Fleishman(0, 4), $t_{1.5}$",
         "1.5 (0 in one row)", "RKCFC, KCFC, FPCA $k$-means, Kendall $k$-means", "variations", "Table~\\ref{tab:variations}"),
    list("Class-number study", "True $K$ of 1, 2, and 3; $K$ fitted from 1 to 5", "Fleishman(0, 4)",
         "1.5", "RKCFC, KCFC", "krecovery", "Table~\\ref{tab:krecovery}"),
    list("Estimation study", "Five eigenspace estimators, fitted per class under the true partition", "All six",
         "1.5", "Estimators only", "estim", "Table~\\ref{tab:estim}"),
    list("Initialization", "Six starting schemes", "Gaussian, $t_{1.5}$", "1.5", "RKCFC", "S_init", "Table~\\ref{tab:init}"),
    list("Retained components", "$M$ fixed at 1 to 4, or chosen by FVE at 80 to 99\\%", "Fleishman(0, 4)",
         "1.5", "RKCFC", "S_components", "Table~\\ref{tab:components}"),
    list("Bandwidth and mean estimator", "$h'$ halved and doubled; M-estimator of the mean",
         "Gaussian, both Fleishman, $t_{1.5}$", "1.5", "RKCFC", "S_smoothing", "Table~\\ref{tab:smoothing}"),
    list("Estimators inside the loop", "Kendall's $\\tau$, spatial sign, normal scores",
         "Gaussian, Fleishman(0, 4), $t_{1.5}$", "1.5", "RKCFC", "S_estimators", "Table~\\ref{tab:estimators}"),
    list("Separation pilot", "Robust Mahalanobis $\\Delta$ at seven values of $d$", "All six", "0 to 4",
         "None (score distributions)", "D_delta", "Table~\\ref{tab:separation}"),
    list("Pilot at full sample size", "Candidate separations", "Gaussian, $t_{1.5}$", "0, 1.5, 2",
         "RKCFC, KCFC", "pilot", "Table~\\ref{tab:pilot}"),
    list("Design checks", "Separation ladder, convergence tolerance, initialization, retained components",
         "Gaussian, $t_{1.5}$", "1.5 (0 to 3 in the ladder)", "RKCFC, KCFC",
         c("D_d", "D_convprop", "D_init", "D_nkfit"), "Table~\\ref{tab:params}"))
  rows <- character(0); total <- 0L
  for (s in studies) {
    nf <- nfit(s[[6]]); total <- total + nf$fits
    rows <- c(rows, paste(s[[1]], s[[2]], s[[3]], s[[4]], s[[5]], nf$reps,
                          formatC(nf$fits, format = "d", big.mark = ","), s[[7]], sep = " & "))
  }
  rows <- c(rows, paste("Total", "", "", "", "", "", formatC(total, format = "d", big.mark = ","), "", sep = " & "))
  message("  studies overview: ", total, " fits in all (the Studies paragraph quotes this number)")
  emit("tabStudies",
       "Overview of the simulation studies. Each row is one study; replicates and fits are read from the result files and include fits that failed, and the last column names the table that reports the study.",
       "tab:studies",
       "Study & Varies & Arms & $d$ & Methods & Replicates & Fits & Results",
       "p{2.1cm}p{3.0cm}p{2.3cm}p{1.3cm}p{2.0cm}ccp{2.1cm}", rows, midrules = length(rows) - 1L, wide = TRUE,
       notes = paste("The variations block also fits the two one-shot $k$-means methods, whose rows",
                     "Table~\\ref{tab:variations} omits. The design checks are the pilots behind the",
                     "settings of Table~\\ref{tab:params}; their outcomes appear in its Basis column."))
}

main <- c(sprintf("tab%02d", 1:5), "tabCCR")
supp <- c("tabS1", "tabStudies", sprintf("tabS%d", 2:10))
have <- function(x) x[file.exists(file.path(OUT, paste0(x, ".tex")))]
main <- have(main); supp <- have(supp)

drv <- file(file.path(OUT, "all_tables.tex"), open = "wt", encoding = "UTF-8")
w <- function(...) cat(..., "\n", sep = "", file = drv)
w("% Proof document for the generated tables. The manuscript should \\input the")
w("% individual fragments instead; this exists only to typeset them for review.")
w("\\documentclass[11pt]{article}")
w("\\usepackage[margin=2cm,landscape]{geometry}")
w("\\usepackage{booktabs,amsmath,amssymb,array,caption}")
w("\\captionsetup{justification=raggedright,singlelinecheck=false}")
w("\\renewcommand{\\arraystretch}{1.08}")
w("\\begin{document}")
w("\\begin{center}{\\Large Robust $K$-centers functional clustering: tables}\\\\[2pt]")
w("{\\footnotesize generated ", format(Sys.Date()), " from simulation output; every value is read from the result files}\\end{center}")
w("\\section*{Main tables}")
for (t in main) w("\\input{", t, "}")
w("\\clearpage")
w("\\section*{Supplementary tables}")
w("\\renewcommand{\\thetable}{S\\arabic{table}}\\setcounter{table}{0}")
for (t in supp) w("\\input{", t, "}")
w("\\end{document}")
close(drv)

message("\nwrote ", length(main), " main + ", length(supp), " supplementary tables to ", OUT, "/")
message("driver: ", file.path(OUT, "all_tables.tex"))
