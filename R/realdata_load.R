# THE COHORT: CEA over months 0-36, loaded once, defensively.
#
#   source(file.path("R", "realdata_load.R"))
#   dat <- load_cea36()
#
# ONE SOURCE, STATED ONCE. The file holds seven objects and three of them could
# plausibly be called "the CEA data". The cohort is:
#
#   marker_yn_list_36$CEA      CEA measurements, 0-36 months, 1508 subjects
#   restricted to dt_base_36   the 1405 with a baseline and outcome record
#
# Verified 2026-09-09: that restriction reproduces `CEA_199_242_36m[[1]]` exactly,
# row for row, so the two routes into the cohort agree and neither is a silent
# alternative to the other. The 103 excluded subjects have CEA trajectories but no
# baseline record and no survival time -- in `dt_base_12` either -- so they could
# be clustered but could never appear in F2, F3 or any table that follows. A
# partition fitted on subjects who then vanish would be inconsistent, and the
# exclusions are accounted for regardless.
#
# THE FILE CONTAINS DIRECT IDENTIFIERS AND THIS LOADER MUST NEVER RETURN THEM.
# The raw table in the source file carries the hospital admission number and the
# patient name in
# its first two columns. The analysis sets carry only an integer `ID`. This loader
# touches the analysis sets alone and asserts, below, that nothing it returns
# holds a non-ASCII column name -- a cheap structural check that would catch the
# raw table being wired in by mistake.

# WHERE THE COHORT LIVES, on either machine.
#
# THE DATA ARE NOT DISTRIBUTED with this code. Set the environment variable
# RKCFC_DATA to the path of the source file; nothing here copies or writes
# patient data.
DATA_CANDIDATES <- c(
  Sys.getenv("RKCFC_DATA", ""))

default_data_file <- function() {
  hit <- DATA_CANDIDATES[nzchar(DATA_CANDIDATES) & file.exists(DATA_CANDIDATES)]
  if (!length(hit))
    stop("cohort file not found. Set RKCFC_DATA to its path.\n  tried:\n    ",
         paste(DATA_CANDIDATES[nzchar(DATA_CANDIDATES)], collapse = "\n    "),
         call. = FALSE)
  hit[1]
}

DATA_FILE_DEFAULT <- NULL   # resolved at call time, not at source time

load_cea36 <- function(path = NULL, scale = c("raw", "log")) {
  scale <- match.arg(scale)
  if (is.null(path)) path <- default_data_file()
  if (!file.exists(path))
    stop("cohort file not found: ", normalizePath(path, mustWork = FALSE),
         "\n  pass `path=` or set RKCFC_DATA.", call. = FALSE)

  x <- readRDS(path)
  need <- c("marker_yn_list_36", "dt_base_36")
  miss <- need[!need %in% names(x)]
  if (length(miss))
    stop("the file does not hold ", paste(miss, collapse = " and "),
         ". Is this the right extract?", call. = FALSE)

  cea <- x$marker_yn_list_36[["CEA"]]
  base <- x$dt_base_36
  if (is.null(cea)) stop("marker_yn_list_36 has no CEA element.", call. = FALSE)

  # NO IDENTIFIERS, ASSERTED RATHER THAN ASSUMED. Column names in the analysis
  # sets are ASCII; the raw table's are not. This is a structural guard against
  # someone later pointing this loader at the raw table.
  for (nm in list(names(cea), names(base)))
    if (any(grepl("[^\\x01-\\x7F]", nm, perl = TRUE)))
      stop("a returned table carries non-ASCII column names, which in this file ",
           "means DIRECT IDENTIFIERS. Refusing to return it.", call. = FALSE)

  n_all <- length(unique(cea$ID))
  keep  <- cea$ID %in% base$ID
  cea   <- cea[keep, , drop = FALSE]
  n_kept <- length(unique(cea$ID))

  cea <- cea[order(cea$ID, cea$time), , drop = FALSE]
  ids <- unique(cea$ID)

  # The measurement on the requested scale. `yit_lg` was verified to be the
  # NATURAL log of `yit` -- not log1p, not log10 -- so the back-transform is exp().
  y <- if (scale == "raw") cea$yit else cea$yit_lg

  Lt <- split(cea$time, factor(cea$ID, levels = ids))
  Ly <- split(y,        factor(cea$ID, levels = ids))
  names(Lt) <- names(Ly) <- NULL

  ni <- lengths(Ly)
  # The floor the simulation established: the method never converged at n_i = 2
  # and behaved normally at 3 (variations B10/B11/B12, 2026-09-08). Assert it
  # here so a future extract that violates it is caught before a fit is attempted
  # rather than diagnosed from a wall of non-convergence.
  if (min(ni) < 3L)
    stop(sprintf(paste("%d subjects have fewer than 3 observations (minimum %d).",
                       "The simulation shows the method does not converge below 3."),
                 sum(ni < 3L), min(ni)), call. = FALSE)

  base <- base[match(ids, base$ID), , drop = FALSE]

  list(Ly = Ly, Lt = Lt, ids = ids, base = base, scale = scale,
       interval = c(0, max(cea$time)),
       excluded = n_all - n_kept, n_source = n_all,
       raw = cea)
}

# Descriptives of the sparse design: the distribution of the number of
# measurements per patient and of the observation times.
cohort_descriptives <- function(d) {
  ni <- lengths(d$Ly); tt <- unlist(d$Lt)
  list(n_subjects = length(ni), n_obs = length(tt),
       n_source = d$n_source, excluded = d$excluded,
       ni_median = stats::median(ni), ni_iqr = stats::quantile(ni, c(.25, .75)),
       ni_range = range(ni), ni_table = table(ni),
       time_median = stats::median(tt), time_iqr = stats::quantile(tt, c(.25, .75)),
       time_range = range(tt))
}

# ===========================================================================
# THE REAL-DATA ANALYSIS SETTINGS, CHOSEN BY SWEEP AND STATED ONCE
#
# Every one of these was a measured choice, made 2026-09-09 after the first
# Study F run failed completely. They live here rather than in config/ for a
# specific reason: `source_hash()` hashes config/*.csv, so adding them there
# would change the provenance stamp of every SIMULATION result and make new
# runs refuse to combine with the banked ones. This file is outside that set.
# Each real-data runner therefore passes THIS FILE to source_hash() as `extra`,
# so real-data results still stamp the settings that produced them.
#
# REAL_TMAX = 18. Tmax 36 and 24 could not be fitted AT ANY bandwidth tried
#   (F0B: 1/3 and 0/3 seeds), and 12 failed too. 18 fitted 3/3. It also costs
#   almost nothing -- 10 subjects of 1405, 93.1 percent of observations -- and
#   it RAISES F2's power, from 127 events after a month-36 landmark to 187
#   after a month-18 one.
#
# REAL_BWK = 6, and this is the parameter that actually binds. GCV selects
#   1.4-2.4 here and produces a DEGENERATE 5 percent class; a fixed 1, 2, 3 or 4
#   fails outright; 6 and 8 converge on every seed and split 14/86. The error
#   that killed the first run names `userBwCov`, which on the robust branch is
#   bwK -- handed to fdapace::Lwls2D at RKCFC/R/KFPCA.R:188. The same failure
#   costs 2.92 percent of simulated fits in the variations block, so it is a
#   property of the estimator, not of this cohort.
#
#   IT MUST BE REPORTED AS A CHOSEN VALUE. bwK = 6 is a third of the domain and
#   3-4x GCV's selection, and in-sample error is 30 percent worse than at GCV.
#   "GCV" is no longer an available answer to how the bandwidth was selected,
#   because GCV is what failed.
#
# REAL_BW = 3.6 unchanged. F0 crossed it over 3.6, 5.4, 7.2 and 9.0 and
#   feasibility did not move by a single cell. It is the subject-level
#   Nadaraya-Watson bandwidth and never enters the failing call.
REAL_TMAX <- 18
REAL_BW   <- 3.6
REAL_BWK  <- 6

# TIED, NOT COINCIDENTALLY EQUAL. The landmark must sit at the END of the
# exposure window: the grouping is not known until the trajectory is complete,
# which is the whole point of F2. Writing 18 twice would let one be changed
# without the other and reintroduce the guarantee-time bias.
REAL_LANDMARK <- REAL_TMAX

# THE PRIMARY PARTITION, named once. F2 and F3 both consume it and used to spell
# it separately; a switch applied to one file and not the other would put the
# landmark analysis and the Cox models on different partitions without any
# error.
#
# "log_M1", DECIDED 2026-09-10 on read_realdata_partitions.R, and on criteria
# fixed before any survival number was seen:
#
#   * REPRODUCIBLE. Twenty log-scale starts agree with the best at aRI >= 0.949
#     (median 0.997). The raw-scale starts do not (0.534), and KCFC's on either
#     scale are a lottery (0.35, 0.22).
#   * STABLE, OR AT LEAST NOT UNSTABLE. The raw objective's minimum, 58/1337,
#     dissolves under resampling (Jaccard 0.28); the raw best-of-20, 94/1301,
#     scores 0.48. Squared error on CEA with excess kurtosis 501 is minimised by
#     isolating the extreme patients, and that group is not a cluster.
#   * A SHAPE SPLIT. Binned log-CEA class means: level share 9 percent, one class
#     high-then-rebounding, the other flat. Every raw candidate is a level split
#     (55-76 percent) that a threshold on peak CEA reproduces at AUC 0.96-0.98;
#     log_M1's AUC is 0.87.
#
# The raw-scale analysis stays as the SENSITIVITY, its null survival result
# included. Changing the scale after seeing results would be indefensible;
# changing it on these three criteria is what the criteria are for, and the
# paper must say which it was.
#
# F4 was fitted on the raw scale and must be re-run on the log scale before its
# Jaccard describes this partition; run_realdata_f4.R now follows this constant.
REAL_PRIMARY <- "log_M1"

# Truncate to [0, Tmax] and RE-APPLY the n_i >= 3 floor, because cutting the
# domain is exactly what pushes a subject below it.
truncate_cea36 <- function(d, Tmax = REAL_TMAX) {
  k  <- lapply(d$Lt, function(x) x <= Tmax)
  Ly <- Map(function(y, kk) y[kk], d$Ly, k)
  Lt <- Map(function(t, kk) t[kk], d$Lt, k)
  ok <- lengths(Ly) >= 3L
  d$Ly <- Ly[ok]; d$Lt <- Lt[ok]; d$ids <- d$ids[ok]
  d$base <- d$base[ok, , drop = FALSE]
  if (!is.null(d$raw)) d$raw <- d$raw[d$raw$ID %in% d$ids & d$raw$time <= Tmax, , drop = FALSE]
  d$interval <- c(0, Tmax)
  d$n_dropped <- sum(!ok)
  d
}

# The robust options, built once so no runner can drift from another.
realdata_optns <- function(cfg, interval = c(0, REAL_TMAX),
                           bw = REAL_BW, bwK = REAL_BWK) {
  list(interval = interval, nK = gint(cfg, "nK_fit"), bw = bw,
       bwK = bwK, bwmean = cfg$global$bwmean,
       nRegGrid = gint(cfg, "nRegGrid"),
       basis = fda::create.bspline.basis(
         interval, nbasis = gint(cfg, "basis_nbasis"),
         norder = gint(cfg, "basis_norder"),
         breaks = seq(interval[1], interval[2],
                      length.out = gint(cfg, "basis_nbreaks"))))
}

# PEAK CEA, the comparator F0C made necessary. A single summary reproduced the
# partition at AUC 0.961, so the open question is not whether the clusters are
# real -- F0C showed they differ in shape -- but whether they BEAT this number
# at predicting survival. F2 and F3 carry it as a comparator, not a covariate.
peak_cea <- function(d) vapply(d$Ly, max, numeric(1))
