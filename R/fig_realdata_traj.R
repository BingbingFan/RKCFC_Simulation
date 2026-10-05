# FIGURE 1 (and Figure S2) -- smoothed class trajectories of log CEA, by partition.
#
#   Rscript R/fig_realdata_traj.R [f1dir] [outfile]
#
# Two switches, read from the environment so that the manuscript command stays
# the same. Of the four versions rendered on 2026-10-03 (log or original y-axis,
# binned means or smooth), the manuscript uses the log y-axis with the smooth, and the
# defaults are that choice; the others remain available.
#   RKCFC36_FIG_Y      "log" (default) draws log CEA (natural logarithm) on the
#                      y-axis; "original" draws CEA in ng/mL.
#   RKCFC36_FIG_CURVE  "smooth" (default) draws what ggplot2's geom_smooth()
#                      draws by default for more than 1000 points: a GAM with a
#                      cubic shrinkage spline, y ~ s(t, bs = "cs"), fitted by
#                      REML to the class's observations (mgcv::gam), with its
#                      pointwise 95 percent band (fit +/- 1.96 SE); "mean" draws
#                      the mean of the observations in 2-month bins within each
#                      class, with pointwise 95 percent bands (mean +/- 1.96 SE).
#
# WHAT IS DRAWN, AND WHAT IS NOT. The curves are EMPIRICAL summaries of the
# observations in each class, computed identically for both methods, so the
# panels are comparable. They are not the fitted class mean functions: F1.rds
# stores partitions only (no fitted objects, no patient-level data), and
# refitting here would draw a model quantity beside a data quantity without
# saying which was which. The caption says so.
#
# The dotted reference line is 5 ng/mL, the usual clinical cut-off, drawn on
# whichever scale the axis uses, so a reader can see which class sits above it
# and when. Panel titles name the method only; the caption gives the scale of
# the fit. The x-axis is labeled at 0, 6, 12, 18, 24 and 36 months.
#
# NO IDENTIFIERS. Data come through load_cea36(), which returns the analysis set
# only and asserts no identifier column is present; everything drawn is a
# class-level summary.

source(file.path("R", "realdata_load.R"))

args  <- commandArgs(trailingOnly = TRUE)
args  <- args[args != "--args"]
f1dir <- if (length(args) >= 1) args[1] else file.path("results_real", "F1")
out   <- if (length(args) >= 2) args[2] else file.path("..", "doc", "V2", "figures", "fig1_trajectories.png")
yscale <- match.arg(Sys.getenv("RKCFC36_FIG_Y", "log"), c("log", "original"))
curve  <- match.arg(Sys.getenv("RKCFC36_FIG_CURVE", "smooth"), c("smooth", "mean"))

f1   <- readRDS(file.path(f1dir, "F1.rds"))
if (isTRUE(any(f1$table$smoke))) stop("F1.rds is a smoke run.", call. = FALSE)
Tmax <- if (!is.null(f1$Tmax)) f1$Tmax else REAL_TMAX
d    <- truncate_cea36(load_cea36(scale = if (yscale == "log") "log" else "raw"), Tmax)
n    <- length(d$Ly)
tt   <- unlist(d$Lt); yy <- unlist(d$Ly); who <- rep(seq_len(n), lengths(d$Lt))
# Two-month bins for the "mean" variant: the SAME binning read_realdata_partitions.R
# uses, so a binned value quoted anywhere is the number on the page.
br   <- seq(0, Tmax, by = 2); mids <- (br[-1] + br[-length(br)]) / 2
bin  <- cut(tt, br, include.lowest = TRUE)
grid <- seq(0, Tmax, length.out = 241)
xbr  <- c(0, 6, 12, 18, 24, 36); xbr <- xbr[xbr <= Tmax]

class_curves <- function(g) {
  cls <- g[who]
  lapply(sort(unique(g)), function(k) {
    if (curve == "mean") {
      m <- tapply(yy[cls == k], bin[cls == k], mean)[levels(bin)]
      s <- tapply(yy[cls == k], bin[cls == k], function(v) stats::sd(v) / sqrt(length(v)))[levels(bin)]
      list(x = mids, m = as.numeric(m), lo = as.numeric(m - 1.96 * s), hi = as.numeric(m + 1.96 * s))
    } else {
      dk  <- data.frame(t = tt[cls == k], y = yy[cls == k])
      fit <- mgcv::gam(y ~ s(t, bs = "cs"), data = dk, method = "REML")
      p   <- stats::predict(fit, newdata = data.frame(t = grid), se.fit = TRUE)
      list(x = grid, m = as.numeric(p$fit), lo = as.numeric(p$fit - 1.96 * p$se.fit),
           hi = as.numeric(p$fit + 1.96 * p$se.fit))
    }
  })
}
cols <- c("#c0392b", "#2e86c1", "#27ae60", "#8e44ad")
# THE PARTITIONS DRAWN are the primary fit's and its KCFC counterpart on the same
# scale (the record's `primary`: "raw_M1" since 2026-10-03, "log_M1" before),
# whichever scale the y-axis uses, so the two scales' partitions are drawn alike.
prim  <- if (!is.null(f1$primary)) f1$primary else "log_M1"
prim2 <- sub("M1$", "M2", prim)
scale_lab <- if (grepl("^raw", prim)) "original scale" else "log scale"
g1 <- f1$fits[[prim]]; g2 <- f1$fits[[prim2]]
if (is.null(g1) || is.null(g2)) stop("F1.rds lacks ", prim, " or ", prim2, call. = FALSE)
curves <- list(RKCFC = class_curves(g1), KCFC = class_curves(g2))
yl <- range(unlist(lapply(curves, function(cc) lapply(cc, function(z) c(z$lo, z$hi)))), na.rm = TRUE)
yl <- if (yscale == "log") {
  c(floor(yl[1] * 10) / 10 - 0.05, ceiling(yl[2] * 10) / 10 + 0.15)
} else {
  c(min(0, floor(yl[1])), yl[2] + 0.12 * diff(yl))
}
ref  <- if (yscale == "log") log(5) else 5
ylab <- if (yscale == "log") "log CEA (ng/mL)" else "CEA (ng/mL)"

panel <- function(cc, g, labels, title) {
  plot(NA, xlim = c(0, Tmax), ylim = yl, xlab = "Months after surgery",
       ylab = ylab, main = title, las = 1, xaxt = "n")
  axis(1, at = xbr)
  abline(h = ref, lty = 3, col = "grey40")
  # left-anchored just inside the axis: with pos = 3 the label was centred on
  # x = 0.2 and the "5" fell off the left edge of the plot
  text(0.3, ref, "5 ng/mL", adj = c(0, -0.35), cex = 0.75, col = "grey40")
  for (k in seq_along(cc)) {
    z <- cc[[k]]; ok <- is.finite(z$m) & is.finite(z$lo)
    polygon(c(z$x[ok], rev(z$x[ok])), c(z$lo[ok], rev(z$hi[ok])),
            col = grDevices::adjustcolor(cols[k], 0.18), border = NA)
    lines(z$x[ok], z$m[ok], col = cols[k], lwd = 2.2)
    if (curve == "mean") points(z$x[ok], z$m[ok], col = cols[k], pch = 16, cex = 0.7)
  }
  legend("topright", legend = sprintf("%s (n = %d)", labels, as.integer(table(g))),
         col = cols[seq_along(cc)], lwd = 2.2, bty = "n", cex = 0.85)
}

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
grDevices::png(out, width = 2400, height = 1050, res = 300)
graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.5, 1))
# Class names from the record when it has them (the 36-month one does, written by
# R/realdata36_select.R); the 18-month record keeps the names it was published with.
cn <- if (!is.null(f1$class_names)) f1$class_names else c("rebounding", "flat")
panel(curves$RKCFC, g1, cn, "(a) RKCFC")
# KCFC's classes have no names; the text calls them the smaller and the larger class,
# so the legend does the same, in class order.
size_lab <- function(g) { n <- as.integer(table(g)); if (n[1] <= n[2]) c("smaller class", "larger class") else c("larger class", "smaller class") }
panel(curves$KCFC,  g2, size_lab(g2), "(b) KCFC")
invisible(grDevices::dev.off())

# The numbers the caption or the text may quote, so nothing is typed from memory.
cat(sprintf("partitions drawn: %s and %s (%s); y-axis %s; curves are %s\n", prim, prim2, scale_lab,
            if (yscale == "log") "log CEA" else "CEA (ng/mL)",
            if (curve == "mean") "means in 2-month bins" else "GAM smooths of the observations (geom_smooth default)"))
for (nm in names(curves)) for (k in seq_along(curves[[nm]])) {
  z <- curves[[nm]][[k]]
  imin  <- which.min(z$m)
  cross <- z$x[seq_along(z$x) > imin & z$m > ref][1]
  cat(sprintf("%-6s class %d  n = %4d  at 0: %.2f  min %.2f (month %4.1f)  max after month 6: %.2f (month %4.1f)  above %s after the minimum from month %s  at %d: %.2f\n",
              nm, k, sum((if (nm == "RKCFC") g1 else g2) == k),
              z$m[1], z$m[imin], z$x[imin],
              max(z$m[z$x > 6], na.rm = TRUE), z$x[z$x > 6][which.max(z$m[z$x > 6])],
              if (yscale == "log") "log(5)" else "5 ng/mL", if (is.na(cross)) "never" else sprintf("%.1f", cross),
              as.integer(Tmax), z$m[length(z$m)]))
}
cat(sprintf("wrote %s\n", out))
