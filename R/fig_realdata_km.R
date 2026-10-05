# FIGURE 2 -- Kaplan-Meier survival by class, with numbers at risk.
#
#   Rscript R/fig_realdata_km.R [f1dir] [outfile]
#
# THE TIME ORIGIN. For a 36-month record it is REAL_COX_ORIGIN of
# R/realdata36_settings.R: "surgery" since 2026-10-03 (the
# reported estimate is the from-surgery Cox model, with the month-36 landmark as the
# sensitivity analysis), "landmark" before. For the 18-month record it is the
# landmark, as published. The landmark version plots patients alive and under
# follow-up at the landmark from that point, with the axis labelled in months after
# surgery so the offset is visible; the from-surgery version plots every patient
# from surgery. Membership is defined over the clustering window, so a from-surgery
# curve carries guarantee-time bias, which is why both are reported; the caption
# names the origin. Numbers at risk are printed beneath.
#
# The log-rank P in each panel is computed here from the same cohort F2 uses, so it
# matches F2's output by construction rather than by transcription.
#
# PANEL TITLES name the method only ("(a) RKCFC", "(b) KCFC"); the caption gives
# the scale and the classes. The x-axis is labeled every 24 months from surgery,
# which also keeps the numbers at risk from running into each other (with 12-month
# ticks the low-CEA row overlapped).
#
# NO IDENTIFIERS: survival time and event come from the analysis set that
# load_cea36() returns; nothing patient-level is written.

source(file.path("R", "realdata_load.R"))
if (!requireNamespace("survival", quietly = TRUE)) stop("package `survival` is required.", call. = FALSE)

args  <- commandArgs(trailingOnly = TRUE)
args  <- args[args != "--args"]
f1dir <- if (length(args) >= 1) args[1] else file.path("results_real", "F1")
out   <- if (length(args) >= 2) args[2] else file.path("..", "doc", "V2", "figures", "fig2_km_landmark.png")

f1 <- readRDS(file.path(f1dir, "F1.rds"))
if (isTRUE(any(f1$table$smoke))) stop("F1.rds is a smoke run.", call. = FALSE)
Tmax <- if (!is.null(f1$Tmax)) f1$Tmax else REAL_TMAX
LM   <- if (!is.null(f1$landmark)) f1$landmark else REAL_LANDMARK
ORIGIN <- "landmark"
if (isTRUE(Tmax == 36)) { source(file.path("R", "realdata36_settings.R")); ORIGIN <- REAL_COX_ORIGIN }
prim  <- if (!is.null(f1$primary)) f1$primary else "log_M1"
prim2 <- sub("M1$", "M2", prim)
d  <- truncate_cea36(load_cea36(), Tmax); b <- d$base
if (ORIGIN == "landmark") {
  at <- b$Survival_Time >= LM; off0 <- LM
  xlab <- sprintf("Months after surgery, origin month %d", LM)
} else {
  at <- rep(TRUE, nrow(b)); off0 <- 0
  xlab <- "Months after surgery"
}
tt <- b$Survival_Time[at] - off0; ev <- b$Death[at]
ticks <- seq(0, 24 * ceiling(max(tt) / 24), by = 24)
cols  <- c("#c0392b", "#2e86c1")

panel <- function(g, labels, title) {
  gg  <- factor(g[at])
  fit <- survival::survfit(survival::Surv(tt, ev) ~ gg)
  sdf <- survival::survdiff(survival::Surv(tt, ev) ~ gg)
  p   <- stats::pchisq(sdf$chisq, df = length(sdf$n) - 1, lower.tail = FALSE)
  plot(fit, col = cols, lwd = 2.2, conf.int = FALSE, mark.time = TRUE, xaxt = "n", las = 1,
       xlim = c(0, max(ticks)), ylim = c(0.5, 1), xlab = "", ylab = "Survival probability", main = title)
  axis(1, at = ticks, labels = ticks + off0)
  title(xlab = xlab, line = 2.3)
  legend("bottomleft", legend = sprintf("%s (n = %d)", labels, as.integer(table(gg))),
         col = cols, lwd = 2.2, bty = "n", cex = 0.85, inset = c(0, 0.08))
  text(max(ticks), 0.98, if (p < 0.001) "log-rank P < 0.001" else sprintf("log-rank P = %s", formatC(p, format = "f", digits = 3)),
       adj = 1, cex = 0.85)
  nr <- summary(fit, times = ticks, extend = TRUE)
  off <- -0.07 * max(ticks)
  mtext("At risk", side = 1, line = 3.4, at = off, adj = 1, cex = 0.65)
  for (k in seq_along(levels(gg))) {
    v <- nr$n.risk[as.integer(nr$strata) == k]
    mtext(v, side = 1, line = 3.4 + 0.9 * k, at = ticks, cex = 0.65, col = cols[k])
    mtext(labels[k], side = 1, line = 3.4 + 0.9 * k, at = off, adj = 1, cex = 0.65, col = cols[k])
  }
  cat(sprintf("%s: at origin %s, events %s, log-rank chisq %.2f p %.4g\n", title,
              paste(as.integer(table(gg)), collapse = "/"),
              paste(as.integer(tapply(ev, gg, sum)), collapse = "/"), sdf$chisq, p))
}

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
grDevices::png(out, width = 2400, height = 1250, res = 300)
graphics::par(mfrow = c(1, 2), mar = c(7.2, 4.2, 2.5, 1))
cn <- if (!is.null(f1$class_names)) f1$class_names else c("rebounding", "flat")
cat(sprintf("origin: %s; partitions %s and %s\n", ORIGIN, prim, prim2))
panel(f1$fits[[prim]],  cn, "(a) RKCFC")
# KCFC's classes have no names; the text calls them the smaller and the larger class,
# so the legend and the at-risk rows do the same.
size_lab <- function(g) { n <- as.integer(table(g)); if (n[1] <= n[2]) c("smaller class", "larger class") else c("larger class", "smaller class") }
panel(f1$fits[[prim2]], size_lab(f1$fits[[prim2]]), "(b) KCFC")
invisible(grDevices::dev.off())
cat(sprintf("wrote %s\n", out))
