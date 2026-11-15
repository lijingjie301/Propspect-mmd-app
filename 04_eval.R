# 04 -- Evaluation: AUC, ROC, calibration, DCA, risk strata -------------------
# Produces the performance tables (csv), the discrimination/calibration/DCA
# figures and the two-group cumulative-incidence curves.

source("R/utils.R")
suppressPackageStartupMessages({
  library(survival); library(prodlim); library(riskRegression)
  library(cmprsk); library(ggplot2); library(patchwork)
})

DIR_OUT <- "output"
DIR_FIG <- "output/figures"
DIR_TAB <- "output/tables"
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_TAB, showWarnings = FALSE, recursive = TRUE)

fit <- readRDS(file.path(DIR_OUT, "03_models.rds"))
valid <- fit$valid; ext <- fit$external; train <- fit$train
horizons <- fit$horizons; vars <- fit$vars
labels <- c(CSC = "CSC", FGR = "Fine-Gray", RSF = "RSF")
cols <- c(CSC = "#01539E", `Fine-Gray` = "#73B748", RSF = "#EE2E22")
cohorts <- list(`Internal validation` = valid, `External test` = ext)
B <- 400

## ---- discrimination and calibration table ---------------------------------
rows <- list()
for (h in horizons) {
  for (m in names(labels)) {
    for (cn in names(cohorts)) {
      dat <- cohorts[[cn]]
      r <- fit$risk[[if (cn == "Internal validation") "valid" else "external"]][[m]][, paste0("t", h)]
      cf <- km_cens(dat)
      a <- auc_ipcw(r, dat, h, cf)
      set.seed(11)
      bs <- replicate(B, {
        i <- sample(length(r), replace = TRUE)
        auc_ipcw(r[i], dat[i, ], h, km_cens(dat[i, ]))
      })
      cal <- calibration_slope(dat, r, h)
      cp <- calibration_points(dat, r, h)
      oe <- aj_cif(dat, h) / mean(r)
      rows[[length(rows) + 1]] <- data.frame(
        horizon = h, model = labels[[m]], cohort = cn,
        auc = a, auc_lo = quantile(bs, .025, na.rm = TRUE), auc_hi = quantile(bs, .975, na.rm = TRUE),
        oe = oe, slope = cal[["slope"]], intercept = cal[["intercept"]],
        brier = brier_ipcw(dat, r, h), ici = ici_deciles(dat, r, h)
      )
    }
  }
}
perf <- do.call(rbind, rows)
write.csv(perf, file.path(DIR_TAB, "performance.csv"), row.names = FALSE)
print(perf[, c("horizon", "model", "cohort", "auc", "oe", "brier")], digits = 3)

## ---- ROC panels ------------------------------------------------------------
roc_panel <- function(dat, h, show_y) {
  cf <- km_cens(dat)
  d_all <- do.call(rbind, lapply(names(labels), function(m) {
    r <- fit$risk[[if (nrow(dat) == nrow(valid)) "valid" else "external"]][[m]][, paste0("t", h)]
    cbind(model = labels[[m]], roc_ipcw(r, dat, h, cf))
  }))
  ggplot(d_all, aes(fpr, tpr, colour = model)) +
    geom_abline(linetype = "dashed", colour = "grey65", linewidth = 0.3) +
    geom_line(linewidth = 0.55) +
    scale_colour_manual(values = cols, name = NULL, limits = names(cols)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, .25)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, .25)) +
    labs(x = "1 - Specificity", y = if (show_y) "Sensitivity" else NULL,
         title = sprintf("%d months", h)) +
    theme_pub()
}

ps <- list()
for (cn in names(cohorts)) for (h in horizons) {
  ps[[length(ps) + 1]] <- roc_panel(cohorts[[cn]], h, h == horizons[1])
}
p_roc <- wrap_plots(ps, nrow = 2) + plot_layout(guides = "collect") &
  theme(legend.position = "top")
save_fig(p_roc, file.path(DIR_FIG, "fig4_roc"), 8.4, 4.6)

## ---- calibration panels ----------------------------------------------------
cal_panel <- function(dat, h, m) {
  r <- fit$risk[[if (nrow(dat) == nrow(valid)) "valid" else "external"]][[m]][, paste0("t", h)]
  cp <- calibration_points(dat, r, h)
  lim <- c(0, max(c(cp$pred, cp$obs)) * 1.1)
  ggplot(cp, aes(pred, obs)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
    geom_point(aes(colour = obs >= pred), size = 1.6) +
    scale_colour_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#1F6FB4"), guide = "none") +
    coord_equal(xlim = lim, ylim = lim) +
    labs(x = "Predicted", y = "Observed", title = sprintf("%s, %d months", labels[[m]], h)) +
    theme_pub(7)
}
ps <- list()
for (cn in names(cohorts)) for (m in names(labels)) {
  ps[[length(ps) + 1]] <- cal_panel(cohorts[[cn]], 12, m)
}
p_cal <- wrap_plots(ps, nrow = 2)
save_fig(p_cal, file.path(DIR_FIG, "fig5_calibration_12m"), 7.2, 4.6)

## ---- DCA -------------------------------------------------------------------
th <- seq(0.02, 0.60, by = 0.01)
dca_rows <- list()
for (h in c(12, 24)) {
  for (cn in names(cohorts)) {
    dat <- cohorts[[cn]]
    ev <- dat$status
    for (m in names(labels)) {
      r <- fit$risk[[if (cn == "Internal validation") "valid" else "external"]][[m]][, paste0("t", h)]
      dca_rows[[length(dca_rows) + 1]] <- data.frame(
        horizon = h, cohort = cn, strategy = labels[[m]], threshold = th,
        nb = net_benefit(dat, r, h, th))
    }
    # reference strategies: treat everyone and treat nobody
    cf <- km_cens(dat)
    g_t <- gfun(cf, pmax(dat$time - 1e-6, 1e-3))
    g_h <- gfun(cf, rep(h, nrow(dat)))
    ev <- as.integer(dat$time <= h & dat$status == 1)
    w <- ifelse(ev == 1, 1 / g_t, ifelse(dat$time > h, 1 / g_h, 0))
    prev <- sum(w * ev) / sum(w)
    dca_rows[[length(dca_rows) + 1]] <- data.frame(
      horizon = h, cohort = cn, strategy = "Treat all", threshold = th,
      nb = prev - (1 - prev) * th / (1 - th))
    dca_rows[[length(dca_rows) + 1]] <- data.frame(
      horizon = h, cohort = cn, strategy = "Treat none", threshold = th, nb = 0)
  }
}
dca <- do.call(rbind, dca_rows)
write.csv(dca, file.path(DIR_TAB, "dca.csv"), row.names = FALSE)

p_dca <- ggplot(subset(dca, horizon == 12), aes(threshold, nb, colour = strategy)) +
  geom_line(linewidth = 0.6) +
  geom_hline(yintercept = 0, colour = "grey20", linewidth = 0.3) +
  facet_wrap(~cohort, nrow = 1) +
  scale_colour_manual(values = c(cols, `Treat all` = "grey55"), name = NULL) +
  labs(x = "Threshold probability", y = "Net benefit") +
  theme_pub()
save_fig(p_dca, file.path(DIR_FIG, "fig6_dca_12m"), 6.6, 3.2)

## ---- two-group risk stratification -----------------------------------------
# Cut-off: the 80th percentile of the training predicted risk at 24 months.
strat_rows <- list()
for (m in names(labels)) {
  r_tr <- fit$risk$train[[m]][, "t24"]
  cut <- quantile(r_tr, 0.8)
  for (cn in names(cohorts)) {
    dat <- cohorts[[cn]]
    r <- fit$risk[[if (cn == "Internal validation") "valid" else "external"]][[m]][, "t24"]
    g <- factor(ifelse(r >= cut, "High", "Low"), levels = c("Low", "High"))
    fg <- FGR(Hist(time, status) ~ g, data = dat, cause = 1)
    co <- summary(fg)$coef
    ev <- ifelse(dat$status == 1, 1, ifelse(dat$status == 2, 2, 0))
    gp <- tryCatch(as.numeric(cuminc(dat$time, ev, group = g)$Tests[1, 2]),
                   error = function(e0) NA_real_)
    strat_rows[[length(strat_rows) + 1]] <- data.frame(
      model = labels[[m]], cohort = cn, cut = cut,
      n_high = sum(g == "High"), n_low = sum(g == "Low"),
      shr = exp(co[1, 1]), lo = exp(co[1, 1] - 1.96 * co[1, 3]),
      hi = exp(co[1, 1] + 1.96 * co[1, 3]), p = co[1, 5], gray = gp,
      cif12_high = aj_cif(dat[g == "High", ], 12), cif12_low = aj_cif(dat[g == "Low", ], 12),
      cif24_high = aj_cif(dat[g == "High", ], 24), cif24_low = aj_cif(dat[g == "Low", ], 24))
  }
}
strat <- do.call(rbind, strat_rows)
write.csv(strat, file.path(DIR_TAB, "risk_strata.csv"), row.names = FALSE)
print(strat[, c("model", "cohort", "shr", "gray")], digits = 3)

## cumulative-incidence curves by group (24 months, both cohorts)
cif_curve <- function(dat, group, horizon) {
  sub <- split(dat, group)
  out <- lapply(names(sub), function(g) {
    d <- sub[[g]]
    tt <- sort(unique(d$time[d$status %in% c(1, 2)])); tt <- tt[tt <= horizon]
    S <- 1; cif <- 0; res <- data.frame(time = 0, est = 0, group = g)
    for (t in tt) {
      at <- sum(d$time >= t); d1 <- sum(d$time == t & d$status == 1)
      da <- sum(d$time == t & d$status %in% c(1, 2))
      if (at <= 0) next
      cif <- cif + S * d1 / at; S <- S * (1 - da / at)
      res <- rbind(res, data.frame(time = t, est = cif, group = g))
    }
    res
  })
  do.call(rbind, out)
}

ps <- list()
for (cn in names(cohorts)) {
  dat <- cohorts[[cn]]
  cut <- quantile(fit$risk$train$RSF[, "t24"], 0.8)
  g <- factor(ifelse(fit$risk[[if (cn == "Internal validation") "valid" else "external"]]$RSF[, "t24"] >= cut,
                     "High risk", "Low risk"), levels = c("Low risk", "High risk"))
  cd <- cif_curve(dat, g, 24)
  ps[[length(ps) + 1]] <- ggplot(cd, aes(time, est, colour = group)) +
    geom_step(linewidth = 0.6) +
    scale_colour_manual(values = c(`Low risk` = "#01539E", `High risk` = "#EE2E22"), name = NULL) +
    scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 3), expand = c(0, 0)) +
    labs(x = "Follow-up (months)", y = "Cumulative incidence", title = cn) +
    theme_pub()
}
p_cif <- wrap_plots(ps, nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "top")
save_fig(p_cif, file.path(DIR_FIG, "fig7_cif_riskgroups"), 6.4, 3.0)

cat("evaluation finished\n")
