# Shared helpers -------------------------------------------------------------
# Date parsing, competing-risk utilities, IPCW metrics, plotting theme.
# Everything here is plain base R / survival; no tidyverse dependency.

## ---- date parsing ----------------------------------------------------------
# The source tables store dates either as yyyy-mm-dd strings, as Excel serials,
# or as Chinese free text ("未卒中", "未手术"). This keeps all three cases.
parse_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXct")) return(as.Date(x))
  s <- trimws(as.character(x))
  out <- as.Date(rep(NA_character_, length(s)))

  # dates come as 2025-9-16 as well as 2025-09-16, so let as.Date do the work
  iso <- grepl("^[0-9]{4}[-/][0-9]{1,2}[-/][0-9]{1,2}", s)
  out[iso] <- as.Date(s[iso], tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))

  num <- suppressWarnings(as.numeric(s))
  serial <- is.na(out) & !iso & !is.na(num)
  out[serial] <- as.Date(num[serial], origin = "1899-12-30")
  out
}

## ---- competing-risk outcome ------------------------------------------------
# status: 0 censored, 1 stroke, 2 revascularisation (competing event)
build_outcome <- function(dat, col_base, col_stroke, col_op, col_last) {
  b  <- parse_date(dat[[col_base]])
  ts <- as.numeric(parse_date(dat[[col_stroke]]) - b)
  to <- as.numeric(parse_date(dat[[col_op]])     - b)
  tl <- as.numeric(parse_date(dat[[col_last]])   - b)

  # First event wins: stroke if it happens before surgery, surgery otherwise.
  # Surgery is compared with the stroke date, not with the last follow-up, since
  # in this cohort the two coincide for operated patients.
  status <- integer(nrow(dat))
  time   <- tl

  i1 <- !is.na(ts)
  status[i1] <- 1L
  time[i1]   <- ts[i1]

  i2 <- !is.na(to) & (is.na(ts) | to < ts)
  status[i2] <- 2L
  time[i2]   <- to[i2]

  time[!is.na(time) & time < 0] <- NA

  data.frame(time = time / 30.4375, status = status)   # months
}

## ---- censoring distribution ------------------------------------------------
km_cens <- function(dat) survival::survfit(survival::Surv(dat$time, dat$status == 0) ~ 1)

gfun <- function(fit, t) {
  s <- summary(fit, times = t, extend = TRUE)$surv
  s[is.na(s)] <- 1
  pmax(s, 1e-8)
}

## ---- Aalen-Johansen CIF for cause 1 ---------------------------------------
# Used for observed incidence; hand-rolled so it also works on tiny subgroups
# where cuminc() drops components.
aj_cif <- function(dat, horizon) {
  tt <- sort(unique(dat$time[dat$status %in% c(1, 2)]))
  tt <- tt[tt <= horizon]
  if (!length(tt)) return(0)
  S <- 1; cif <- 0
  for (t in tt) {
    at_risk <- sum(dat$time >= t)
    d1 <- sum(dat$time == t & dat$status == 1)
    da <- sum(dat$time == t & dat$status %in% c(1, 2))
    if (at_risk <= 0) next
    cif <- cif + S * d1 / at_risk
    S   <- S * (1 - da / at_risk)
  }
  cif
}

## ---- IPCW metrics ----------------------------------------------------------
roc_ipcw <- function(risk, dat, horizon, cens_fit) {
  g <- gfun(cens_fit, pmax(dat$time - 1e-6, 1e-3))
  case <- dat$status == 1 & dat$time <= horizon
  ctrl <- dat$time > horizon
  wc <- ifelse(case, 1 / g, 0)
  wn <- as.numeric(ctrl)
  thr <- sort(unique(risk), decreasing = TRUE)
  data.frame(
    fpr = c(0, vapply(thr, function(c0) sum(wn * (risk >= c0)) / sum(wn), numeric(1)), 1),
    tpr = c(0, vapply(thr, function(c0) sum(wc * (risk >= c0)) / sum(wc), numeric(1)), 1)
  )
}

auc_ipcw <- function(risk, dat, horizon, cens_fit) {
  g <- gfun(cens_fit, pmax(dat$time - 1e-6, 1e-3))
  case <- which(dat$status == 1 & dat$time <= horizon)
  ctrl <- which(dat$time > horizon)
  if (!length(case) || !length(ctrl)) return(NA_real_)
  num <- den <- 0
  for (i in case) {
    wi <- 1 / g[i]
    num <- num + wi * sum(risk[ctrl] < risk[i])
    den <- den + wi * length(ctrl)
  }
  num / den
}

brier_ipcw <- function(dat, risk, horizon) {
  fit <- km_cens(dat)
  g_t <- gfun(fit, pmax(dat$time - 1e-6, 1e-3))
  g_h <- gfun(fit, rep(horizon, nrow(dat)))
  ev <- as.integer(dat$time <= horizon & dat$status == 1)
  w <- ifelse(ev == 1, 1 / g_t, ifelse(dat$time > horizon, 1 / g_h, 0))
  sum(w * (ev - risk)^2) / nrow(dat)
}

net_benefit <- function(dat, risk, horizon, thresholds) {
  fit <- km_cens(dat)
  g_t <- gfun(fit, pmax(dat$time - 1e-6, 1e-3))
  g_h <- gfun(fit, rep(horizon, nrow(dat)))
  ev <- as.integer(dat$time <= horizon & dat$status == 1)
  w <- ifelse(ev == 1, 1 / g_t, ifelse(dat$time > horizon, 1 / g_h, 0))
  n <- sum(w)
  vapply(thresholds, function(pt) {
    pos <- risk >= pt
    sum(w * pos * (ev == 1)) / n - sum(w * pos * (ev == 0)) / n * pt / (1 - pt)
  }, numeric(1))
}

## ---- calibration -----------------------------------------------------------
calibration_points <- function(dat, risk, horizon, groups = 10) {
  br <- unique(quantile(risk, probs = seq(0, 1, 1 / groups)))
  if (length(br) < 3) return(NULL)
  cut_g <- cut(risk, breaks = br, include.lowest = TRUE)
  obs <- tapply(seq_along(risk), cut_g, function(i) aj_cif(dat[i, ], horizon))
  data.frame(pred = as.numeric(tapply(risk, cut_g, mean)),
             obs  = as.numeric(obs),
             n    = as.numeric(table(cut_g)))
}

calibration_slope <- function(dat, risk, horizon) {
  ev <- ifelse(dat$time <= horizon & dat$status == 1, 1L,
               ifelse(dat$time > horizon, 0L, NA_integer_))
  ok <- !is.na(ev)
  if (sum(ev[ok]) < 5) return(c(intercept = NA, slope = NA))
  lp <- qlogis(pmin(pmax(risk[ok], 1e-6), 1 - 1e-6))
  fit <- glm(ev[ok] ~ lp, family = binomial)
  c(intercept = unname(coef(fit)[1]), slope = unname(coef(fit)[2]))
}

ici_deciles <- function(dat, risk, horizon, groups = 10) {
  cp <- calibration_points(dat, risk, horizon, groups)
  if (is.null(cp)) return(NA_real_)
  w <- cp$n / sum(cp$n)
  sum(w * abs(cp$obs - cp$pred))
}

## ---- plotting --------------------------------------------------------------
theme_pub <- function(base = 8) {
  ggplot2::theme_classic(base_size = base) +
    ggplot2::theme(
      axis.line  = ggplot2::element_line(linewidth = 0.45, colour = "grey20"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "grey20"),
      plot.title = ggplot2::element_text(size = base, hjust = 0.5),
      legend.position = "top"
    )
}

save_fig <- function(plot, file, width, height) {
  ggplot2::ggsave(paste0(file, ".png"), plot, width = width, height = height, dpi = 600)
  ggplot2::ggsave(paste0(file, ".pdf"), plot, width = width, height = height)
}
