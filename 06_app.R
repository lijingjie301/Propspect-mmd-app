# ==============================================================================
# PROPSPECT-MMD: competing-risks SHAP explorer
#
# Model : competing-risks random survival forest (cause = 1, stroke)
# Output: predicted 12-month cumulative incidence of stroke
# Data  : data/1250.csv  (development cohort, five variables)
#
# Run locally : open this file and source it, or Rscript shiny_app.R
# Deploy      : rsconnect::deployApp("path/to/this/folder")
# ==============================================================================

library(shiny); library(bslib); library(plotly)
library(survival); library(randomForestSRC); library(prodlim)
library(kernelshap); library(shapviz)
library(ggplot2); library(patchwork); library(dplyr)

set.seed(2024)

DATA_FILE  <- if (file.exists("data/1250.csv")) "data/1250.csv" else "data/example_1250.csv"
CACHE_FILE <- "shap_cache.rds"
HORIZON    <- 12                      # months
NTREE      <- 1000
CUTOFF     <- 0.237886                # fixed Youden-optimal 12-month cut-off (training set)
BG_N       <- 40                      # background sample for KernelSHAP
EXPLAIN_N  <- 60                      # patients used for the global SHAP plots

VARS <- c("suzuki", "pci", "symptom", "pinlv", "TTP")

# ---- 1. data ----------------------------------------------------------------
parse_date <- function(x) {
  s <- trimws(as.character(x))
  out <- as.Date(rep(NA_character_, length(s)))
  iso <- grepl("^[0-9]{4}[-/][0-9]{1,2}[-/][0-9]{1,2}", s)
  out[iso] <- as.Date(s[iso], tryFormats = c("%Y-%m-%d", "%Y/%m/%d"))
  num <- suppressWarnings(as.numeric(s))
  serial <- is.na(out) & !iso & !is.na(num)
  out[serial] <- as.Date(num[serial], origin = "1899-12-30")
  out
}

cat("loading data and fitting the model ...\n")
df <- read.csv(DATA_FILE, fileEncoding = "GB18030", check.names = FALSE,
               stringsAsFactors = FALSE)

b  <- parse_date(df[["基线日期"]])
ts <- as.numeric(parse_date(df[["卒中日期"]]) - b)
to <- as.numeric(parse_date(df[["手术日期"]]) - b)
tl <- as.numeric(parse_date(df[["末次随访"]]) - b)

status <- integer(nrow(df)); time <- tl
i1 <- !is.na(ts); status[i1] <- 1L; time[i1] <- ts[i1]
i2 <- !is.na(to) & (is.na(ts) | to < ts); status[i2] <- 2L; time[i2] <- to[i2]
time[!is.na(time) & time < 0] <- NA

df$time <- time / 30.4375
df$status <- status
df <- df[!is.na(df$time) & df$time > 0, c(VARS, "time", "status")]

# ---- 2. 7:3 split and model -------------------------------------------------
idx <- sort(unlist(lapply(split(seq_len(nrow(df)), df$status), function(z)
  sample(z, round(0.7 * length(z))))))
train_data <- df[idx, ]
valid_data <- df[-idx, ]
cat(sprintf("training %d / validation %d\n", nrow(train_data), nrow(valid_data)))

.m <- NULL
MODEL_FILE <- "model.rds"
PRECOMP_FILE <- "precomputed.rds"
rsf <- NULL
if (file.exists(PRECOMP_FILE)) {
  .m <- readRDS(PRECOMP_FILE)
  cat("using precomputed predictions\n")
} else if (file.exists(MODEL_FILE)) {
  .m <- readRDS(MODEL_FILE)
  rsf <- .m$rsf
  cat("using the pre-fitted model\n")
} else {
  set.seed(2024)
  rsf <- rfsrc(Surv(time, status) ~ ., data = train_data,
               ntree = NTREE, nodesize = 10, cause = 1, splitrule = "logrank",
               importance = "none")
}

get_rsf <- function() {
  if (is.null(rsf)) {
    cat("loading the forest on first use ...\n")
    rsf <<- readRDS(MODEL_FILE)$rsf
  }
  rsf
}

risk_at <- function(model, newdata, horizon) {
  pr <- predict(model, newdata = newdata[, VARS, drop = FALSE])
  j <- which.min(abs(pr$time.interest - horizon))
  as.numeric(pr$cif[, j, 1])
}

pred_fun <- function(model, X) {
  X <- as.data.frame(X)
  pr <- predict(model, newdata = X[, VARS, drop = FALSE])
  j <- which.min(abs(pr$time.interest - HORIZON))
  as.numeric(pr$cif[, j, 1])
}

if (!is.null(.m)) {
  train_risk <- .m$train_risk
  valid_risk <- .m$valid_risk
} else {
  train_risk <- risk_at(rsf, train_data, HORIZON)
  valid_risk <- risk_at(rsf, valid_data, HORIZON)
}
cutoff     <- CUTOFF
cutoff_pct <- cutoff * 100

col_low  <- "#04559f"
col_high <- "#bd4b3a"

# ---- 3. labels --------------------------------------------------------------
global_var_labels <- c(
  "suzuki"  = "Suzuki stage",
  "pci"     = "PCA involvement",
  "symptom" = "Symptom type",
  "pinlv"   = "Symptom frequency",
  "TTP"     = "TTP (relative, unitless)"
)

global_value_labels <- list(
  "suzuki"  = c("0" = "Stage 1", "1" = "Stage 2", "2" = "Stage 3",
                "3" = "Stage 4", "4" = "Stage 5", "5" = "Stage 6"),
  "pci"     = c("0" = "None", "1" = "Unilateral", "2" = "Bilateral"),
  "symptom" = c("0" = "TIA", "1" = "Cerebral infarction",
                "2" = "Haemorrhagic stroke", "3" = "Asymptomatic"),
  "pinlv"   = c("0" = "< 2 / month", "1" = "\u2265 2 / month")
)

# ---- 4. background sample and global SHAP -----------------------------------
if (!is.null(.m)) {
  bg <- .m$bg
  base_value <- .m$base_value
} else {
  set.seed(11)
  bg <- as.data.frame(train_data[sample(nrow(train_data), BG_N), VARS])
  base_value <- mean(pred_fun(rsf, bg))   # E[f(x)] over the background sample
}

if (file.exists(CACHE_FILE)) {
  shap_global <- readRDS(CACHE_FILE)
} else {
  set.seed(12)
  ex <- as.data.frame(valid_data[sample(nrow(valid_data), min(EXPLAIN_N, nrow(valid_data))), VARS])
  cat("computing the global KernelSHAP (a few minutes) ...\n")
  ks <- kernelshap(get_rsf(), X = ex, bg_X = bg, pred_fun = pred_fun, verbose = FALSE)
  shap_global <- ks
  try(saveRDS(ks, CACHE_FILE), silent = TRUE)
}

# ---- 5. cumulative-incidence curves ----------------------------------------
aj_curve <- function(dat, group, horizon = 60) {
  out <- lapply(levels(group), function(g) {
    d <- dat[group == g, ]
    tt <- sort(unique(d$time[d$status %in% c(1, 2)])); tt <- tt[tt <= horizon]
    S <- 1; cif <- 0
    res <- data.frame(time = 0, est = 0, group = g)
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

train_data$risk_group <- factor(ifelse(train_risk > cutoff, "High risk", "Low risk"),
                                levels = c("Low risk", "High risk"))
cif_df <- aj_curve(train_data, train_data$risk_group)
atrisk_times <- c(0, 12, 24, 36, 48, 60)
atrisk <- do.call(rbind, lapply(levels(train_data$risk_group), function(g) {
  d <- train_data[train_data$risk_group == g, ]
  data.frame(group = g, time = atrisk_times,
             n = sapply(atrisk_times, function(t) sum(d$time >= t)))
}))

# ---- 6. metrics -------------------------------------------------------------
conc <- concordance(Surv(time, status == 1) ~ valid_risk, data = valid_data)
cindex_valid <- as.numeric(conc$concordance)

brier_at <- function(dat, risk, horizon) {
  fit <- survfit(Surv(dat$time, dat$status == 0) ~ 1)
  g_t <- summary(fit, times = pmax(dat$time - 1e-6, 1e-3), extend = TRUE)$surv
  g_t[is.na(g_t)] <- 1
  g_h <- summary(fit, times = rep(horizon, nrow(dat)), extend = TRUE)$surv
  g_h[is.na(g_h)] <- 1
  ev <- as.integer(dat$time <= horizon & dat$status == 1)
  w <- ifelse(ev == 1, 1 / g_t, ifelse(dat$time > horizon, 1 / g_h, 0))
  sum(w * (ev - risk)^2) / nrow(dat)
}
brier_valid <- brier_at(valid_data, valid_risk, HORIZON)

# ---- 7. UI ------------------------------------------------------------------
ui <- page_sidebar(
  title = "PROPSPECT-MMD: stroke risk and SHAP explanation",
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  tags$head(tags$style(HTML("
    .card-header { font-weight: 600; font-size: 16px; background-color: #f8f9fa; }
    .risk-text { font-size: 48px; font-weight: 800; margin: 0; }
    .level-text { font-size: 40px; font-weight: 800; margin: 0; padding-top: 6px; }
    .low-risk  { color: #04559f; }
    .high-risk { color: #bd4b3a; }
    .footer-text { font-size: 12px; color: #999; text-align: center; padding: 10px; }
  "))),

  sidebar = sidebar(
    title = "Patient characteristics", width = 360,
    sliderInput("TTP", "TTP (relative, unitless)", min = -2, max = 9, value = 3.0, step = 0.01),
    hr(),
    selectInput("suzuki", "Suzuki stage",
                choices = setNames(0:5, paste("Stage", 1:6)), selected = 3),
    selectInput("pci", "PCA involvement",
                choices = c("None" = 0, "Unilateral" = 1, "Bilateral" = 2)),
    selectInput("symptom", "Symptom type",
                choices = c("TIA" = 0, "Cerebral infarction" = 1,
                            "Haemorrhagic stroke" = 2, "Asymptomatic" = 3)),
    selectInput("pinlv", "Symptom frequency",
                choices = c("< 2 per month" = 0, "\u2265 2 per month" = 1)),
    hr(),
    p(class = "text-muted small",
      sprintf("Cut-off: %.1f%% (%d-month risk) | %d variables | %d training patients",
              cutoff_pct, HORIZON, length(VARS), nrow(train_data)))
  ),

  layout_columns(
    col_widths = c(4, 4, 4),
    card(card_header(span(icon("percent"), " Predicted risk")),
         div(style = "text-align:center; padding:10px;", uiOutput("prob_text"))),
    card(card_header(span(icon("circle-info"), " Risk group")),
         div(style = "text-align:center; padding:10px;", uiOutput("level_text"))),
    card(card_header(span(icon("clock"), " Predicted incidence")),
         div(style = "padding:10px;", uiOutput("surv_text")))
  ),

  layout_columns(
    col_widths = c(12, 6, 6),
    card(card_header(span(icon("chart-line"), " Risk gauge")),
         plotlyOutput("gauge_plot", height = "260px")),
    card(card_header(span(icon("wand-magic-sparkles"), " SHAP waterfall")),
         plotOutput("dynamic_shap_waterfall", height = "480px")),
    card(card_header(span(icon("chart-area"), " Cumulative incidence by risk group")),
         plotOutput("cif_curve", height = "480px"))
  ),

  card(card_header(span(icon("chart-bar"), " Global feature importance (mean |SHAP|)")),
       plotOutput("shap_importance", height = "380px")),

  div(class = "footer-text",
      sprintf("PROPSPECT-MMD | competing-risks random survival forest | cut-off %.2f%% | %d variables | %d training patients | point estimates only (no individual confidence interval)",
              cutoff_pct, length(VARS), nrow(train_data)))
)

# ---- 8. server --------------------------------------------------------------
server <- function(input, output, session) {

  get_current_data <- reactive({
    data.frame(suzuki = as.numeric(input$suzuki),
               pci = as.numeric(input$pci),
               symptom = as.numeric(input$symptom),
               pinlv = as.numeric(input$pinlv),
               TTP = as.numeric(input$TTP))[, VARS, drop = FALSE]
  })

  current_risk_prob <- reactive({
    p <- risk_at(get_rsf(), get_current_data(), HORIZON) * 100
    min(max(p, 0), 100)
  })

  output$prob_text <- renderUI({
    prob <- current_risk_prob()
    cls <- if (prob <= cutoff_pct) "low-risk" else "high-risk"
    HTML(sprintf('<p class="risk-text %s">%.2f%%</p>', cls, prob))
  })

  output$level_text <- renderUI({
    prob <- current_risk_prob()
    if (prob <= cutoff_pct) {
      HTML('<p class="level-text low-risk">&#10003; Low risk</p>')
    } else {
      HTML('<p class="level-text high-risk">&#9888; High risk</p>')
    }
  })

  output$surv_text <- renderUI({
    prob <- current_risk_prob()
    grp <- if (prob > cutoff_pct) "High risk" else "Low risk"
    d <- train_data[train_data$risk_group == grp, ]
    c12 <- aj_curve(d, factor(rep(grp, nrow(d)), levels = grp), 12)
    last <- function(x) round(100 * x$est[nrow(x)], 1)
    HTML(sprintf(
      '<div style="font-size:14px;">
        <b>%s</b><br>
        Predicted 12-month risk (point estimate): <b>%.2f%%</b><br>
        Observed 12-month cumulative incidence in this risk group: <b>%.1f%%</b><br>
        <span style="color:#777;">Point estimate only; no individual confidence interval is provided.</span>
      </div>', grp, prob, last(c12)))
  })

  output$gauge_plot <- renderPlotly({
    prob <- current_risk_prob()
    plot_ly(
      type = "indicator", mode = "gauge+number", value = prob,
      number = list(suffix = "%", valueformat = ".2f"),
      gauge = list(
        axis = list(range = list(0, 100)),
        bar = list(color = "rgba(0,0,0,0.45)", thickness = 0.2),
        bgcolor = "white", borderwidth = 2, bordercolor = "gray",
        steps = list(
          list(range = c(0, cutoff_pct), color = col_low),
          list(range = c(cutoff_pct, 100), color = col_high)
        )
      )) %>% layout(margin = list(l = 20, r = 20, t = 30, b = 20))
  })

  output$dynamic_shap_waterfall <- renderPlot({
    X_new <- get_current_data()
    prob <- current_risk_prob()

    ks <- kernelshap(get_rsf(), X = X_new, bg_X = bg, pred_fun = pred_fun, verbose = FALSE)
    sv <- shapviz(ks$S, X = ks$X, baseline = base_value)
    sv$X <- as.data.frame(sv$X)
    sv$S <- as.matrix(sv$S)
    colnames(sv$S) <- global_var_labels[colnames(sv$S)]

    for (v in names(global_value_labels)) {
      if (v %in% colnames(sv$X)) {
        sv$X[[v]] <- global_value_labels[[v]][as.character(sv$X[[v]])]
      }
    }
    colnames(sv$X) <- global_var_labels[colnames(sv$X)]

    tcol <- if (prob > cutoff_pct) col_high else col_low
    sv_waterfall(sv, fill_colors = c("#fb3076", "#258dde")) +
      theme_bw(base_size = 13) +
      labs(title = "Individual risk-factor contributions",
           subtitle = sprintf("Predicted %d-month risk: %.2f%% (%s)",
                              HORIZON, prob, if (prob > cutoff_pct) "high risk" else "low risk")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
            plot.subtitle = element_text(hjust = 0.5, colour = tcol, face = "bold", size = 13),
            axis.text.y = element_text(size = 11, face = "bold", colour = "black"),
            axis.text.x = element_text(size = 11, colour = "black"),
            panel.grid.minor = element_blank())
  })

  output$cif_curve <- renderPlot({
    prob <- current_risk_prob()
    pt_grp <- if (prob > cutoff_pct) "High risk" else "Low risk"

    p1 <- ggplot(cif_df, aes(time, est, colour = group)) +
      geom_step(linewidth = 0.7) +
      scale_colour_manual(values = c(`Low risk` = col_low, `High risk` = col_high), name = NULL) +
      scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 12), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, NA), labels = function(x) sprintf("%.0f%%", 100 * x)) +
      labs(x = "Time (months)", y = "Cumulative incidence of stroke") +
      annotate("text", x = 2, y = 0.02, hjust = 0, size = 4.2, fontface = "bold",
               colour = if (pt_grp == "High risk") col_high else col_low,
               label = paste0("This patient: ", pt_grp)) +
      theme_classic(base_size = 12) +
      theme(legend.position = c(0.75, 0.12), legend.background = element_blank())

    risk_tbl <- ggplot(atrisk, aes(time, group, label = n)) +
      geom_text(size = 3.2) +
      scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 12), expand = c(0, 0)) +
      labs(x = NULL, y = NULL, title = "Number at risk") +
      theme_classic(base_size = 10) +
      theme(axis.line = element_blank(), axis.ticks = element_blank(),
            axis.text.x = element_blank(), plot.title = element_text(size = 9))

    p1 / risk_tbl + plot_layout(heights = c(1, 0.22))
  })

  output$shap_importance <- renderPlot({
    sv <- shapviz(shap_global$S, X = shap_global$X, baseline = base_value)
    colnames(sv$S) <- global_var_labels[colnames(sv$S)]
    colnames(sv$X) <- global_var_labels[colnames(sv$X)]
    sv_importance(sv, kind = "beeswarm", show_numbers = FALSE) +
      scale_colour_gradient(low = "#2C62A6", high = "#D81E5B", name = "Feature value",
                            breaks = c(0, 1), labels = c("Low", "High")) +
      labs(x = sprintf("SHAP value (impact on %d-month risk)", HORIZON), y = NULL) +
      theme_bw(base_size = 12) +
      theme(legend.position = "right", axis.text.y = element_text(size = 10, face = "bold"))
  })


}

cat(sprintf("\napp ready | cut-off %.2f%% | %d training patients\n", cutoff_pct, nrow(train_data)))
shinyApp(ui, server)
