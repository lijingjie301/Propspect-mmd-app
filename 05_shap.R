# 05 -- SHAP for the three models --------------------------------------------
# KernelSHAP on the predicted 12-month absolute risk. CSC and Fine-Gray are
# cheap; the forest uses a wrapper that returns the cause-1 CIF at 12 months.

source("R/utils.R")
suppressPackageStartupMessages({
  library(survival); library(prodlim); library(riskRegression); library(randomForestSRC)
  library(kernelshap); library(shapviz); library(ggplot2); library(patchwork)
})

DIR_OUT <- "output"
DIR_FIG <- "output/figures"
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

fit <- readRDS(file.path(DIR_OUT, "03_models.rds"))
models <- fit$models
train <- fit$train; valid <- fit$valid
vars <- fit$vars
labels <- c(CSC = "CSC", FGR = "Fine-Gray", RSF = "RSF")
H <- 12
B_BG <- 100
N_EXPLAIN <- 120

pred_fun <- function(model, X) {
  if (inherits(model, "rfsrc")) {
    pr <- predict(model, newdata = X)
    as.numeric(pr$cif[, which.min(abs(pr$time.interest - H)), 1])
  } else {
    as.numeric(predictRisk(model, newdata = X, times = H))
  }
}

# background sample from the training set, explained sample from validation
set.seed(11); bg <- train[sample(nrow(train), B_BG), vars]
set.seed(12); ex <- valid[sample(nrow(valid), min(N_EXPLAIN, nrow(valid))), vars]

shap_list <- list()
for (m in names(labels)) {
  t0 <- Sys.time()
  ks <- kernelshap(models[[m]], X = ex, bg_X = bg, pred_fun = pred_fun, verbose = FALSE)
  shap_list[[m]] <- ks
  cat(sprintf("%s done in %.0f s\n", labels[[m]], as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

name_map <- c(age = "Clinical: age", suzuki = "Clinical: Suzuki stage",
              pci = "Clinical: PCA involvement", symptom = "Clinical: symptom type",
              pinlv = "Clinical: symptom frequency", vwi = "Imaging: HR-VWI",
              TTP = "Imaging: TTP", chixushuxian = "Clinical: episode duration",
              chixushijian = "Clinical: episode duration", sex = "Clinical: sex",
              history = "Clinical: family history")

## beeswarm per model
beeswarm_one <- function(ks, title) {
  sv <- shapviz(ks$S, X = ks$X)
  p <- sv_importance(sv, kind = "beeswarm", max_display = length(vars), show_numbers = FALSE) +
    labs(x = sprintf("SHAP value (impact on %d-month risk)", H), y = NULL, title = title) +
    scale_colour_gradient(low = "#2C62A6", high = "#D81E5B", name = "Feature value",
                          breaks = c(0, 1), labels = c("Low", "High")) +
    theme_pub(7)
  p
}

ps <- lapply(names(labels), function(m) beeswarm_one(shap_list[[m]], labels[[m]]))
p_shap <- wrap_plots(ps, nrow = 1) + plot_layout(guides = "collect") &
  theme(legend.position = "right")
save_fig(p_shap, file.path(DIR_FIG, "fig8_shap_beeswarm"), 9.0, 3.4)

## mean |SHAP| bar chart
imp <- do.call(rbind, lapply(names(labels), function(m) {
  a <- colMeans(abs(shap_list[[m]]$S))
  data.frame(model = labels[[m]], variable = names(a), mean_abs = as.numeric(a))
}))
imp$variable <- factor(imp$variable, levels = names(sort(tapply(imp$mean_abs, imp$variable, mean))))
p_bar <- ggplot(imp, aes(mean_abs, variable, fill = model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  scale_fill_manual(values = c(CSC = "#01539E", `Fine-Gray` = "#73B748", RSF = "#EE2E22"), name = NULL) +
  labs(x = "mean |SHAP|", y = NULL) +
  theme_pub()
save_fig(p_bar, file.path(DIR_FIG, "fig9_shap_importance"), 4.8, 3.2)

write.csv(imp, file.path(DIR_OUT, "tables", "shap_importance.csv"), row.names = FALSE)
saveRDS(shap_list, file.path(DIR_OUT, "05_shap.rds"))
cat("shap finished\n")
