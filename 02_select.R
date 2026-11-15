# 02 -- Variable selection: LASSO and Boruta, then the intersection ----------
# LASSO uses a Cox penalty; Boruta uses random-survival-forest importance
# evaluated on the cause-specific hazard. Only variables confirmed by both are
# carried forward. Writes output/02_selection.rds and the Venn figure.

source("R/utils.R")
suppressPackageStartupMessages({
  library(survival); library(glmnet); library(Boruta)
  library(randomForestSRC); library(ggplot2)
})

DIR_OUT <- "output"
DIR_FIG <- "output/figures"
dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)

prep <- readRDS(file.path(DIR_OUT, "01_prepared.rds"))
train <- prep$train; vars <- prep$vars; seed <- prep$seed

x <- as.matrix(train[, vars])
y <- Surv(train$time, as.integer(train$status == 1))

## ---- LASSO ----------------------------------------------------------------
set.seed(seed)
cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10, standardize = TRUE)

pick <- function(fit, lambda) {
  cf <- as.matrix(coef(fit, s = lambda))
  rownames(cf)[abs(cf[, 1]) > 1e-8]
}
lasso_min <- pick(cv_fit, cv_fit$lambda.min)
lasso_1se <- pick(cv_fit, cv_fit$lambda.1se)
cat("LASSO lambda.min:", length(lasso_min), "variables\n")
print(lasso_min)

## ---- Boruta ---------------------------------------------------------------
# getImp returns the cause-specific VIMP of a random survival forest, which is
# the survival analogue of Boruta's usual random-forest importance.
imp_rsf <- function(x, y) {
  ym <- as.matrix(y)
  d <- data.frame(x, time = ym[, "time"], status = ym[, "status"])
  fit <- rfsrc(Surv(time, status) ~ ., data = d, ntree = 300,
               nodesize = 10, cause = 1, splitrule = "logrank",
               importance = "permute")
  v <- if (is.matrix(fit$importance)) fit$importance[, 1] else fit$importance
  v <- as.numeric(v); names(v) <- colnames(x)
  v[!is.finite(v)] <- 0
  v
}

set.seed(seed)
bor <- Boruta(x = as.data.frame(x), y = y, getImp = imp_rsf,
              maxRuns = 100, doTrace = 0)
bor_stats <- attStats(bor)
bor_stats$variable <- rownames(bor_stats)
bor_stats <- bor_stats[order(-bor_stats$meanImp), ]
bor_confirmed <- bor_stats$variable[bor_stats$decision == "Confirmed"]
cat("Boruta confirmed:", length(bor_confirmed), "variables\n")
print(bor_confirmed)

## ---- intersection ----------------------------------------------------------
final_vars <- intersect(lasso_min, bor_confirmed)
cat("intersection:", length(final_vars), "->", paste(final_vars, collapse = ", "), "\n")

## ---- stability check -------------------------------------------------------
# Bootstrap the LASSO to see how often each variable is retained. This is
# reported for transparency; it does not drive the final set by itself.
set.seed(seed)
B <- 200
freq <- setNames(rep(0, length(vars)), vars)
for (b in seq_len(B)) {
  i <- sample(nrow(x), replace = TRUE)
  fit_b <- glmnet(x[i, ], y[i, ], family = "cox", alpha = 1,
                  standardize = TRUE, lambda = cv_fit$lambda.min)
  cf <- as.matrix(coef(fit_b))
  nz <- rownames(cf)[abs(cf[, 1]) > 1e-8]
  freq[nz] <- freq[nz] + 1
}
freq <- freq / B
print(round(sort(freq, decreasing = TRUE), 2))

## ---- final set -------------------------------------------------------------
# The intersection is what the two methods agree on. One variable is dropped
# afterwards: episode duration was retained by both methods but adds nothing to
# discrimination (see the sensitivity check in 04_eval.R) and largely duplicates
# the symptom-frequency item, so the more parsimonious seven-variable set was
# used in the manuscript. Set DROP_FINAL to character(0) to keep all of them.
DROP_FINAL <- c("chixushijian")

final_vars <- setdiff(intersect(lasso_min, bor_confirmed), DROP_FINAL)
cat("final set:", length(final_vars), "->", paste(final_vars, collapse = ", "), "\n")

## ---- figures ---------------------------------------------------------------
# Fig 1: LASSO cross-validation and coefficient path
xlog <- log(cv_fit$lambda)
cv_df <- data.frame(loglam = xlog, cvm = cv_fit$cvm,
                    lo = cv_fit$cvm - cv_fit$cvsd, hi = cv_fit$cvm + cv_fit$cvsd,
                    nz = cv_fit$nzero)
p_cv <- ggplot(cv_df, aes(loglam, cvm)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.03, linewidth = 0.25, colour = "grey75") +
  geom_point(colour = "#C0392B", size = 0.9) +
  geom_vline(xintercept = log(cv_fit$lambda.min), linetype = "dotted", linewidth = 0.4) +
  geom_vline(xintercept = log(cv_fit$lambda.1se), linetype = "dotted", linewidth = 0.4) +
  labs(x = expression(-log(lambda)), y = "Partial likelihood deviance") +
  theme_pub()

beta <- as.matrix(cv_fit$glmnet.fit$beta)      # variables x lambdas
path_long <- data.frame(
  loglam   = rep(xlog, each = nrow(beta)),
  variable = rep(rownames(beta), times = length(xlog)),
  coef     = as.numeric(beta)
)
path_long$selected <- path_long$variable %in% final_vars

p_path <- ggplot(path_long, aes(loglam, coef, group = variable, colour = selected,
                                linewidth = selected)) +
  geom_line() +
  scale_colour_manual(values = c(`TRUE` = "#2E6DA4", `FALSE` = "grey85"), guide = "none") +
  scale_linewidth_manual(values = c(`TRUE` = 0.55, `FALSE` = 0.25), guide = "none") +
  labs(x = expression(-log(lambda)), y = "Coefficient") +
  theme_pub()

save_fig(p_cv,   file.path(DIR_FIG, "fig1a_lasso_cv"),   4.2, 3.2)
save_fig(p_path, file.path(DIR_FIG, "fig1b_lasso_path"), 4.2, 3.2)

# Fig 2: Boruta importance
bor_df <- bor_stats[order(bor_stats$meanImp), ]
bor_df$variable <- factor(bor_df$variable, levels = bor_df$variable)
p_bor <- ggplot(bor_df, aes(meanImp, variable, fill = decision)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = c(Confirmed = "#C0392B", Rejected = "#7FB2D0",
                               Tentative = "#F5DEA8"), name = NULL) +
  labs(x = "Importance (Z-score)", y = NULL) +
  theme_pub()
save_fig(p_bor, file.path(DIR_FIG, "fig2_boruta"), 4.4, 4.2)

# Fig 3: Venn of the two selections
circ <- function(cx, cy, r, n = 400) {
  th <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(th), y = cy + r * sin(th))
}
n_l <- length(setdiff(lasso_min, bor_confirmed))
n_b <- length(setdiff(bor_confirmed, lasso_min))
n_i <- length(final_vars)
r1 <- sqrt(length(lasso_min)) * 0.35
r2 <- sqrt(length(bor_confirmed)) * 0.35
dd <- (r1 + r2) * 0.62
c1 <- circ(-dd / 2, 0, r1); c2 <- circ(dd / 2, 0, r2)
p_venn <- ggplot() +
  geom_polygon(data = c1, aes(x, y), fill = "#C0392B", alpha = 0.35, colour = "#8E2A1C") +
  geom_polygon(data = c2, aes(x, y), fill = "#4E7C72", alpha = 0.35, colour = "#4E7C72") +
  annotate("text", x = -dd / 2 - r1 * 0.45, y = 0, label = n_l, size = 3) +
  annotate("text", x = 0, y = 0, label = n_i, size = 3.2, fontface = "bold") +
  annotate("text", x = dd / 2 + r2 * 0.45, y = 0, label = n_b, size = 3) +
  annotate("text", x = -dd / 2, y = r1 * 0.85, label = "LASSO", size = 2.8, fontface = "bold") +
  annotate("text", x = dd / 2, y = r2 * 0.85, label = "Boruta", size = 2.8, fontface = "bold") +
  coord_equal() + theme_void()
save_fig(p_venn, file.path(DIR_FIG, "fig3_venn"), 3.6, 2.6)

saveRDS(list(lasso_min = lasso_min, lasso_1se = lasso_1se,
             boruta = bor_confirmed, final = final_vars,
             stability = freq, cv_fit = cv_fit, bor_stats = bor_stats),
        file.path(DIR_OUT, "02_selection.rds"))
