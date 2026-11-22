# 03 -- Fit the three competing-risks models ---------------------------------
# CSC (cause-specific Cox), Fine-Gray (subdistribution hazard) and a
# competing-risks random survival forest. Writes output/03_models.rds.

source("R/utils.R")
suppressPackageStartupMessages({
  library(survival); library(riskRegression); library(prodlim)
  library(randomForestSRC)
})

DIR_OUT <- "output"
prep <- readRDS(file.path(DIR_OUT, "01_prepared.rds"))
sel  <- readRDS(file.path(DIR_OUT, "02_selection.rds"))

train <- prep$train; valid <- prep$valid; ext <- prep$external
vars <- sel$final
seed <- prep$seed
horizons <- prep$horizons

cat("final variables:", paste(vars, collapse = ", "), "\n")

fml <- as.formula(paste("Hist(time, status) ~", paste(vars, collapse = " + ")))

m_csc <- CSC(fml, data = train, cause = 1)
m_fgr <- FGR(fml, data = train, cause = 1)

set.seed(seed)
m_rsf <- rfsrc(Surv(time, status) ~ ., data = train[, c("time", "status", vars)],
               ntree = 1000, nodesize = 10, cause = 1, splitrule = "logrank",
               importance = "none")

models <- list(CSC = m_csc, FGR = m_fgr, RSF = m_rsf)

## predicted absolute risk at each horizon
risk_matrix <- function(model, tag, dat, times) {
  if (tag == "RSF") {
    pr <- predict(model, newdata = dat[, vars])
    j <- vapply(times, function(t) which.min(abs(pr$time.interest - t)), integer(1))
    out <- pr$cif[, j, 1]
  } else {
    out <- predictRisk(model, newdata = dat, times = times)
  }
  out <- matrix(out, nrow = nrow(dat))
  colnames(out) <- paste0("t", times)
  out
}

get_risk <- function(dat) {
  out <- lapply(names(models), function(m) risk_matrix(models[[m]], m, dat, horizons))
  names(out) <- names(models)
  out
}
risk <- list(train = get_risk(train), valid = get_risk(valid), external = get_risk(ext))

saveRDS(list(models = models, risk = risk, vars = vars,
             train = train, valid = valid, external = ext,
             horizons = horizons, seed = seed),
        file.path(DIR_OUT, "03_models.rds"))

cat("models fitted and risks stored\n")
