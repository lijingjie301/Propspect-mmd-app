# 01 -- Load, split 7:3 first, then fit the imputer in the training set and apply it unchanged to validation ----
# Reads data/1250.csv (development cohort) and data/229.xlsx (external cohort),
# builds the competing-risk outcome, splits the development cohort, imputes missing
# TTP with a random-forest imputer fitted only in the training set, and writes
# output/01_prepared.rds.
#
# The imputer is fitted in the training set, used with 5-fold out-of-fold predictions
# for training patients with missing TTP, and then applied unchanged to the internal
# validation and external cohorts. Outcome and follow-up time are never used.

source("R/utils.R")
suppressPackageStartupMessages({
  library(survival); library(readxl); library(randomForestSRC)
})

DATA_1250 <- "data/1250.csv"
DATA_229  <- "data/229.xlsx"
if (!file.exists(DATA_1250) || !file.exists(DATA_229)) {
  stop("Place the study data at data/1250.csv and data/229.xlsx; patient data are not distributed with this repository.")
}
DIR_OUT   <- "output"
dir.create(DIR_OUT, showWarnings = FALSE)

SEED <- 42
HORIZONS <- c(3, 6, 9, 12)
VARS <- c("age", "suzuki", "pci", "symptom", "pinlv",
          "chixushijian", "sex", "history", "hypertension", "tdm",
          "hyperli", "hhcy", "cari", "thy", "ane", "vwi", "TTP")

## ---- date columns ----------------------------------------------------------
# The four follow-up dates are not always in the same position, so they are
# located by name (Chinese headers) with a positional fallback.
find_date_cols <- function(dat) {
  nm <- colnames(dat)
  grab <- function(pat) {
    hit <- grep(pat, nm)
    if (length(hit)) nm[hit[1]] else NA_character_
  }
  out <- c(grab("基线"), grab("卒中"), grab("手术"), grab("随访|末次"))
  if (any(is.na(out)) && ncol(dat) >= 4) out <- tail(nm, 4)
  out
}

## ---- development cohort ----------------------------------------------------
# The csv is GB18030 encoded.
d <- read.csv(DATA_1250, fileEncoding = "GB18030", check.names = FALSE,
              stringsAsFactors = FALSE)
date_cols <- find_date_cols(d)
stopifnot(length(date_cols) == 4, !any(is.na(date_cols)))
cat("date columns:", paste(date_cols, collapse = ", "), "\n")

oc <- build_outcome(d, date_cols[1], date_cols[2], date_cols[3], date_cols[4])
d$time <- oc$time
d$status <- oc$status
d <- d[!is.na(d$time) & d$time > 0, ]

cat(sprintf("development cohort: %d patients, %d strokes, %d operations\n",
            nrow(d), sum(d$status == 1), sum(d$status == 2)))

## ---- external cohort -------------------------------------------------------
e <- as.data.frame(read_excel(DATA_229))
edate <- find_date_cols(e)
oc_e <- build_outcome(e, edate[1], edate[2], edate[3], edate[4])
e$time <- oc_e$time
e$status <- oc_e$status
e <- e[!is.na(e$time) & e$time > 0, ]

cat(sprintf("external cohort: %d patients, %d strokes\n",
            nrow(e), sum(e$status == 1)))

## ---- imputation utilities --------------------------------------------------
# A random-forest regression imputer is used for TTP. We use randomForestSRC
# because it returns a fitted forest that can be applied to new patients; the
# missForest package does not provide an apply/predict interface for new data.
rf_formula <- function(target, predictors) {
  as.formula(paste0("`", target, "` ~ ",
                    paste0("`", predictors, "`", collapse = " + ")))
}

fit_imputer <- function(train, vars, seed, targets = NULL, ntree = 100, nodesize = 5) {
  if (is.null(targets)) {
    targets <- vars[vapply(train[vars], anyNA, logical(1))]
  }
  targets <- intersect(vars, targets)
  if (!length(targets)) {
    cat("no missing values in the training set; imputation skipped\n")
    return(list())
  }
  # This script's imputation step is written for TTP; other predictors are complete.
  other_missing <- setdiff(vars[vapply(train[vars], anyNA, logical(1))], targets)
  other_missing <- union(other_missing, setdiff(targets, "TTP"))
  if (length(other_missing)) {
    stop("Unexpected missing values outside TTP: ",
         paste(other_missing, collapse = ", "),
         ". Extend fit_imputer before running.")
  }

  imputers <- list()
  for (v in targets) {
    xv <- setdiff(vars, v)
    obs <- !is.na(train[[v]]) & complete.cases(train[, xv, drop = FALSE])
    if (sum(obs) < 25) {
      stop("Too few observed values to fit the imputer for ", v, ".")
    }
    fit <- rfsrc(
      rf_formula(v, xv),
      data = train[obs, c(v, xv), drop = FALSE],
      ntree = ntree,
      nodesize = nodesize,
      mtry = max(1L, floor(sqrt(length(xv)))),
      seed = seed
    )
    imputers[[v]] <- list(fit = fit, predictors = xv)
    cat(sprintf("fitted training-set imputer for %s using %d observed patients\n",
                v, sum(obs)))
  }
  imputers
}

impute_training_oof <- function(train, imputers, seed, k = 5) {
  if (!length(imputers)) return(train)
  set.seed(seed)
  folds <- sample(rep(seq_len(k), length.out = nrow(train)))
  for (v in names(imputers)) {
    miss <- which(is.na(train[[v]]))
    if (!length(miss)) next
    pred <- rep(NA_real_, length(miss))
    for (f in seq_len(k)) {
      rows_missing <- miss[folds[miss] == f]
      if (!length(rows_missing)) next
      rows_fit <- which(!is.na(train[[v]]) & folds != f)
      if (length(rows_fit) < 25) {
        stop("Too few training rows to fit fold ", f, " for ", v, ".")
      }
      fit_f <- rfsrc(
        rf_formula(v, imputers[[v]]$predictors),
        data = train[rows_fit, c(v, imputers[[v]]$predictors), drop = FALSE],
        ntree = 100, nodesize = 5,
        mtry = max(1L, floor(sqrt(length(imputers[[v]]$predictors)))),
        seed = seed + f
      )
      pred[match(rows_missing, miss)] <-
        as.numeric(predict(fit_f, newdata = train[rows_missing, imputers[[v]]$predictors,
                                                   drop = FALSE])$predicted)
    }
    if (anyNA(pred)) stop("Out-of-fold imputation failed for ", v, ".")
    train[[v]][miss] <- pred
    cat(sprintf("filled %d training-set missing values for %s using out-of-fold predictions\n",
                length(miss), v))
  }
  # Refit the final imputer on all observed training values for applying to validation.
  train
}

apply_imputer <- function(dat, imputers) {
  if (!length(imputers)) return(dat)
  for (v in names(imputers)) {
    miss <- which(is.na(dat[[v]]))
    if (!length(miss)) next
    pred <- as.numeric(predict(imputers[[v]]$fit,
                               newdata = dat[miss, imputers[[v]]$predictors, drop = FALSE])$predicted)
    if (anyNA(pred) || any(!is.finite(pred))) {
      stop("Imputation produced non-finite values for ", v, ".")
    }
    dat[[v]][miss] <- pred
    cat(sprintf("imputed %d values for %s in validation/external data using the training-set imputer\n",
                length(miss), v))
  }
  dat
}

## ---- 7:3 stratified split --------------------------------------------------
# Stratified on the three-level status so the event rate is balanced.
set.seed(SEED)
idx <- sort(unlist(lapply(split(seq_len(nrow(d)), d$status), function(z)
  sample(z, round(0.7 * length(z))))))

train_raw <- d[idx, , drop = FALSE]
valid_raw <- d[-idx, , drop = FALSE]
ext_raw   <- e

## ---- train-only imputation, then apply to validation -----------------------
# Determine which variables need an imputer from all cohorts, but fit only on training.
missing_targets <- unique(c(
  VARS[vapply(train_raw[VARS], anyNA, logical(1))],
  VARS[vapply(valid_raw[VARS], anyNA, logical(1))],
  VARS[vapply(ext_raw[VARS],   anyNA, logical(1))]
))
imputers <- fit_imputer(train_raw, VARS, SEED, targets = missing_targets)
train <- impute_training_oof(train_raw, imputers, SEED)
valid <- apply_imputer(valid_raw, imputers)
ext   <- apply_imputer(ext_raw, imputers)

if (anyNA(train[, VARS]) || anyNA(valid[, VARS]) || anyNA(ext[, VARS])) {
  stop("Missing values remain after imputation.")
}

cat(sprintf("training %d (%d strokes) / internal validation %d (%d strokes)\n",
            nrow(train), sum(train$status == 1),
            nrow(valid), sum(valid$status == 1)))

saveRDS(list(train = train, valid = valid, external = ext,
             vars = VARS, horizons = HORIZONS, seed = SEED,
             imputation = list(method = "training-set random-forest imputation",
                               target = "TTP", ntree = 100,
                               training = "5-fold out-of-fold predictions",
                               validation = "training-set imputer applied unchanged")),
        file.path(DIR_OUT, "01_prepared.rds"))
