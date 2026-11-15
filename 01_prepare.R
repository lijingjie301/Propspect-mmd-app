# 01 -- Load, split 7:3, impute ----------------------------------------------
# Reads data/1250.csv (development cohort) and data/229.xlsx (external cohort),
# builds the competing-risk outcome, splits the development cohort, and imputes
# missing values. Writes output/01_prepared.rds.

source("R/utils.R")
suppressPackageStartupMessages({
  library(survival); library(readxl); library(mice)
})

DATA_1250 <- if (file.exists("data/1250.csv")) "data/1250.csv" else "data/example_1250.csv"
DATA_229  <- "data/229.xlsx"
DIR_OUT   <- "output"
dir.create(DIR_OUT, showWarnings = FALSE)

SEED <- 42
HORIZONS <- c(3, 6, 12, 24)
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

## ---- imputation ------------------------------------------------------------
# The current tables are complete. The block below only fires if something is
# missing, so the script stays correct for earlier exports that had blanks.
impute_block <- function(dat, vars) {
  na_counts <- colSums(is.na(dat[, vars, drop = FALSE]))
  if (sum(na_counts) == 0) {
    cat("no missing values; imputation skipped\n")
    return(dat)
  }
  cat("missing values detected:\n"); print(na_counts[na_counts > 0])

  # numeric variables -> predictive mean matching; factors stay categorical
  work <- dat[, c(vars, "time", "status")]
  work$status <- factor(work$status)
  set.seed(SEED)
  imp <- mice(work, m = 5, maxit = 5, printFlag = FALSE)
  filled <- complete(imp, 1)
  dat[, vars] <- filled[, vars]
  dat
}

d <- impute_block(d, VARS)

## ---- 7:3 stratified split --------------------------------------------------
# Stratified on the three-level status so the event rate is balanced.
set.seed(SEED)
idx <- sort(unlist(lapply(split(seq_len(nrow(d)), d$status), function(z)
  sample(z, round(0.7 * length(z))))))

train <- d[idx, ]
valid <- d[-idx, ]

cat(sprintf("training %d (%d strokes) / internal validation %d (%d strokes)\n",
            nrow(train), sum(train$status == 1),
            nrow(valid), sum(valid$status == 1)))

saveRDS(list(train = train, valid = valid, external = e,
             vars = VARS, horizons = HORIZONS, seed = SEED),
        file.path(DIR_OUT, "01_prepared.rds"))
