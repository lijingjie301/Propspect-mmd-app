# PROPSPECT-MMD — competing-risks prediction of preoperative stroke in paediatric moyamoya disease

Analysis code and Shiny explorer for the study. The outcome is preoperative stroke
(cause 1), with revascularisation treated as a competing event (cause 2). The final
model uses five predictors — time to peak (relative, unitless), PCA involvement,
symptom type, symptom frequency and Suzuki stage — and provides individualised risk
estimates at 3, 6, 9 and 12 months.

This repository contains **code and analysis parameters only**. Patient-level data,
fitted model objects and precomputed outputs are not distributed; they stay on the
study server because they contain or derive from clinical records.

## Pipeline

| script | what it does | output |
|---|---|---|
| `00_run_all.R` | runs everything in order | — |
| `01_prepare.R` | outcome from dates, 7:3 stratified split, training-set imputation applied to validation | `output/01_prepared.rds` |
| `02_select.R` | LASSO + Boruta, intersection, bootstrap stability | `output/02_selection.rds`, figures |
| `03_models.R` | cause-specific Cox, Fine–Gray and competing-risks RSF | `output/03_models.rds` |
| `04_eval.R` | time-dependent AUC, ROC, calibration, DCA, risk strata, CIF | `output/tables/*.csv`, figures |
| `05_shap.R` | KernelSHAP at 3, 6 and 12 months | `output/05_shap.rds`, figures |
| `06_app.R` | Shiny app for a single patient | local web app |

Helpers live in `R/utils.R` (date parsing, Aalen–Johansen CIF, IPCW versions of
AUC, Brier score and net benefit, plotting theme).

## Analysis parameters

| setting | value |
|---|---|
| split | 70:30 stratified on three-level outcome status, seed 42 |
| imputation | random-forest imputation fitted in the training set (100 trees), with outcome and follow-up excluded; training missing TTP filled by 5-fold out-of-fold predictions, and the fitted imputer applied unchanged to the internal validation and external cohorts |
| horizons | 3, 6, 9 and 12 months |
| bootstrap | 300 patient-level resamples, seed 2024 |
| RSF | ntree = 1,000, nodesize = 10, splitrule = "logrank", seed = 2024 |
| 12-month cut-off | 0.237886 (23.79%) |
| 6-month cut-off | 0.151276 |

The evaluation target is the observable outcome under the current care pathway:
cases are preoperative strokes by the horizon; controls are children known to be
stroke-free, including those operated on before the horizon without a preceding
stroke; only true loss to follow-up receives inverse-probability-of-censoring
weights.

## Data

Patient data are not included. To run the pipeline, place the development cohort at
`data/1250.csv` and the external cohort at `data/229.xlsx` on the local machine.
The expected columns of `data/1250.csv` are `event`, `age`, `suzuki`, `pci`,
`symptom`, `pinlv`, `chixushijian`, `sex`, `history`, `hypertension`, `tdm`,
`hyperli`, `hhcy`, `cari`, `thy`, `ane`, `TTP`, four date columns (baseline, stroke,
operation, last follow-up) and, when available, the admission date; the file is
GB18030 encoded. RNF213 genotype and HR-VWI score are not used by the model.

Fitted objects (`model.rds`, `precomputed.rds`, `shap_cache.rds`) are generated
locally by the pipeline and are deliberately excluded from this repository.

## Requirements

R >= 4.3

```r
install.packages(c("survival", "readxl", "glmnet", "Boruta",
                   "randomForestSRC", "riskRegression", "prodlim", "cmprsk",
                   "pROC", "ggplot2", "patchwork", "kernelshap", "shapviz",
                   "shiny", "bslib", "plotly"))
```

## Run

```bash
Rscript 00_run_all.R          # full pipeline
Rscript 06_app.R              # or source("06_app.R") in RStudio
```

Figures are written as PNG and PDF into `output/figures/`.

## Shiny app

Deployed version: <https://lijingjie301.shinyapps.io/PRESTO-pMMD/>

The calculator returns the predicted 12-month cumulative incidence of preoperative
stroke, the risk stratum based on the fixed training cut-off (23.79%), and KernelSHAP
values for the five predictors. The deployed application is hosted separately and
keeps its model objects server-side; they are not part of this repository.

To deploy with your local model objects:

```r
rsconnect::deployApp(appDir = ".", appName = "PRESTO-pMMD",
                     account = "lijingjie301", server = "shinyapps.io")
```
