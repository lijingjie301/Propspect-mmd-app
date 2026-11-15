# PROPSPECT-MMD — competing-risks prediction of stroke in paediatric moyamoya disease

Analysis code and the Shiny explorer for the study. The outcome is stroke
(cause 1) with revascularisation handled as a competing event (cause 2);
prediction is done at 3, 6, 12 and 24 months.

## Pipeline

| script | what it does | output |
|---|---|---|
| `00_run_all.R` | runs everything in order | — |
| `01_prepare.R` | outcome from dates, imputation, 7:3 stratified split | `output/01_prepared.rds` |
| `02_select.R` | LASSO + Boruta, intersection, bootstrap stability | `output/02_selection.rds`, fig 1-3 |
| `03_models.R` | CSC, Fine-Gray, competing-risks RSF | `output/03_models.rds` |
| `04_eval.R` | time-dependent AUC, ROC, calibration, DCA, risk strata, CIF | `output/tables/*.csv`, fig 4-7 |
| `05_shap.R` | KernelSHAP at 12 months | `output/05_shap.rds`, fig 8-9 |
| `06_app.R` | Shiny app for a single patient | local web app |

Helpers live in `R/utils.R` (date parsing, Aalen-Johansen CIF, IPCW versions of
AUC / Brier score / net benefit, plotting theme).

## Data

`data/example_1250.csv` is a synthetic cohort of the same size and structure
(covariates bootstrapped from the real file, event times re-simulated from a
fitted Cox model), so the pipeline runs end to end without patient data.
Event rates and predicted risks are close to, but not identical with, the study
cohort.

To reproduce the paper, put the real files in `data/1250.csv` and
`data/229.xlsx`; `01_prepare.R` picks them up automatically and only falls back
to the synthetic file when they are absent.

Columns expected: `event`, `age`, `suzuki`, `pci`, `symptom`, `pinlv`,
`chixushijian`, `sex`, `history`, `hypertension`, `tdm`, `hyperli`, `hhcy`,
`cari`, `thy`, `ane`, `rnf213`, `TTP`, `vwi`, plus four date columns
(baseline, stroke, operation, last follow-up). The csv is GB18030 encoded.

## Requirements

R >= 4.3

```r
install.packages(c("survival", "readxl", "mice", "glmnet", "Boruta",
                   "randomForestSRC", "riskRegression", "prodlim", "cmprsk",
                   "ggplot2", "patchwork", "kernelshap", "shapviz",
                   "shiny", "bslib", "plotly"))
```

## Run

```bash
Rscript 00_run_all.R          # full pipeline
Rscript 06_app.R              # or source("06_app.R") in RStudio
```

Figures are written as PNG (600 dpi) and PDF into `output/figures/`.

## Shiny app

Deployed version: <https://lijingjie301.shinyapps.io/PRESTO-pMMD/>

The app fits the forest at start-up, then explains a single patient with exact
KernelSHAP (7 features, 2^7 coalitions). If `shap_cache.rds` is present it is
loaded instead of recomputing the global explanation.

Deploy with:

```r
rsconnect::deployApp(appDir = ".", appName = "PRESTO-pMMD",
                     account = "lijingjie301", server = "shinyapps.io")
```

## Seeds

42 for the split, 2024 for the forest and the selection, 11/12 for the SHAP
background and explained samples. All set explicitly in the scripts.
