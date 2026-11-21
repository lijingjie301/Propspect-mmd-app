# PROPSPECT-MMD — competing-risks prediction of preoperative stroke in paediatric moyamoya disease

Analysis code and the Shiny explorer for the study. The outcome is preoperative
stroke (cause 1) with revascularisation handled as a competing event (cause 2).
The final model uses five predictors — time to peak (relative, unitless),
PCA involvement, symptom type, symptom frequency and Suzuki stage — and provides
individualised risk estimates at 3, 6, 9 and 12 months, with 12 months as the
primary horizon.

## Pipeline

| script | what it does | output |
|---|---|---|
| `00_run_all.R` | runs everything in order | — |
| `01_prepare.R` | outcome from dates, imputation, 7:3 stratified split | `output/01_prepared.rds` |
| `02_select.R` | LASSO + Boruta, intersection, bootstrap stability | `output/02_selection.rds`, fig 1-3 |
| `03_models.R` | CSC, Fine-Gray, competing-risks RSF | `output/03_models.rds` |
| `04_eval.R` | time-dependent AUC, ROC, calibration, DCA, risk strata, CIF | `output/tables/*.csv`, fig 4-7 |
| `05_shap.R` | KernelSHAP at 3, 6 and 12 months | `output/05_shap.rds`, fig 8-9 |
| `06_app.R` | Shiny app for a single patient | local web app |

Helpers live in `R/utils.R` (date parsing, Aalen-Johansen CIF, IPCW versions of
AUC / Brier score / net benefit, plotting theme).

## Model files used by the app

| file | content |
|---|---|
| `model.rds` | fitted competing-risks random survival forest (ntree = 1,000, nodesize = 10, splitrule = "logrank", seed = 2024) |
| `precomputed.rds` | training/validation risks, SHAP background sample and base value (loaded at start-up) |
| `shap_cache.rds` | cached global KernelSHAP values for the validation subset |

The app reads `precomputed.rds` at start-up and defers loading the forest in
`model.rds` until the first prediction, which keeps start-up short. If these
files are absent the app refits the forest from `data/1250.csv`.

The 12-month Youden-optimal cut-off in the training set is **0.237886**
(23.79%); the 6-month cut-off is **0.151276**.

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
`cari`, `thy`, `ane`, `TTP`, plus four date columns (baseline, stroke,
operation, last follow-up) and, when available, the admission date. The csv is
GB18030 encoded. RNF213 genotype and HR-VWI score are not used by the model.

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

The app returns the predicted 12-month cumulative incidence of preoperative
stroke, the risk stratum based on the fixed training cut-off (23.79%), and an
exact KernelSHAP explanation for the five predictors. The calculator is freely
accessible without registration or a password.

Deploy with:

```r
rsconnect::deployApp(appDir = ".", appName = "PRESTO-pMMD",
                     account = "lijingjie301", server = "shinyapps.io")
```

## Reproducing the released model

`model.rds` is the frozen fitted random survival forest used for the paper and by the deployed calculator (ntree = 1,000, nodesize = 10, splitrule = "logrank", seed = 2024; MD5 `6d22e0b0512c4393c9b35fcea1b025f6`). `precomputed.rds` and `shap_cache.rds` hold the precomputed training/validation risks and the global SHAP values used by the app, so individual predictions can be reproduced without refitting.

`benchmark.csv` lists three standard inputs spanning the risk range and the 12-month risk returned by the frozen model (0.49%, 18.18% and 27.07%). The deployed calculator (<https://lijingjie301.shinyapps.io/PRESTO-pMMD/>) uses the same `model.rds` and returns the same values.

The synthetic example cohort (`data/example_1250.csv`) verifies that the pipeline runs end to end; it does not reproduce the clinical-cohort forest, which is provided as the frozen object above.

## Seeds

42 for the split, 2024 for the forest and the selection, 11/12 for the SHAP
background and explained samples. All set explicitly in the scripts.
