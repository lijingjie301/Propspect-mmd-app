# Run the whole pipeline ------------------------------------------------------
# Usage:  Rscript 00_run_all.R
# Each step is also runnable on its own once the previous objects exist.

t0 <- Sys.time()

source("01_prepare.R")
source("02_select.R")
source("03_models.R")
source("04_eval.R")
source("05_shap.R")

cat(sprintf("\nall steps finished in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
