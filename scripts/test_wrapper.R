library(rxode2)
library(data.table)
# Assumes 'parameter_estimation.R' loads 'future' and 'DEoptim'

argnum    <- 1
model_name <- "model_B"

cat("--- Job Index:", argnum, "---\n")
cat("--- Model:", model_name, "---\n")


# --- 3. Load and Prepare Data ---

library(data.table)

#' @title Prepare Joint Dataset for Model Fitting
#' @description Merges time-course cell count data with sparse glucose measurements.
#' @param cell_count_dt A data.table with cell counts, including columns:
#'   `cellLine`, `ploidy`, `glucose` (initial concentration), and `hours`.
#' @param glucose_dt A data.table with glucose measurements, including columns:
#'   `cellLine`, `ploidy`, `G0` (initial concentration), `Day`, and `G` (measured concentration).
#' @return A single data.table with all columns from `cell_count_dt` plus a
#'   new `glucose_measured` column. Rows without a direct glucose measurement
#'   at a specific time point will have `NA` in this column.

prepare_joint_data <- function(cell_count_dt, glucose_dt) {
  
  # --- 1. Prepare the glucose data for merging ---
  glucose_to_merge <- copy(glucose_dt)
  
  # Rename columns to match the cell count data's keys
  setnames(glucose_to_merge, old = c("G0", "G","CellLine"), new = c("glucose", "glucose_measured","cellLine"))
  
  # Convert 'Day' to 'hours'
  glucose_to_merge[, hours := Day * 24]
  
  # --- 2. Perform a left merge ---
  # This keeps all rows from the cell count data. Where keys match,
  # 'glucose_measured' is filled. Otherwise, it becomes NA.
  merged_dt <- merge(
    cell_count_dt,
    glucose_to_merge[, .(cellLine, ploidy, glucose, hours, glucose_measured)],
    by = c("cellLine", "ploidy", "glucose", "hours"),
    all.x = TRUE
  )
  
  return(merged_dt)
}

setwd("/share/lab_crd/lab_crd/HighPloidy_CostBenefits/data/GlucoseStarvation/")

smm <- readRDS("data/counts/corrected.Rds")
smm$glucose <- as.numeric(smm$glucose)
smm$hours <- 24 * smm$time_numeric
smm$alive <- smm$N - smm$D
smm$dead <- smm$D

setDT(smm)
smm[, c("time_numeric", "N", "D") := NULL]


y <- data.table::fread("data/glucose/processed/250925_Glucose_Promega_chemicalSUM159-SNU668_Andor.csv")

tmp <- y[is.na(y$ploidy),]
y$ploidy[is.na(y$ploidy)]<-"2N"
tmp$ploidy <- "4N"
y <- rbind(tmp,y)
y$CellLine[y$CellLine=="SUM159"] <- "SUM-159-chem"
y$ploidy[y$ploidy=="4N"&y$CellLine=="SNU668"] <- "high"
y$ploidy[y$ploidy=="2N"&y$CellLine=="SNU668"] <- "low"

smm <- smm[smm$cellLine%in%y$CellLine,]

smm <- prepare_joint_data(smm,y)


# --- 4. Split Data and Select Subset for This Job ---
groups_to_fit <- split(smm, by = c("cellLine", "ploidy"))

if (argnum > length(groups_to_fit) || argnum < 1) {
  stop(paste("SLURM JOB ERROR: Job index", argnum, "is out of bounds. There are only", length(groups_to_fit), "groups."), call. = FALSE)
}

smm_subset <- groups_to_fit[[argnum]]
group_name <- names(groups_to_fit)[argnum]
cat("--- Starting SLURM job for group:", group_name, "---\n")


# --- 5. Source Fitting Functions & Define Parameter Bounds ---
source("scripts/parameter_estimation.R") # Contains the DEoptim fitting functions and models

# Define parameter bounds for intra_model
par_bounds_intra <- data.table(
  par   = c("Vmax_uptake", "Km_uptake", "kp", "Y_xr", "m_r", "kdStarv", "kw", "r_half_g", "nr_g", "r_half_d", "nr_d","r_cell"),
  lower = c(1e-6, 1e-2, 1e-3, 1e2, 1e-8, 1e-4, 1e-10, 1e-6, 1, 1e-8, 1,1e-10),
  upper = c(1e-2, 1, 0.2, 1e6, 1e-4, 0.5, 1e-4, 1e-2, 8, 1e-3, 8,1)
)

# Define global parameter bounds for model_B
par_bounds_b <- data.table(
  par   = c("theta", "kp", "kd", "kd2", "g50a", "na", "g50d", "nd", "v1", "v2"),
  lower = c(1e3, 1e-3, 1e-5, 1e-5, 1e-5, 1, 1e-5, 1, 1e-8, 1e-8),
  upper = c(1e5, 1, 10, 10, 5, 10, 1, 10, 1e-3, 1e-3)
)


# --- 6. Set Model-Specific Variables based on Argument ---
if (model_name == "intra_model") {
  model_to_fit <- intra_model
  par_bounds_to_use <- par_bounds_intra
  output_prefix <- "data/models/model_intra_deoptim_"
} else if (model_name == "model_B") {
  if (!exists("model_B")) stop("Model 'model_B' is not defined in the script.")
  model_to_fit <- model_B
  par_bounds_to_use <- par_bounds_b
  output_prefix <- "data/models/model_B_deoptim_"
} else {
  stop("Invalid model name provided. Choose 'intra_model' or 'model_B'.")
}

log_lower_bounds <- log(pmax(par_bounds_to_use$lower, 1e-12))
log_upper_bounds <- log(pmax(par_bounds_to_use$upper, 1e-12))
# ---------------------------------------------------------------

data_subsets <- split(smm_subset, by = "glucose")

tst <- lapply(1:10,function(i){
  pars <- runif(length(log_lower_bounds),min = log_lower_bounds,max = log_upper_bounds)
  res <- deoptim_objective_function(pars,model = model_to_fit,data_subsets,par_names=par_bounds_to_use$par)
  list(pars=pars,res=res)
})


