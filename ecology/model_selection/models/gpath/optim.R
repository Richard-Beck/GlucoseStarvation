library(cmdstanr)
library(posterior)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
MODEL_NAME      <- if (length(args) >= 1) args[1] else "gpath"
R_val           <- if (length(args) >= 2) as.integer(args[2]) else 1L
P_val           <- if (length(args) >= 3) as.integer(args[3]) else 1L
W_val           <- if (length(args) >= 4) as.integer(args[4]) else 0L
CONSTRAINT_FLAG <- if (length(args) >= 5) as.integer(args[5]) else 0L 
WASTE_MECH_FLAG <- if (length(args) >= 6) as.integer(args[6]) else 1L 
NUM_THREADS     <- if (length(args) >= 7) as.integer(args[7]) else 4L
TASK_ID         <- if (length(args) >= 8) as.integer(args[8]) else 1L

HOLDOUT <- FALSE
tmp <- strsplit(MODEL_NAME,split="_")[[1]]

if(length(tmp)>1){
  if(tmp[[2]]=="holdout"){
    HOLDOUT <- TRUE
    MODEL_NAME <- tmp[[1]]
  }
}

run_id <- sprintf("%dR_%dP_%dW_C%d_M%d", R_val, P_val, W_val, CONSTRAINT_FLAG, WASTE_MECH_FLAG)
source("ecology/model_selection/models/gpath/gpath.R")

cat(sprintf(">>> RUNNING: %s | Config: %s | task=%d | threads=%d\n", MODEL_NAME, run_id, TASK_ID, NUM_THREADS))

stan_data <- readRDS("data/stan_ready_data.Rds")
grp_keys <- paste(stan_data$line_id, stan_data$ploidy_metric, sep = "_")
unique_grps <- unique(grp_keys)
stan_data$N_groups <- length(unique_grps)
stan_data$group_id <- match(grp_keys, unique_grps)
stan_data$is_train <- rep(1L,stan_data$N_wells)
if(HOLDOUT) stan_data$is_train[stan_data$ploidy_metric>0] <- 0
config <- generate_stan_config(R = R_val, P = P_val, W = W_val, strict_spec = (CONSTRAINT_FLAG == 1L), M = WASTE_MECH_FLAG, base_priors = base_priors)
for (nm in names(config)) stan_data[[nm]] <- config[[nm]]

stan_data$waste_mech <- if (W_val > 0) rep(as.numeric(WASTE_MECH_FLAG), W_val) else numeric(0)

#stan_data$ploidy_metric <- log2(stan_data$ploidy_abs) - 1

R_init_base <- matrix(0.0, nrow = stan_data$N_wells, ncol = R_val)
if (R_val > 1) for (c in 2:R_val) R_init_base[, c] <- 1.0 
stan_data$R_init_base <- R_init_base
stan_data <- stan_data[!sapply(stan_data, is.character)]

stan_file <- file.path("ecology/model_selection/models", MODEL_NAME, paste0(MODEL_NAME, "_hier.stan"))
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE),force_recompile=T)

cat(sprintf("[Task %d]: Optimizing...\n", TASK_ID))
outDir <- file.path("ecology/model_selection/data", MODEL_NAME, run_id, ifelse(HOLDOUT,"holdout","hier"))
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

res <- tryCatch({
  opt <- mod$optimize(data = stan_data, algorithm = "bfgs", init = 2, refresh = 0, threads = NUM_THREADS,
                      save_latent_dynamics = T,jacobian = TRUE,
                      show_messages = TRUE)
  rc <- opt$return_codes()[1]
  ##if opt worked, use it:
  if (rc == 0L) list(rc = rc, lp = opt$lp(), draws = opt$draws(format = "matrix")) else{
    ## if it failed default to text
    txt <- capture.output(opt$output())
    csv_line <- grep("^\\s*file = ", txt, value = TRUE)[2]
    csv <- sub("^\\s*file =\\s*", "", csv_line)
    
    opt <- read_cmdstan_csv(csv)
    list(rc = rc, lp = c(opt$point_estimates[,"lp__"]), draws = opt$point_estimates)
  }})
  
#saveRDS(stan_data, file.path(outDir, "stan_data.Rds")) 
saveRDS(res$draws, file.path(outDir, sprintf("optim_draws_%d.Rds", TASK_ID)))
saveRDS(res$lp, file.path(outDir, sprintf("optim_lp_%d.Rds", TASK_ID)))
#saveRDS(res$rc, file.path(outDir, sprintf("optim_rc_%d.Rds", TASK_ID)))

cat(">>> Done.\n")

