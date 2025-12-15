## assumes this script is run from GlucoseStarvation project root.
gen_gluc_points <- function(data, probs=c(0.05,0.95), Gmin=exp(-14), Gmax=exp(3.5), n=2000){
  y <- data$glucose
  y <- split(y, f=interaction(y$G0, y$hours, drop=TRUE))
  y <- y[sapply(y, nrow) > 0]
  out <- lapply(y, function(yi){
    Ggrid <- exp(seq(log(Gmin), log(Gmax), length.out=n))
    ll <- sapply(Ggrid, function(G) sum(log(data$ll_lum(G, yi))))          # sum log-likelihoods
    w  <- exp(ll - max(ll))                                                # stabilize
    dG <- c(diff(Ggrid), tail(diff(Ggrid),1))                               # bin widths
    w  <- w * dG; w <- w / sum(w)                                          # proper normalization
    cdf <- cumsum(w)
    lower <- approx(x=cdf, y=Ggrid, xout=probs[1], rule=2, ties="ordered")$y
    upper <- approx(x=cdf, y=Ggrid, xout=probs[2], rule=2, ties="ordered")$y
    est   <- Ggrid[which.max(ll)]
    data.frame(G0=yi$G0[1],hours=yi$hours[1],lower=lower, upper=upper, est=est)
  })
  res <- do.call(rbind, out); rownames(res) <- names(out); res
}

source("scripts/parameter_estimation.R")
source("scripts/rxode_models.R")

input_dir <- "data/model_inputs/"
output_dir <- "data/model_outputs/"

#args <- commandArgs(trailingOnly = TRUE)

#input_file <- as.character(args[1])
#n_cores   <- as.numeric(args[2])
#model_name <- as.character(args[3])

n_cores   <- 32
plan(sequential)
model_name <- "model_B"


par_info <- data.table(
  par   = c("theta", "kp", "kd", "kd2", "g50a", "na", "g50d", "nd", "v1", "v2"),
  lower = c(1e3, 1e-3, 1e-5, 1e-5, 1e-5, 1, 1e-5, 1, 1e-8, 1e-8),
  upper = c(1e5, 1, 10, 10, 5, 10, 1, 10, 1e-3, 1e-3)
)

model <- load_model("model_B")


for(i in 1:4){
  
  input_file <- file.path(input_dir,list.files(input_dir)[i])
  exp_id <- gsub(".Rds","",basename(input_file))
  output_path <- file.path(output_dir,exp_id)
  dir.create(output_path,recursive = T)
  
  data <- readRDS(input_file)
  
  opt <- fit_model_deoptim(data, model,par_info,
                           n_pop = 10 * length(par_info$par),
                           total_gens = 200,
                           cores = n_cores)
  saveRDS(opt,file.path(output_path,"opt.Rds"))
  x <- plot_objective_function(opt$optim$bestmem,model,data,par_info$par,sample_freq = 0.1)
  library(ggplot2)
  p <- ggplot(x,aes(x=hours))+
    facet_grid(rows=vars(G0),scales="free")+
    geom_line(aes(y=NL,color="alive"))+
    geom_line(aes(y=ND,color="dead"))+
    geom_point(data=data$cells,aes(y=N,color="alive"))+
    geom_point(data=data$cells,aes(y=D,color="dead"))+
    ggtitle(exp_id)
  p
  ggsave(file.path(output_path,"cells.png"),p)
  
  gdat <- gen_gluc_points(data)
  
  p <- ggplot(x,aes(x=hours))+
    facet_grid(rows=vars(G0),scales="free")+
    geom_line(aes(y=G))+
    geom_errorbar(data=gdat,aes(ymin=lower,ymax=upper))+
    geom_point(data=gdat,aes(y=est))+
    ggtitle(exp_id)
  p
  ggsave(file.path(output_path,"glucose.png"),p)
  
}



