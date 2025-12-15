library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)
library(ggplot2)

# ---- compile / data ----
stan_data <- readRDS("data/stan_ready_data.Rds")
stan_data$mode <- 1

model <- cmdstan_model("STAN/model_B.stan", quiet = TRUE)

# ---- sample from priors via HMC (mode=1 => priors only; no ODE in GQ) ----
fit_prior <- model$sample(
  data = stan_data,
  chains = 1,
  parallel_chains = 1,
  iter_sampling = 2000,
  iter_warmup = 50,
  fixed_param = FALSE,
  seed = 123,
  refresh = 0
)

draws_df <- as_draws_df(fit_prior$draws())

# ---- helpers matching Stan ----
softcap <- function(x, cap) cap - log1p(exp(cap - x))
cap_log_main <- 40.0
cap_log_hill <- 6.0
par10_names <- c("theta","kp","kd","kd2","g50a","na","g50d","nd","v1","v2")

get_col <- function(nm) {
  if (!nm %in% names(draws_df)) stop("Missing in draws: ", nm)
  draws_df[[nm]]
}

# ---- collapse wells -> 10 groups by (line_id, is_high_ploidy) ----
groups <- tibble(
  well = seq_len(stan_data$N_wells),
  line_id = stan_data$line_id,
  high = stan_data$is_high_ploidy
) %>%
  distinct(line_id, high) %>%
  arrange(line_id, high) %>%
  mutate(group = row_number(),
         group_label = paste0("line ", line_id, " | high=", high))

N_draws <- nrow(draws_df)

# ---- build ODE-parameter draws per group ----
out <- vector("list", nrow(groups))

for (gi in seq_len(nrow(groups))) {
  l <- groups$line_id[gi]
  h <- groups$high[gi]
  
  p_mat <- matrix(NA_real_, nrow = N_draws, ncol = 10)
  colnames(p_mat) <- par10_names
  
  for (pp in 1:10) {
    raw <- get_col(sprintf("mu_global[%d]", pp)) +
      get_col(sprintf("sigma_line[%d]", pp)) * get_col(sprintf("z_line[%d,%d]", pp, l)) +
      get_col(sprintf("beta_high[%d]", pp)) * h
    
    if (pp %in% c(6, 8)) p_mat[, pp] <- 1.0 + exp(softcap(raw, cap_log_hill))
    else                 p_mat[, pp] <- exp(softcap(raw, cap_log_main))
  }
  
  df <- as.data.frame(p_mat) %>%
    mutate(.draw = seq_len(N_draws),
           group = groups$group_label[gi])
  
  out[[gi]] <- pivot_longer(df, cols = all_of(par10_names), names_to = "param", values_to = "value")
}

p_long <- bind_rows(out) %>%
  mutate(
    group = factor(group, levels = groups$group_label),
    param = factor(param, levels = par10_names)
  )

# ---- plot: rows=group, cols=param ----
ggplot(p_long, aes(x = value)) +
  geom_histogram(bins = 50) +
  facet_grid(group ~ param, scales = "free_x") +
  scale_x_log10() +
  theme_bw() +
  labs(x = NULL, y = "count", title = "ODE parameter values seen by solver (mode=1 priors), by (line_id, high)")
