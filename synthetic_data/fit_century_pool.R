source("simulate_data_century.R")
num_rep <- 3;
t_meas <- c(seq(from=1/360, to=7/360, by=1/360), seq(from=14/360, to=28/360, by=7/360),
            seq(from=60/360, to=360/360, by=30/360));
t_cap <- c(seq(from=0.5/360, to=6.5/360, by=1/360), seq(from=10.5/360, to=24.5/360, by=7/360),
           seq(from=45/360, to=345/360, by=30/360));
init_C <- 1E3*c(.1, .1, .8);
data <- simulate_data_century(t_meas, t_cap, init_C, num_rep);
data$CO2_flux <- data$CO2_flux_mat;

library(rstan);
rstan_options(auto_write = TRUE);
options(mc.cores = parallel::detectCores());
fit <- stan("century_pool.stan", data=data, iter=500, seed=1234);
