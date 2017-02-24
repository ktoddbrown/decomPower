.libPaths('/vega/stats/users/mk3971/rpackages/')
library(rstan)
num_rep <- 5
t_meas_all <- c(seq(from=1/360, to=14/360, by=0.5/360), seq(from=15/360, to=30/360, by=1/360),
            seq(from=32/360, to=60/360, by=2/360), seq(from=67/360, to=360/360, by=7/360))
t_cap_all <-  head(t_meas_all, -1) + 0.85 * (tail(t_meas_all, -1) - head(t_meas_all, -1))
t_cap_all <- c(t_meas_all[1] - (t_meas_all[2] - t_cap_all[1]), t_cap_all)
sm <- stan_model("century_hier_nonstiff.stan")
downsample_set <- c(1,2,4,5,10,20,50)
for (rep in 1:3) {
  for (i in 1:length(downsample_set)) {
    i
    n_downsample <- downsample_set[i]
    t_meas <- t_meas_all[seq(1,length(t_meas_all), n_downsample)]
    t_cap <- t_cap_all[seq(1,length(t_cap_all), n_downsample)]
    gamma_global <- 100*c(.1, .1, .8); # global initializations of carbons in the three pools
    source("simulate_data_century.R")
    data <- simulate_data_century(t_meas, t_cap, gamma_global, num_rep);
    rstan_options(auto_write = TRUE);
    options(mc.cores = parallel::detectCores());
    fit <- sampling(sm, data=data, iter=3000, seed=1234, control = list(adapt_delta=0.95));
    saveRDS(c(data=data, fit=fit), file=paste("FITS/fit_",i,"_",rep,".rds",sep=""))
  }
}

