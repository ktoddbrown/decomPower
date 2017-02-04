source("simulate_data_century.R")
num_rep <- 5;
t_meas <- c(seq(from=1/360, to=7/360, by=1/360), seq(from=14/360, to=28/360, by=7/360),
            seq(from=60/360, to=360/360, by=30/360));
t_cap <- c(seq(from=0.5/360, to=6.5/360, by=1/360), seq(from=10.5/360, to=24.5/360, by=7/360),
           seq(from=45/360, to=345/360, by=30/360));
init_C <- 1E3*c(.1, .1, .8);
data <- simulate_data_century(t_meas, t_cap, init_C, num_rep);

init <- function()
  list(gamma = c(0.1, 0.1, 0.8),
       turnover = c(1.5, 25, 1000),
       a21 = rep(0.5, 5),
       a31 = rep(0.004, 5),
       a12 = rep(0.42, 5),
       a32 = rep(0.03, 5),
       a13 = rep(0.45, 5))

library(rstan);
rstan_options(auto_write = TRUE);
options(mc.cores = parallel::detectCores());
fit_hier_bdf <- stan("fixed_param/century_hier.stan", data = data, init = init,
                    iter = 1, seed = 1234, chains = 1, algorithm = "Fixed_param")
fit_hier_me <- stan("fixed_param/century_hier_me.stan", data = data, init = init,
                 iter = 1, seed = 1234, chains = 1, algorithm = "Fixed_param")
fit_hier_rk45 <- stan("fixed_param/century_hier_nonstiff.stan", data = data, init = init,
                     iter = 1, seed = 1234, chains = 1, algorithm = "Fixed_param")


## Plot simulated data
fit_me <- extract(fit_hier_me)
fit_bdf <- extract(fit_hier_bdf)
fit_rk45 <- extract(fit_hier_rk45)
library(ggplot2)
CO2_flux_hat <- matrix(nrow = data$N_t, ncol = data$num_rep)
CO2_flux_hat_me <- CO2_flux_hat
CO2_flux_hat_bdf <- CO2_flux_hat
CO2_flux_hat_rk45 <- CO2_flux_hat
for (i in 1:data$num_rep) {
  for (t in 1:data$N_t) {
    CO2_flux_hat_bdf[t, i] <- fit_bdf$CO2_flux_hat_pred[ ,t, i]
    CO2_flux_hat_me[t, i] <- fit_me$CO2_flux_hat_pred[ , t, i]
    CO2_flux_hat_rk45[t, i] <- fit_rk45$CO2_flux_hat_pred[ , t, i]
  }
}

diff_me_bdf <- CO2_flux_hat_bdf - CO2_flux_hat_me
diff_me_rk45 <- CO2_flux_hat_rk45 - CO2_flux_hat_me
diff_bdf_rk45 <- CO2_flux_hat_rk45 - CO2_flux_hat_bdf

max(abs(diff_me_bdf))
max(abs(diff_me_rk45))
max(abs(diff_bdf_rk45))
