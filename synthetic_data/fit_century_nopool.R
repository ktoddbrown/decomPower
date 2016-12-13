source("simulate_data_century.R")
num_rep <- 1;
ts <- c(seq(from=1/360, to=7/360, by=1/360), seq(from=14/360, to=28/360, by=7/360),
        seq(from=30/360, to=360/360, by=30/360));
init_co2 <- c(100, 100, 800);
data<-simulate_data_century(ts, init_co2, num_rep);
data$eCO2 <- data$evolved_co2_mat[,1];

library(rstan);
rstan_options(auto_write = TRUE);
options(mc.cores = parallel::detectCores());
fit <- stan("century_nopool.stan", data=data, iter=1000, seed=1234);
