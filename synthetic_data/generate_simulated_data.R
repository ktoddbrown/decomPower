.libPaths('/vega/stats/users/mk3971/rpackages/')
args = commandArgs(trailingOnly=TRUE); # command-line inputs
n_downsample = as.integer(args[1]);    # factor by which the data is down-sampled
num_rep <- as.integer(args[2]);        # number of replicates
t_meas <- c(seq(from=1/360, to=14/360, by=0.5/360), seq(from=15/360, to=30/360, by=1/360),
            seq(from=32/360, to=60/360, by=2/360), seq(from=67/360, to=360/360, by=7/360));
t_cap <-  head(t_meas, -1) + 0.85 * (tail(t_meas, -1) - head(t_meas, -1));
t_cap <- c( t_meas[1] - (t_meas[2] - t_cap[1]), t_cap);
t_meas <- t_meas[seq(1,length(t_meas), n_downsample)]
t_cap <- t_cap[seq(1,length(t_cap), n_downsample)]
gamma_global <- 100*c(.1, .1, .8); # global initializations of carbons in the three pools
source("simulate_data_century.R")
data <- simulate_data_century(t_meas, t_cap, gamma_global, num_rep);
