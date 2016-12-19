simulate_data_century <- function(t_meas, t_cap, init_C, num_rep) {
library(deSolve)
library(gtools)
library(VGAM)  
  
genDerivs <-function(t, Ct, params) {
  # General diff eq model: dC_dt = I(t) + A(t)*C(t)
  # INPUTS:
  #   t: time
  #   Ct: the value of the vector C at time t, C(t)
  #   params: it has two fields, params$I and params$A
  dC_dt = params$I + params$A %*% Ct;
  return(list(dC_dt));
}

# Generating parameters: 
m <- 3
C_t0 <- matrix(init_C, nrow=m);
turnover <- c(1.5, 25, 1000);
K <- 1/turnover;
I <- rep(0, m);
N_t <- length(t_meas);
CO2_flux_mat <- matrix(NA, nrow = N_t, ncol = num_rep);
for (rep in 1:num_rep) {
  Alpha <- matrix(0, m, m);
  Alpha_rep <- matrix(0, m, m);
  Alpha[2, 1] = 0.5;
  Alpha[3, 1] = 0.004;
  Alpha[1, 2] = 0.42;
  Alpha[3, 2] = 0.03;
  Alpha[1, 3] = 0.45;
  Alpha[1, 1] = 1 - Alpha[2, 1] - Alpha[3, 1];
  Alpha[2, 2] = 1 - Alpha[1, 2] - Alpha[3, 2];
  Alpha[3, 3] = 1 - Alpha[1, 3] - Alpha[2, 3];
  kappa <- c(100,100,100);
  Alpha_rep[,1] <-  rdirichlet(1, Alpha[,1] * kappa[1]);
  Alpha_rep[,2] <-  rdirichlet(1, Alpha[,2] * kappa[2]);
  Alpha_rep[,3] <-  rdirichlet(1, Alpha[,3] * kappa[3]);
  Alpha[1, 1] <- 0;
  Alpha[2, 2] <- 0;
  Alpha[3, 3] <- 0;
  A <- Alpha * matrix(rep(K, m), nrow = m, byrow = TRUE) - diag(K);
  params <- list(I=I, A=A);
  t0 <- 0;
  # Solving the ODE system for given parameters:
  meas_data<-ode(y = C_t0, func = genDerivs,
          times = c(t0,t_meas), parms = params);
  cap_data<-ode(y = C_t0, func = genDerivs,
                 times = c(t0,t_cap), parms = params);
  # Computing evolved CO2 and plotting it:
  totalC_t0 = sum(meas_data[1,2:(m+1)]);
  CO2_t_meas <- totalC_t0 - rowSums(meas_data[2:nrow(meas_data), 2:(m+1)]);
  CO2_t_cap <- totalC_t0 - rowSums(cap_data[2:nrow(cap_data),2:(m+1)]);
  CO2_flux <- (CO2_t_meas - CO2_t_cap)/(t_meas-t_cap); 
  CO2_flux_mat[,rep] <- exp(log(CO2_flux) + rnorm(length(CO2_flux),0,.1));
}
simulated_data <- list(N_t = N_t, t_meas = t_meas, t_cap = t_cap, 
                       num_rep = num_rep, totalC_t0 = totalC_t0,
                       t0=t0, CO2_flux=CO2_flux_mat);
return(simulated_data);
}

