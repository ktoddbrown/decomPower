simulate_data_century <- function(ts, init_co2, num_rep) {
library(deSolve)
genDerivs <-function(t, Ct, params) {
  # General diff eq model: dC_dt = I(t) + A(t)*C(t)
  # INPUTS:
  #   t: time
  #   Ct: the value of the vector C at time t, C(t)
  #   params: it has two fields, params$I and params$A
  dC_dt = params$I + params$A %*% Ct;
  return(list(dC_dt));
}
calc_co2_flux <- function(co2_t, A) {
  r1 <- sum(A[,1]);
  r2 <- sum(A[,2]);
  r3 <- sum(A[,3]);
  co_flux <- -r1 * co2_t[, 1] - r2 * co2_t[, 2] - r3 * co2_t[, 3];
  return(co_flux);
}
# Generating parameters: 
m <- 3
C0 <- matrix(init_co2, nrow=m);
turnover <- c(1.5, 25, 1000);
K <- 1/turnover;
I <- rep(0, m);
N_t <- length(ts);
co2_flux_mat <- matrix(NA, nrow = N_t, ncol = num_rep);
evolved_co2_mat <- matrix(NA, nrow = N_t, ncol = num_rep);
for (rep in 1:num_rep) {
  Alpha <- matrix(0, m, m);
  Alpha[2, 1] = 0.5;
  Alpha[3, 1] = 0.004;
  Alpha[1, 2] = 0.42;
  Alpha[3, 2] = 0.03;
  Alpha[1, 3] = 0.45;
  Alpha <- Alpha * matrix(rnorm(9, 1, 0.05), nrow=3); 
  A <- Alpha * matrix(rep(K, m), nrow = m, byrow = TRUE) - diag(K);
  params <- list(I=I, A=A);
  t0 <- 0;
  # Solving the ODE system for given parameters:
  ode_data<-ode(y = C0, func = genDerivs,
          times = c(t0,ts), parms = params);
  # Computing evolved CO2 and plotting it:
  total_co2_t0 = sum(ode_data[1,2:(m+1)]);
  co2_t <- ode_data[2:nrow(ode_data),2:(m+1)];
  evolved_co2 <- total_co2_t0 - rowSums(co2_t);
  evolved_co2_mat[,rep] <- exp(log(evolved_co2) + rnorm(length(evolved_co2),0,.1));
  co2_flux <- calc_co2_flux(co2_t, A);
  co2_flux_mat[,rep] <- co2_flux + rnorm(length(co2_flux),0,1);
}
simulated_data <- list(N_t = N_t, ts = ts, num_rep = num_rep, totalC_t0 = total_co2_t0,
                       t0=t0,evolved_co2_mat=evolved_co2_mat, co2_flux_mat=co2_flux_mat);
return(simulated_data);
}

