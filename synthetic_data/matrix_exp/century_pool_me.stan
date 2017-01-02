// Fitting the century model with complete pooling
// This version uses a matrix exponential solution instead of
// an ODE integrator.
functions { 
  /** 
  * Compute evolved CO2 from the system given the specified 
  * parameters and times. This is done by solving the century
  * model ODE system with a matrix exponential solution and
  * then calculating the rate CO2 is emmited.
  * 
  * @param N_t number of times
  * @param t0 initial time
  * @param ts times
  * @param gamma partitioning coefficient
  * @param k decomposition rates 
  * @param ajk transfer rates
  * @return evolved CO2 for times ts 
  */ 
    vector evolved_CO2(int N_t, real t0, vector ts, 
                       vector gamma, real totalC_t0, 
                       vector k, real a21, real a31, real a12, 
                       real a32, real a13) { 
      vector[3] C_t0;        // initial state 
      matrix[3, 3] A;        // ODE matrix
      vector[3] C_t[N_t];    // predicted pool content
      vector[N_t] CO2_t;     // evolved CO2 at times ts

      A[1, 1] = -k[1];
      A[1, 2] = a12 * k[2];
      A[1, 3] = a13 * k[3];
      A[2, 1] = a21 * k[1];
      A[2, 2] = -k[2];
      A[2, 3] = 0;
      A[3, 1] = a31 * k[1];
      A[3, 2] = a32 * k[2];
      A[3, 3] = -k[3];

      C_t0 = gamma * totalC_t0;

      for (t in 1:N_t) {
        C_t[t] = matrix_exp(ts[t] * A) * C_t0;
        CO2_t[t] = totalC_t0 - sum(C_t[t]); 
      }

      return CO2_t; 
    } 
} 
data { 
  real<lower=0> totalC_t0;      // initial total carbon 
  real t0;                      // initial time 
  int<lower=0> N_t;             // number of measurement times 
  int<lower=0> num_rep;         // number of replicates 
  vector<lower=t0>[N_t] t_meas; // measurement times 
  vector<lower=t0>[N_t] t_cap;  // cap times 
  matrix<lower=0>[N_t, num_rep] CO2_flux; // measured carbon fluxes 
} 
transformed data { 
  real x_r[0];  // no real data for ODE system 
  int x_i[0];   // no integer data for ODE system 
}  
parameters { 
  vector<lower=0>[3] turnover;  // turnover rates
  simplex[3] gamma;             // partitioning coefficients (a simplex) 
  vector<lower=0>[3] sigma;     // turnover standard deviation 
  real<lower=0> sigma_obs;      // observation standard deviation
  simplex[3] A1;                // output rates from pool 1
  simplex[3] A2;                // output rates from pool 2
  simplex[3] A3;                // output rates from pool 3
} 
transformed parameters { 
  vector<lower=0>[3] k;         // decomposition rates (1/turnover)
  vector[N_t] CO2_meas;         // evolved CO2 at measurement times
  vector[N_t] CO2_cap;          // evolved CO2 at cap times
  vector[N_t] CO2_flux_hat;     // CO2 flux (average evolved CO2 between t_cap & t_meas)
  real<lower=0, upper=1> a21;   // transfer rates
  real<lower=0, upper=1> a31;
  real<lower=0, upper=1> a12;
  real<lower=0, upper=1> a32;
  real<lower=0, upper=1> a13;
  k = 1 ./ turnover;
  // transfer rates are the same for all replications:
  a21 = A1[2];
  a31 = A1[3];
  a12 = A2[1];
  a32 = A2[3];
  a13 = A3[1];
  CO2_meas = evolved_CO2(N_t, t0, t_meas, gamma, totalC_t0, 
                          k, a21, a31, a12, a32, a13); 
  CO2_cap = evolved_CO2(N_t, t0, t_cap, gamma, totalC_t0, 
                          k, a21, a31, a12, a32, a13);
  CO2_flux_hat = (CO2_meas - CO2_cap)./(t_meas - t_cap);
} 
model { 
  // priors 
  turnover[1] ~ normal(1.5, 0.15 * sigma[1]);
  turnover[2] ~ normal(25, 2.5 * sigma[2]);
  turnover[3] ~ normal(1000, 100 * sigma[3]);
  sigma ~ cauchy(0,1); 
  sigma_obs ~ cauchy(0,1);
 
  // likelihood     
  to_vector(CO2_flux) ~ lognormal(to_vector(rep_matrix(log(CO2_flux_hat),num_rep)), sigma_obs);
} 
