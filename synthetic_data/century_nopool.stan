// Fitting the century model with no pooling
functions { 
  /** 
  * ODE system for the Century model with no input fluxes. 
  * @param t time at which derivatives are evaluated. 
  * @param C system state at which derivatives are evaluated. 
  * @param theta parameters for system. 
  * @param x_r real constants for system (empty). 
  * @param x_i integer constants for system (empty). 
  */ 
    real[] century_model(real t, real[] C, real[] theta, 
                         real[] x_r, int[] x_i) { 
      real k[3]; 
      real a21;
      real a31;
      real a12;
      real a32;
      real a13;
      real dC_dt[3]; 
      k = theta[1:3];
      a21 = theta[4];
      a31 = theta[5];
      a12 = theta[6];
      a32 = theta[7];
      a13 = theta[8];
      dC_dt[1] = -k[1] * C[1] + a12 * k[2] * C[2] + a13 * k[3] * C[3];
      dC_dt[2] = -k[2] * C[2] + a21 * k[1] * C[1];
      dC_dt[3] = -k[3] * C[3] + a31 * k[1] * C[1] + a32 * k[2] * C[2]; 
      return dC_dt; 
    } 
  
  /** 
  * Compute evolved CO2 from the system given the specified 
  * parameters and times. This is done by simulating the system 
  * defined by the ODE function century_model and then 
  * calculating the rate CO2 is emmited.
  * 
  * @param N_t number of times
  * @param t0 initial time
  * @param ts times
  * @param gamma partitioning coefficient
  * @param k decomposition rates 
  * @param ajk transfer rates
  * @param x_r real data (empty) 
  * @param x_i integer data (empty) 
  * @return evolved CO2 for times ts 
  */ 
    vector evolved_CO2(int N_t, real t0, vector ts, 
                       vector gamma, real totalC_t0, 
                       vector k, real a21, real a31, real a12, 
                       real a32, real a13, real[] x_r, int[] x_i) { 
      
      real C_t0[3];         // initial state 
      real theta[8];        // ODE parameters 
      real C_t[N_t,3];      // predicted pool content 
      vector[N_t] CO2_t;    // evolved CO2 at times ts
      
      C_t0 = to_array_1d(gamma*totalC_t0);
      theta[1:3] = to_array_1d(k);
      theta[4] = a21;
      theta[5] = a31;
      theta[6] = a12;
      theta[7] = a32;
      theta[8] = a13;
      C_t = integrate_ode_bdf(century_model,  
                              C_t0, t0, to_array_1d(ts), theta, x_r, x_i); 
      for (t in 1:N_t) 
        CO2_t[t] = totalC_t0 - sum(C_t[t]); 
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
  simplex[3] A1[num_rep];       // output rates from pool 1
  simplex[3] A2[num_rep];       // output rates from pool 2
  simplex[3] A3[num_rep];       // output rates from pool 3
} 
transformed parameters { 
  vector<lower=0>[3] k;             // decomposition rates (1/turnover)
  matrix[N_t, num_rep] CO2_meas;    // evolved CO2 at measurement times
  matrix[N_t, num_rep] CO2_cap;     // evolved CO2 at cap times
  matrix[N_t, num_rep] CO2_flux_hat;// CO2 flux (average evolved CO2 between t_cap & t_meas)
  real<lower=0, upper=1> a21[num_rep];  // transfer rates 
  real<lower=0, upper=1> a31[num_rep];
  real<lower=0, upper=1> a12[num_rep];
  real<lower=0, upper=1> a32[num_rep];
  real<lower=0, upper=1> a13[num_rep];
  k = 1 ./ turnover;
  // transfer rates are different for each replication:
  for (i in 1:num_rep) {
    a21[i] = A1[i, 2];
    a31[i] = A1[i, 3];
    a12[i] = A2[i, 1];
    a32[i] = A2[i, 3];
    a13[i] = A3[i, 1];
  }
  for (i in 1:num_rep) {
    CO2_meas[,i] = evolved_CO2(N_t, t0, t_meas, gamma, totalC_t0, 
                         k, a21[i], a31[i], a12[i], a32[i], a13[i], 
                         x_r, x_i); 
    CO2_cap[,i] = evolved_CO2(N_t, t0, t_cap, gamma, totalC_t0, 
                        k, a21[i], a31[i], a12[i], a32[i], a13[i], 
                        x_r, x_i);
    CO2_flux_hat[,i] = (CO2_meas[,i] - CO2_cap[,i])./(t_meas - t_cap);
  }
} 
model { 
  // priors 
  turnover[1] ~ normal(1.5, 0.15 * sigma[1]);
  turnover[2] ~ normal(25, 2.5 * sigma[2]);
  turnover[3] ~ normal(1000, 100 * sigma[3]);
  sigma ~ cauchy(0,1); 
  sigma_obs ~ cauchy(0,.1);
  
  // likelihood     
  to_vector(CO2_flux) ~ lognormal(to_vector(log(CO2_flux_hat)), sigma_obs);   
} 
