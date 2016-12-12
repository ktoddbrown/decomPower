functions { 
 
  /** 
   * ODE system for thr Century model with no inputs. 
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
   * parameters and times.  This is done by simulating the system 
   * defined by the ODE function century_model and then 
   * calculating the rate CO2 is emmited
   * 
   * @param T number of times. 
   * @param t0 initial time. 
   * @param ts observation times. 
   * @param gamma partitioning coefficient. 
   * @param k decomposition rates  
   * @param ajk transfer rates 
   * @param x_r real data (empty) 
   * @param x_i integer data (empty) 
   * @return evolved CO2 for times ts 
   */ 
  vector evolved_CO2(int N_t, real t0, real[] ts, 
                     vector gamma, real totalC_t0, 
                     vector k, real a21, real a31, real a12, 
                     real a32, real a13, real[] x_r, int[] x_i) { 
 
    real C_t0[3];               // initial state 
    real theta[8];              // ODE parameters 
    real C_hat[N_t,3];          // predicted pool content 
    vector[N_t] eCO2_hat; 
 
    C_t0 = to_array_1d(gamma*totalC_t0);
    theta[1:3] = to_array_1d(k);
    theta[4] = a21;
    theta[5] = a31;
    theta[6] = a12;
    theta[7] = a32;
    theta[8] = a13;
    C_hat = integrate_ode_bdf(century_model,  
                           C_t0, t0, ts, theta, x_r, x_i); 
    for (t in 1:N_t) 
      eCO2_hat[t] = totalC_t0 - sum(C_hat[t]); 
    return eCO2_hat; 
  } 
} 
data { 
  real<lower=0> totalC_t0;     // initial total carbon 
  real t0;                     // initial time 
  int<lower=0> N_t;            // number of measurement times 
  real<lower=t0> ts[N_t];      // measurement times 
  real<lower=0> eCO2[N_t];     // measured cumulative evolved carbon 
} 
transformed data { 
  real x_r[0];                 // no real data for ODE system 
  int x_i[0];                  // no integer data for ODE system 
} 
parameters { 
  vector<lower=0>[3] turnover; // normalized turnover rates
  real<lower=0, upper=1> a21;  // transfer rates
  real<lower=0, upper=1> a31;
  real<lower=0, upper=1> a12;
  real<lower=0, upper=1> a32;
  real<lower=0, upper=1> a13;
  simplex[3] gamma;           // partitioning coefficients (a simplex) 
  vector<lower=0>[3] sigma;   // observation standard deviation 
  real<lower=0> sigma_obs;
} 
transformed parameters { 
  vector[N_t] eCO2_hat; 
  vector<lower=0>[3] k;
  k = 1 ./ turnover;
  eCO2_hat = evolved_CO2(N_t, t0, ts, gamma, totalC_t0, 
                          k, a21, a31, a12, a32, a13, 
                          x_r, x_i); 
} 
model { 
  // priors 
  turnover[1] ~ normal(1.5, 0.15 * sigma[1]);
  turnover[2] ~ normal(25, 2.5 * sigma[2]);
  turnover[3] ~ normal(1000, 100 * sigma[3]);
  sigma ~ cauchy(0,1); 
  sigma_obs ~ cauchy(0,1);
 
  // likelihood     
  eCO2 ~ lognormal(log(eCO2_hat), sigma_obs);   // lognormal error 
} 