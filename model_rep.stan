// This model introduces variation within tumor tissue

data {
  
  // This part defines variables used in model section
  
  int<lower=1> N;   // number of observations
  int<lower=1> P;   // number of patients
  
  
  real y[N];   // transformed beta values
  int<lower=0, upper=1> tInd[N];   // tumor indicator
  
  
  int<lower=1, upper=P> pID[N];   // Patient ID
}



parameters {
  
  vector[P] b_pat_raw;
  vector[P] c_patT_raw;
  vector[N] d_T_raw;
  
  real betaT_raw;
  real mu_raw;
  
  real<lower=0,upper=1> sigma_e_raw; 
  real<lower=0,upper=1> sigma_p_raw;  
  real<lower=0,upper=1> sigma_pt_raw;
  real<lower=0,upper=1> sigma_t_raw;  
}


transformed parameters {
  
  
  vector[P] b_pat;   /// random  effect of patients
  vector[P] c_patT;  // random effect of tumor tissue (varying between patients)
  vector[N] d_T; // random effect of tumor tissue (intratumor)
  real betaT;   // fixed effect of tumor tissue
  real mu;    // grand mean
  
  real<lower=0> sigma_e;  // error sd
  real<lower=0> sigma_p;  // patient sd
  real<lower=0> sigma_pt; // tissue sd at patient level
  real<lower=0> sigma_t;  // intratumor sd
  
  
  
  mu <- 100000 * mu_raw;
  betaT <- 100000 * betaT_raw;
  
  sigma_e <- 100 * sigma_e_raw;
  sigma_p <- 100 * sigma_p_raw;
  sigma_pt <- 100 * sigma_pt_raw;
  sigma_t <- 100 * sigma_t_raw;
  
  b_pat <- sigma_p * b_pat_raw;
  c_patT <- sigma_pt * c_patT_raw;
  d_T <- sigma_t * d_T_raw;
}



model {
  //prior
    
  mu_raw ~ normal(0,1);
  betaT_raw ~ normal(0,1);
  
  sigma_t_raw ~ uniform(0,1);
  sigma_e_raw ~ uniform(0,1);
  sigma_p_raw ~ uniform(0,1);
  sigma_pt_raw ~ uniform(0,1);

  b_pat_raw ~ normal(0,1);
  c_patT_raw ~ normal(0,1);
  d_T_raw ~ normal(0,1);
  

  
  
  //posterior
  
  for (n in 1:N){
    if (tInd[n]==0) {
      y[n] ~ normal(b_pat[pID[n]] + mu, sigma_e);
    } else {
      y[n] ~ normal(b_pat[pID[n]] + betaT + c_patT[pID[n]] + d_T[n] + mu, 
      sigma_e);
    }
  }
  
}
