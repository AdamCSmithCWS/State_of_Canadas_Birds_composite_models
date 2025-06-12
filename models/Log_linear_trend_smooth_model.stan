//Measurement error log-linear trend

data {
  int<lower=0> n_years; // number of years of full time-series
  int<lower=0> n_indices; // number of years of full time-series
  vector[n_years] year_centered;
  array[n_indices] int<lower=0> year; // years with indices (allows for missing years)
  array[n_indices] real ln_index; // log scale annual indices of abundance or annual population size
  array[n_indices] real ln_index_sd; // SD of the log scale annual indices

}

parameters {
  real MU; // intercept

  real BETA;// slope
  real<lower=0> sigma;//
  vector[n_years] noise_raw;

}

transformed parameters {
  vector[n_years] mu; // vector of true annual indices after accounting for uncertainty
  vector[n_years] noise;

  noise = noise_raw*sigma;
  mu = (year_centered * BETA) + MU + noise;

}

model {
  MU ~ student_t(3,0,2);
  for(i in 1:n_indices){
  ln_index[year[i]] ~ normal(mu[year[i]], ln_index_sd[year[i]]);
  }
   //priors
  BETA ~ std_normal();
  noise_raw ~ std_normal();
  sigma ~ normal(0,0.1);
}

generated quantities {
vector[n_years] smooth_inds = (year_centered * BETA) + MU;


}
