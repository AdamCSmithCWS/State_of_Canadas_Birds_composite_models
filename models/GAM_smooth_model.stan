//Measurement error GAM

data {
  int<lower=0> n_years; // number of years of full time-series
  int<lower=0> n_indices; // number of years of full time-series
  array[n_indices] int<lower=0>  year; // years with indices (allows for missing years)
  array[n_indices] real ln_index; // log scale annual indices of abundance or annual population size
  array[n_indices] real ln_index_sd; // SD of the log scale annual indices

    // data for spline s(year)
  int<lower=1> n_knots_year;  // number of knots in the basis function for year
  matrix[n_years, n_knots_year] year_basis; // basis function matrix

}


parameters {
  real MU;

  vector[n_knots_year] BETA_raw;//_raw;
  real<lower=0> sdBETA;    // sd of GAM coefficients
  real<lower=0> sigma;//
  vector[n_years] noise_raw;

}

transformed parameters {
  vector[n_years] smooth_pred;
  vector[n_knots_year] BETA;
  vector[n_years] mu; // vector of true annual indices after accounting for uncertainty
  vector[n_years] noise;


  BETA = sdBETA*BETA_raw;
  noise = noise_raw*sigma;

  smooth_pred = year_basis * BETA;
  mu = MU + smooth_pred + noise;

}

model {
  MU ~ student_t(3,0,2);
  for(i in 1:n_indices){
  ln_index[i] ~ normal(mu[year[i]], ln_index_sd[i]);
  }
   //priors
  sdBETA ~ std_normal();
  BETA_raw ~ std_normal();
  noise_raw ~ std_normal();
  sigma ~ normal(0,0.1);
}

generated quantities {
vector[n_years] annual_diffs;
vector[n_years] scaled_status;
vector[n_years] scaled_log_status;
vector[n_years] mu_smooth = MU + smooth_pred;

annual_diffs[1] = 0;
scaled_status[1] = 0;
scaled_log_status[1] = 1;


for(y in 2:n_years){
  scaled_status[y] = mu_smooth[y]-mu_smooth[1];
  scaled_log_status[y] = exp(mu_smooth[y])/exp(mu_smooth[1]);
  annual_diffs[y] = mu_smooth[y]-mu_smooth[y-1];

}


}
