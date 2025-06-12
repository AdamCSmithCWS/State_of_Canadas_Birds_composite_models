//All years estimated independently for State of Birds model

data {
  int<lower=0> n_years; // number of years of full time-series
  int<lower=0> n_species; // number of species in the group
  array[n_years] int n_species_year; // vector of number of species to include in each year
  array[n_years,n_species] int species; // matrix of species indicators included in each year allows for missing species
  array[n_years,n_species] real ln_index; // log scale annual indices of abundance or annual population size
  array[n_years,n_species] real ln_index_sd; // SD of the log scale annual indices

}


parameters {
  //real<lower=0> sigma;
  array[n_years] real MU;// mean annual status values
  array[n_years] real<lower=0> sigma;    // sd each yearly summary
  array[n_years,n_species] real noise_raw;
}

transformed parameters{
    array[n_years,n_species] real ln_index_true; // true index value for each species and year

    for(i in 1:n_years){

      for(s in 1:n_species){

    ln_index_true[i,s] = MU[i] + (sigma[i])*noise_raw[i,s];


}
}

}
model {
  MU ~ student_t(3,0,1);

  for(i in 1:n_years){
 noise_raw[i,] ~ normal(0,1);

   for(s in species[i,1:n_species_year[i]]){ // stepping through species with data
  ln_index[i,s] ~ normal(ln_index_true[i,s], ln_index_sd[i,s]);
    }

  }

  sigma ~ student_t(3,0,1);
}

 generated quantities {
 vector[n_years+1] annual_status;

 annual_status[1] = 0;

 for(y in 1:n_years){
   annual_status[y+1] = MU[y];

 }
}

