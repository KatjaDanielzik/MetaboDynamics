data {
  int <lower=0 > N; // number of data points maven values
  real maven[N]; // LCMS measurements processed by maven software
  int <lower=1 > M; // number of metabolites
  int <lower=1 > t; // number of time points
  int <lower=1 > D; // number of radiation doses (experimental conditions)
  int <lower=1 > Me[N]; // metabolite predictor
  int <lower=1 > X[N]; // time predictor
  int <lower=1 > Do[N]; // dose predictor
  int <lower=0 > Nc; // number of data points cell counts
  int <lower=0 > Cc[Nc]; // cell counts
  int <lower=1 > X_c[Nc]; // time predictor for cell counts
  int <lower=1 > Do_c[Nc]; // dose predictor for cell counts
  
  # Priors
  real<lower=0> prior_mean_abundance[2]; // prior for mu_maven (normal distribution): mu and sigma, expected range for mean of metabolite abundances
  real<lower=0> prior_sd_abundance; // expected value for standard deviation of metabolite abundances per metabolite/condition
  real<lower=0> prior_counts; // expected value of cell counts (rate of 
                                // exponential distribution) -> will be divided by one  
}

parameters {
  real mu_maven[M, t, D]; // mean of maven values for each metabolite, time, dose, and cell line
  real <lower=0> sd_maven[M, D]; // standard deviation of maven values for each metabolite, and dose
  real <lower=0> lambda_maven[M]; // metabolite specific lambda (hyperprior for sigma)
  
  real <lower=0> mu_counts[t, D]; // mean of cell counts for each time, dose, and cell line
}

transformed parameters {
  real mean_natural[M,t,D];
  real cpc[M, t, D];
  real log_cpc[M, t, D];
  for (m in 1:M) {
    for (i in 1:t) {
      for (d in 1:D) {
          mean_natural[m,i,d] = exp(mu_maven[m, i, d]+(sd_maven[m, d]/2)); # expected value for log-normal distribution 2= exp(mu+(sigma²/2))
          cpc[m, i, d] = mean_natural[m,i,d]/ mu_counts[i, d];
          log_cpc[m, i, d] = log10(cpc[m, i, d]);
      }
    }
  }

  // scaling of log_cpc with z-transformation?
  // Calculate mean and standard deviation for z-transformation
  real mean_log_cpc[M, D];
  real sd_log_cpc[M, D];
  for (m in 1:M) {
    for (d in 1:D) {
        mean_log_cpc[m, d] = mean(log_cpc[m, :, d]);
        sd_log_cpc[m, d] = sd(log_cpc[m, :, d]);
    }
  }

  // Apply z-transformation
  real mu[M,t,D];
  for (m in 1:M) {
    for (d in 1:D) {
      for (i in 1:t) {
          mu[m, i, d] = (log_cpc[m, i, d] - mean_log_cpc[m, d]) / sd_log_cpc[m, d];
      }
    }
  }

}

model {
  // priors
  // maven
  for (m in 1:M) {
    lambda_maven[m] ~ exponential(1/prior_sd_abundance);
      for (d in 1:D) {
        for (i in 1:t) {
          sd_maven[m, d] ~ exponential(lambda_maven[m]); // hierarchy of sd_maven: pooling of sd for all measurements of one metabolite
          mu_maven[m, i, d] ~ normal(prior_mean_abundance[1], prior_mean_abundance[2]);
      }
    }
  }

  // cell counts
  for (i in 1:t) {
    for (d in 1:D) {
        mu_counts[i,d]~exponential(1/prior_counts); // expects values between 0 and 2e6 cells (highest probability mass)
    }
  }

  // single models
  // LC-MS
  for (n in 1:N) {
    maven[n] ~ lognormal(mu_maven[Me[n], X[n], Do[n]], sd_maven[Me[n], Do[n]]);
  }

  // cell counts
   for (n in 1:Nc) {
    Cc[n] ~ poisson(mu_counts[X_c[n], Do_c[n]]);
  }
}

generated quantities {
  // yrep and log-lik
  real maven_rep[N];
  real counts_rep[Nc];
  real log_lik[N + Nc];

  // differences between time points and euclidean distances between z-scaled vectors
  real delta_mu[M,D,t,t];
  real euclidean_distance[M,D,D]; # euclidean distance between metabolite and cell line specific longitudinal vectors of different doses

  // Prior predictive check
  real mu_counts_prior = exponential_rng(1/prior_counts);
  real counts_prior = poisson_rng(mu_counts_prior);

  real lambda_maven_prior = exponential_rng(1/prior_sd_abundance);
  real sigma_maven_prior = exponential_rng(lambda_maven_prior);
  real mu_maven_prior = normal_rng(prior_mean_abundance[1],prior_mean_abundance[2]);
  real maven_prior = lognormal_rng(mu_maven_prior,sigma_maven_prior);


  // y_rep
  for (n in 1:N) {
    maven_rep[n] = lognormal_rng(mu_maven[Me[n], X[n], Do[n]], sd_maven[Me[n], Do[n]]); // back transformation of maven values
  }


  for (n in 1:Nc) {
    counts_rep[n] = poisson_rng(mu_counts[X_c[n], Do_c[n]]);
  }

  // log_lik
  for (n in 1:N) {
    log_lik[n] = lognormal_lpdf(maven[n] | mu_maven[Me[n], X[n], Do[n]], sd_maven[Me[n], Do[n]]);
  }

  for (n in 1:Nc) {
    log_lik[N + n] = poisson_lpmf(Cc[n] | mu_counts[X_c[n], Do_c[n]]);
  }


  // differences between time points per cell line and experimental condition
  for (m in 1:M) {
    for (j in 1:D) {
      for (t1 in 1:t) {
        for(t2  in 1:t){
          if(t1 < t2){
          delta_mu[m,j,t1,t2] = mu[m,t2,j] - mu[m,t1,j];
          }
        else {
          delta_mu[m,j,t1,t2] = 0;
        }
        }
      }
    }
  }


  // euclidean distances between scaled dose specific vectors
  for (m in 1:M) {
    for (d1 in 1:D) {
      for (d2 in 1:D) {
        if (d1 < d2) {
            vector[t] mu_d1;
            vector[t] mu_d2;
            mu_d1 = to_vector(mu[m,:,d1]); # time point length vectors per metabolite dose and cell line
            mu_d2 = to_vector(mu[m,:,d2]);
            euclidean_distance[m,d1,d2] = distance(mu_d1,mu_d2); # euclidean distance between vectors
          } else {
            euclidean_distance[m,d1,d2] = 0;
          }
        }
      }
    }
  }
