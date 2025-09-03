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
}

parameters {
  real mu_maven[M, t, D]; // mean of maven values for each metabolite, time, dose, and cell line
  real <lower=0 > sd_maven[M, t, D]; // standard deviation of maven values for each metabolite, time, dose, and cell line
  real <lower=0> lambda_maven[M, D]; // metabolite specific lambda (hyperprior for sigma)
  
  real <lower=0 > mu_counts[t, D]; // mean of cell counts for each time, dose, and cell line
}

transformed parameters {
  real mean_natural[M,t,D];
  real cpc[M, t, D];
  real log_cpc[M, t, D];
  for (m in 1:M) {
    for (i in 1:t) {
      for (d in 1:D) {
          mean_natural[m,i,d] = exp(mu_maven[m, i, d]+(sd_maven[m, i, d]/2)); # expected value for log-normal distribution 2= exp(mu+(sigma²/2))
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
    lambda_maven[m,] ~ exponential(2);
      for (d in 1:D) {
        for (i in 1:t) {
          sd_maven[m, i, d] ~ exponential(lambda_maven[m,d]); // hierarchy of sd_maven: pooling of sd for all measurements of one metabolite
          mu_maven[m, i, d] ~ normal(12, 5);
      }
    }
  }

  // cell counts
  for (i in 1:t) {
    for (d in 1:D) {
        mu_counts[i,d]~exponential(1/1e6); // expects values between 0 and 2e6 cells (highest probability mass)
    }
  }

  // single models
  // LC-MS
  for (n in 1:N) {
    maven[n] ~ lognormal(mu_maven[Me[n], X[n], Do[n]], sd_maven[Me[n], X[n], Do[n]]);
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
  real delta_mu[M,t-1,D];
  real euclidean_distance[M,D,D]; # euclidean distance between metabolite and cell line specific longitudinal vectors of different doses

  // Prior predictive check
  real mu_counts_prior = exponential_rng(1/1e6);
  real counts_prior = poisson_rng(mu_counts_prior);
  
  real lambda_maven_prior = exponential_rng(2);
  real sigma_maven_prior = exponential_rng(lambda_maven_prior);
  real mu_maven_prior = normal_rng(12,5);
  real maven_prior = normal_rng(mu_maven_prior,sigma_maven_prior);


  // y_rep
  for (n in 1:N) {
    maven_rep[n] = lognormal_rng(mu_maven[Me[n], X[n], Do[n]], sd_maven[Me[n], X[n], Do[n]]); // back transformation of maven values
  }


  for (n in 1:Nc) {
    counts_rep[n] = poisson_rng(mu_counts[X_c[n], Do_c[n]]);
  }

  // log_lik
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(log10(maven[n]) | mu_maven[Me[n], X[n], Do[n]], sd_maven[Me[n], X[n], Do[n]]);
  }

  for (n in 1:Nc) {
    log_lik[N + n] = poisson_lpmf(Cc[n] | mu_counts[X_c[n], Do_c[n]]);
  }


  // delta mu between time points at one radiation dose of one metabolite
  for (m in 1:M) {
    for (i in 2:t) {
      for (d in 1:D) {
          delta_mu[m,i-1,d] = mu[m,i,d] - mu[m,i-1,d];
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
