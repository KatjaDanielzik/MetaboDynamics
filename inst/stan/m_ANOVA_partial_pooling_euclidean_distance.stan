data {
  int<lower=0> N; // total number of data points
  int<lower=0> M; // number of metabolites
  int<lower=0> t; // number of timesteps
  int<lower=0> d; // number of doses (conditions)
  real y [N]; // outcome value: log(cpc)
  int<lower=1> Me[N]; //metabolite predictor
  int<lower=1> X[N]; //time predictor
  int<lower=1> Do[N]; //dose predictor
}

parameters {
  real mu[M,t,d];  // metabolite, time, dose
  real <lower=0> sigma[M,t,d]; // metabolite, time, dose, and cell line specific sigma
  real <lower=0> lambda[M,d]; // metabolite specific lambda (hyperprior for sigma)
}

model {
  for (m in 1:M) {
    lambda[m,] ~ exponential(2);
    for (i in 1:t) {
      for (j in 1:d) {
          sigma[m,i,j] ~ exponential(lambda[m,j]);
          mu[m,i,j] ~ normal(0, 2); // prior for mu
        }
      }
    }
  for (n in 1:N){
    y[n]~normal(mu[Me[n],X[n],Do[n]],sigma[Me[n],X[n],Do[n]]);
  }
}

generated quantities {
  real log_lik[N];
  real y_rep[N];
  real y_prior;
  real mu_prior;
  real <lower=0> sigma_prior;
  real <lower=0> lambda_prior;
  real delta_mu[M,d,t,t];
  real euclidean_distance[M,d,d]; # euclidean distance between metabolite and cell line specific longitudinal vectors of different doses 

  for (n in 1:N){
    y_rep[n]=normal_rng(mu[Me[n],X[n],Do[n]],sigma[Me[n],X[n],Do[n]]);
    log_lik[n] = normal_lpdf(y[n]|mu[Me[n],X[n],Do[n]],sigma[Me[n],X[n],Do[n]]);
  }

  lambda_prior = exponential_rng(2);
  sigma_prior = exponential_rng(lambda_prior);
  mu_prior = normal_rng(0, 2);
  y_prior = normal_rng(mu_prior, sigma_prior);


  // differences between time points per cell line and experimental condition
  for (m in 1:M) {
    for (j in 1:d) {
      for (t1 in 1:t) {
        for(t2 in 1:t){
          if(t1 < t2){
          delta_mu[m,j,t1,t2] = mu[m,t2,j] - mu[m,t1,j]; # t2 - t1 -> positive estimates mean increase
          }
          else {
          delta_mu[m,j,t1,t2] = 0;
        }
        }
      }
    }
  }
  

  // euclidean distances
  for (m in 1:M) {
      for (d1 in 1:d) {
        for (d2 in 1:d) {
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
