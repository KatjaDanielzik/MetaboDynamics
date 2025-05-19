/// The input data is a vector 'y' of length 'N' ("log cpc") and a vector m of length M
// (metabolites)
data {
  int<lower=0> N; // total number of data points
  int<lower=0> M; // number of metabolites
  int<lower=0> t; // number of timesteps
  real y [N]; // outcome value: log(cpc) //each Metabolite*timepoints*3 replicates
  int<lower=0> X[N]; //time predictor
  int<lower=0> Me[N]; //metabolite predictor
}

// The parameters accepted by the model. Our model
// accepts three parameters 'mu', 'sigma' and lambda with mu [-inf;+inf] and sigma>0/lambda>0
parameters {
 real mu[M,t];  // metabolite and time specific mu
 real <lower=0> sigma [M,t]; // metabolite and time specific sigma
 real <lower=0> lambda[M]; //hyperprior for time specific sigma
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma' allowing for individual mus and sigmas per metabolite and timepoint
// we are pooling the metabolite specific timepoint sigmas to learn a regularizing prior lambda.
model {
  lambda~exponential(2);
  for (m in 1:M){
    sigma[m,]~exponential(lambda[m]);
    mu[m,]~normal(0,2); // we standardized to sd=1 and mu=0, so at
    //most all values at one predictor are 2 standard deviations away from 0
  }
for (i in 1:N){
    y[i]~normal(mu[Me[i],X[i]],sigma[Me[i],X[i]]);
  }
}


// model outputs additional to parameters
generated quantities {
  //loo
  real log_lik[N];
  // ppc
  real y_rep[M,t];
  // prior predictive check
  real y_prior;
  real mu_prior;
  real <lower=0> sigma_prior;
  real <lower=0> lambda_prior;
  //effect size:
  real eff_size[M,t-1]; //calculate effect size for every metabolite and comparison
                        // between timesteps
  //calcualte delta_mu between timesteps
  real delta_mu[M,t-1];

  //calculate posterior prediction
for (m in 1:M){
  for (i in 1:t){
   y_rep[m,i]= normal_rng(mu[m,i], sigma[m,i]);
   }
  }

   //Prior predictive check
   lambda_prior=abs(exponential_rng(2));
   sigma_prior=abs(exponential_rng(lambda_prior));
   mu_prior=normal_rng(0,2);
   y_prior=normal_rng(mu_prior,sigma_prior);

   //LOO
   for (n in 1:N){
   log_lik[n] = normal_lpdf(y[n]|mu[Me[n],X[n]],sigma[Me[n],X[n]]);
   }

   // calculate delta_mu
  for (m in 1:M){
    for (i in 2:t){
      delta_mu[m,i-1] = mu[m,i]-mu[m,i-1];
   }
  }

   // calculate effect size (Krutschke et al. 2013)
  for (m in 1:M){
    for (i in 2:t){
      eff_size[m,i-1] = (mu[m,i-1]-mu[m,i])/sqrt(((sigma[m,i-1]*sigma[m,i-1])+(sigma[m,i]*sigma[m,i]))/2);
   }
  }
}
