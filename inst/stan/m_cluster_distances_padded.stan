
data {
  int <lower=1> C; // number of cluster comparisons
  array[C] int<lower=1> N; //number of observations
  int <lower=0> M; //maximum number of observations
  matrix [C,M] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  array[C] real<lower=0> mu;
  array[C] real<lower=0> sigma;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
mu~normal(2,2);
sigma~exponential(1);
  for (c in 1:C){
    y[c,1:N[c]]~normal(mu[c],sigma[c]);
  }
}
// 
// generated quantities{
//   matrix [C,M] y_rep;
//     for(c in 1:C){
//       for(n in 1:M){
//     y_rep[c,n]=normal_rng(mu[c],sigma[c]);
//       }
//     }
// }
