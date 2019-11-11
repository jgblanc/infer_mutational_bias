functions {
  real my_raf_lpdf(real[] x, real W, int N, real mu, real B){
    vector[num_elements(x)] prob;
    real lprob;
    for (i in 1:num_elements(x)) {
      prob[i] = ((x[i] * (1- x[i]))^((4*N*mu) - 1)) * exp(-4 * N * B * W * x[i]) ;
    }
    lprob = sum(log(prob));
    return lprob;
  }
}
data {
  int<lower=0> m;         // Number snps
  real y[m];              // RAF
  real<lower=0> mu; // mutation rate
  real B; // effect size
  int N; // population size
}
parameters {
  real<lower=1e-20, upper=1e-3> W;                // selection
}
model {
  y ~ my_raf_lpdf(W, N, mu, B);       // prior log-density
}
