functions {
  real my_derived_risk_lpdf(real f, real W, int N, real mu, real B){
    real prob;
    prob = (exp(4 * N * B * W) - exp(4 * N * B * W * f)) / (-1 + exp(4 * N * B * W)) ;
    return prob;
  }
  real my_ll_risk_lpmf(int[] s, real[] x, real W, int N, real mu, real B){
    vector[num_elements(x)] prob;
    real lprob;
    for (i in 1:num_elements(x)) {
      if (x[i] == 0) {
        prob[i] = 1 - my_derived_risk_lpdf(x[i] | W, N, mu, B);
      }
      else {
        prob[i] = my_derived_risk_lpdf(x[i] | W, N, mu, B);
      }
    }
    lprob = sum(log(prob));
    return lprob;
  }
}
data {
  int<lower=0> m;         // Number snps
  real x[m];              // RAF
  int s[m];
  real<lower=0> mu; // mutation rate
  real B; // effect size
  int N; // population size
}
parameters {
  real<lower=1e-20, upper=1e-3> W;                // selection
}
model {
  s ~ my_ll_risk_lpmf(x, W, N, mu, B);
}
