data {
    int<lower = 0> N;
    int<lower = 0> N_obs;
    int<lower = 0> N_mis;
    int<lower = 0> ii_obs[N_obs];
    int<lower = 0> ii_mis[N_mis];
    real y_obs[N_obs];
    real c[N];
}
parameters {
    real y_mis[N_mis];
    vector[N] x;
    real<lower=0> sigma_p;
    real<lower=0> sigma_o;
    real x0;
    real beta;
}
transformed parameters {
    real y[N];
    y[ii_obs] = y_obs;
    y[ii_mis] = y_mis;
}
model {
    // priors
    sigma_o ~ cauchy(0, 3);
    sigma_p ~ cauchy(0, 3);
    x0 ~ normal(0, 10);
    beta ~ normal(0, 5);

    for(t in 1:N) {
        // state-model
        if(t == 1) {
            x[t] ~ normal(x0, sigma_p);
        } else {
            x[t] ~ normal(x[t-1] + beta * c[t-1], sigma_p);
        }
        // observation-model
        y[t] ~ normal(x[t], sigma_o);
    }
}
