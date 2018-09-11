// - county RE
// - Bernoulli on including county RE
// - error in variables
data {
    int<lower=0> N; // number of observations / locations
    int<lower=0> p; // number of covariates
    matrix[N, p] X; // covariate matrix, w/o intercept
    vector<lower=0>[N] y; // response
    // matrix[N, N] W; // spatial weights matrix (non-sparse version)
    int Wn; // num. non-zero entries in spatial weights
    int Wi[Wn]; // i index of spatial weights
    int Wj[Wn]; // j index of spatial weights
    vector[Wn] Ww; // spatial weight
    
}
parameters {
    real<lower = -1, upper = 1> rho; // spatial autoregression coef.
    real alpha; // intercept
    vector[p] beta; // covariate effects
    vector[p] gamma; // lagged covariate effects
    real<lower = 0> sigma; // sd of noise
}
transformed parameters {
    vector[N] y_hat;
    
    // Non-sparse version
    // y_hat = rho*W*y + alpha*rep_vector(1, N) + X*beta + W*X*gamma;
    // Sparse version
    {
        vector[N] y_hat_local;
        vector[N] y_hat_lags;
        
        y_hat_local = alpha*rep_vector(1, N) + X*beta;
        y_hat_lags = rep_vector(0, N);
        for (w_idx in 1:Wn) {
            y_hat_lags[Wi[w_idx]] += rho*Ww[w_idx]*y[Wj[w_idx]];
            for (col_idx in 1:p) {
            y_hat_lags[Wi[w_idx]] += Ww[w_idx]*X[Wj[w_idx], col_idx]*gamma[col_idx];
            }
        }
        
        y_hat = y_hat_local + y_hat_lags;
    }
}
model {
    y ~ normal(y_hat, sigma);
    
    // priors
    alpha ~ normal(0, 3);
    beta ~ normal(0, 3);
    gamma ~ normal(0, 3);
    rho ~ normal(0, .5);
    sigma ~ cauchy(0, 3);
}
generated quantities {
}
