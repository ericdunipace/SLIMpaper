functions{
  real std_cauchy_lpdf(vector Y) {
    return - sum(log1p(square(Y)));
  }
  real std_cauchy_real_lpdf(real Y) {
    return (- log1p(square(Y)));
  }
  real std_Normal_lpdf(vector Y) {
    return - 0.5 * dot_self(Y);
  }
  real std_Normal_real_lpdf(real Y) {
    return - 0.5 * Y*Y;
  }
}
data {
  int N;
  int P;
  int Y[N];
  matrix[N,P] X;
  real m0;
  real scale_intercept;
}

transformed data {
  matrix[N,P] X_std;
  vector[P] mean_x;
  vector[P] sd_x;
  real m_scale = m0 / (P - m0);
  real inv_sqrt_N = 1/sqrt(1.0 * N);
  real slab_scale = 2;    // Scale for large slopes
  real slab_scale2 = square(slab_scale);
  real slab_df = 25;      // Effective degrees of freedom for large slopes
  real half_slab_df = 0.5 * slab_df;

  {
    matrix[N,P] mu_x;
    matrix[N,P] sigma_x;

    for(p in 1:P) {
      mean_x[p] = mean(col(X,p));
      sd_x[p] = sd(col(X,p));
      for(n in 1:N) {
        mu_x[n,p] =  mean_x[p];
        sigma_x[n,p] = sd_x[p];
      }
    }
    X_std = (X - mu_x) ./ sigma_x;
  }
  // Q_x = qr_Q(X_std)[, 1:P] * sqrt(N - 1);
  // R_x = qr_R(X_std)[1:P, ] / sqrt(N - 1);
  // R_inv_x = inverse(R_x);

}

parameters {
  vector[P] beta_tilde;
  real beta0_tilde;
  vector<lower=0>[P] sd_param;
  real<lower=0> tau_tilde;
  real<lower=0> c2_tilde;
}

transformed parameters{
  vector[N] eta;
  vector[P] beta;
  real intercept;


  {
    real tau_0 = m_scale * inv_sqrt_N;
    real tau = tau_tilde * tau_0;
    real c2 = slab_scale2 * c2_tilde;
    vector[P] lambda = sqrt(c2) * sd_param ./ sqrt(c2 + square(tau) * square(sd_param));

    beta = beta_tilde .* lambda * tau;
    intercept = beta0_tilde * scale_intercept;

    eta = X_std * beta + intercept;
  }

}

model {
  beta_tilde ~ std_Normal();
  beta0_tilde ~ std_Normal_real();
  sd_param ~ std_cauchy();
  tau_tilde ~ std_cauchy_real();
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);

  Y ~ bernoulli_logit(eta);


}

generated quantities {
  vector[P + 1] theta;
  vector[N] prob;
  {
    vector[P] beta_temp;

    beta_temp = beta;

    theta[1] = intercept  - (beta_temp)' * (mean_x ./ sd_x);

    for(i in 1 : (P)){
      theta[i+1] = beta_temp[i]/sd_x[i];
    }
    // eta = theta[1] + X * beta_temp ./ sd_x
    prob = inv_logit(eta);
  }

}
