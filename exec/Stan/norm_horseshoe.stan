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
  vector[N] Y;
  vector[P] mean_x;
  matrix[N,P] Q_x;
  matrix[P,P] R_inv_x;
  matrix[N,P] X;
  real m0;
  real scale_intercept;
}

transformed data {
  vector[N] Y_std;
  real sd_y;
  real mu_y;
  real m_scale = m0 / (P - m0);
  real inv_sqrt_N = 1/sqrt(1.0 * N);
  real slab_scale = 2;    // Scale for large slopes
  real slab_scale2 = square(slab_scale);
  real slab_df = 25;      // Effective degrees of freedom for large slopes
  real half_slab_df = 0.5 * slab_df;

  mu_y = mean(Y);
  sd_y = sd(Y);
  Y_std = (Y - mu_y)/sd_y;

}

parameters {
  vector[P] beta_tilde;
  real beta0_tilde;
  vector<lower=0>[P] sd_param;
  real<lower=0> tau_tilde;
  real<lower=0> c2_tilde;
  real<lower=0> sigma;
}

transformed parameters{
  vector[N] y_hat_std;
  vector[P] beta;
  real intercept;

  {
    real tau_0 = m_scale * sigma * inv_sqrt_N;
    real tau = tau_tilde * tau_0;
    real c2 = slab_scale2 * c2_tilde;
    vector[P] lambda = sqrt(c2) * sd_param ./ sqrt(c2 + square(tau) * square(sd_param));

    beta = beta_tilde .* lambda * tau;
    intercept = beta0_tilde * scale_intercept;

    y_hat_std = Q_x *  beta + intercept;
  }

}

model {
  beta_tilde ~ std_Normal();
  beta0_tilde ~ std_Normal_real();
  sd_param ~ std_cauchy();
  tau_tilde ~ std_cauchy_real();
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);
  sigma ~ std_Normal_real();

  Y_std ~ normal(y_hat_std, sigma);


}

generated quantities {
  vector[P + 1] theta;
  vector[N] y_hat;
  {
    vector[P] beta_temp;

    beta_temp = R_inv_x * beta;

    theta[1] = (intercept - (beta_temp)' * (mean_x ) ) * sd_y + mu_y ;

    for(i in 1 : (P)){
      theta[i+1] = beta_temp[i] * sd_y;
    }

    y_hat = y_hat_std * sd_y + mu_y;
  }

}
