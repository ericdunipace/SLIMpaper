// adappted from https://mc-stan.org/docs/2_18/stan-users-guide/fit-gp-section.html
// functions {
//   matrix L_cov_exp_quad_ARD(vector[] x,
//                             real alpha,
//                             vector rho,
//                             real delta) {
//     int N = size(x);
//     matrix[N, N] K;
//     real sq_alpha = square(alpha);
//     for (i in 1:(N-1)) {
//       K[i, i] = sq_alpha + delta;
//       for (j in (i + 1):N) {
//         K[i, j] = sq_alpha
//                       * exp(-0.5 * dot_self((x[i] - x[j]) ./ rho));
//         K[j, i] = K[i, j];
//       }
//     }
//     K[N, N] = sq_alpha + delta;
//     return cholesky_decompose(K);
//   }
// }
// data {
//   int<lower=1> N;
//   int<lower = 0> N_test;
//   int<lower=1> P;
//   vector[P] x[N];
//   int y[N];
//   vector[P] x_test[N_test];
// }
// transformed data {
//   real delta = 1e-9;
//   int<lower = N> N_total = N+N_test;
//   vector[P] X[N_total];
//
//   for(i in 1:N) X[i] = x[i];
//   if(N_total > N) {
//     for(i in 1:N_test) X[i + N] = x_test[i];
//   }
// }
// parameters {
//   vector<lower=0>[P] rho;
//   real<lower=0> alpha;
//   vector[N_total] eta;
//   real bias;
// }
// model {
//   vector[N_total] f;
//   {
//     matrix[N_total, N_total] L_K = L_cov_exp_quad_ARD(X, alpha, rho, delta);
//     f = L_K * eta;
//   }
//
//   rho ~ inv_gamma(5., 5.);
//   alpha ~ std_normal();
//   sigma ~ std_normal();
//   eta ~ std_normal();
//   bias ~ std_normal();
//
//   y ~ bernouli_logit(bias + f[1:N]);
// }
//
// generated quantities {
//   vector[N_total] pred_eta;
//
//   {
//     matrix[N_total, N_total] L_K = L_cov_exp_quad_ARD(X, alpha, rho, delta);
//     pred_eta = L_K * eta + bias;
//   }
//
//
// }

functions {
  matrix L_cov_exp_quad_ARD(vector[] x,
                            real alpha,
                            vector rho,
                            real delta) {
    int N = size(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha
                      * exp(-0.5 * dot_self((x[i] - x[j]) ./ rho));
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha + delta;
    return cholesky_decompose(K);
  }
}
data {
  int<lower=1> N;
  int<lower=0> N_test;
  int<lower=1> P;
  vector[P] x[N];
  int y[N];
  vector[P] x_test[N_test];
}
transformed data {
  real delta = 1e-9;
  int<lower = N> N_total = N + N_test;
  vector[P] X[N_total];

  for(i in 1:N) X[i] = x[i];
  if(N_total > N) {
    for(i in 1:N_test) X[i + N] = x_test[i];
  }
}
parameters {
  vector<lower=0>[P] rho;
  real<lower=0> alpha;
  vector[N_total] eta;
  real bias;
}
model {
  vector[N_total] f;
  {
    matrix[N_total, N_total] L_K = L_cov_exp_quad_ARD(X, alpha, rho, delta);
    f = L_K * eta;
  }

  rho ~ inv_gamma(5, 5);
  alpha ~ std_normal();
  eta ~ std_normal();
  bias ~ std_normal();

  y ~ bernoulli_logit(bias + f[1:N]);
}

generated quantities {
  vector[N_total] pred_eta;

  {
    matrix[N_total, N_total] L_K = L_cov_exp_quad_ARD(X, alpha, rho, delta);
    vector[N_total] f = L_K * eta;
    pred_eta = f + bias;
  }

}
