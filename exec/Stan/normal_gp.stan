functions {
  matrix[] L_cov_exp_quad_ARD(matrix cost,
                            vector alpha,
                            vector rho,
                            real delta) { //squared exponential kernel
    int N = rows(cost);
    int P = rows(alpha);
    matrix[N, N] K;
    matrix[N, N] L[P];
    vector[P] sq_alpha;
    sq_alpha = alpha .* alpha;
    for (p in 1:P) {
      for (i in 1:(N-1)) {
        K[i, i] = sq_alpha[p] + delta;
        for (j in (i + 1):N) {
          K[i, j] = sq_alpha[p]
                        * exp(-0.5 * cost[i,j] / rho[p]);
          K[j, i] = K[i, j];
        }
      }
      K[N, N] = sq_alpha[p] + delta;
      L[p] = cholesky_decompose(K);
    }

    return L;
  }
  vector[] gen_mean(int N,
                    int P,
                    vector[] eta,
                    vector bias,
                    matrix cost,
                    vector alpha,
                    vector rho,
                    real delta ) { //generate mean function
    vector[P] f[N];
    vector[N] f_temp[P];
    matrix[N, N] L_K[P] = L_cov_exp_quad_ARD(cost, alpha, rho, delta);
    for(p in 1:P) f_temp[p] = L_K[p] * eta[p] + bias[p];
    for(p in 1:P) for(n in 1:N) f[n,p] = f_temp[p,n];

    return f;
  }
}
data {
  int<lower=1> N;
  int<lower=1> nT;
  int<lower=nT> N_total;
  int<lower=1> P;
  vector[P] y[N];
  matrix[N_total, N_total] cost;
  int time_idx[N];
}
transformed data {
  real delta = 1e-9; //positive term to keep pos def covar kernel
}
parameters {
  vector<lower=0>[P] rho;
  vector<lower=0>[P] alpha;
  vector[N_total] eta[P];
  vector[P] bias;
  cholesky_factor_corr[P] L_corr;
  vector<lower=0>[P] sds;
}
model {
  vector[P] mu[N]; //mean of data
  matrix[P,P] L; //cholesky of covariance of data

  //priors
  L_corr ~ lkj_corr_cholesky(2.0);
  rho ~ inv_gamma(5., 5.);
  alpha ~ std_normal();
  for(p in 1:P) eta[p] ~ std_normal();
  bias ~ std_normal();
  sds ~ std_normal();


  //suff stat for likelihood
  {
    vector[P] f[N_total] = gen_mean(N_total, P, eta, bias, cost, alpha, rho, delta);
    for(n in 1:N) mu[n] = f[time_idx[n]]; //replicate mean for various obs
  }
  L = diag_pre_multiply(sds, L_corr); //cholesky factor of covariance

  //likelihood
  y ~ multi_normal_cholesky(mu, L);
}

generated quantities {
  vector[P] pred_eta[N_total] = gen_mean(N_total, P, eta, bias, cost, alpha, rho, delta);

}
