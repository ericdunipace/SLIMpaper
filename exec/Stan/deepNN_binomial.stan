functions{
  /** from https://discourse.mc-stan.org/t/calculating-layers-in-a-neural-net/1399/10
   * @param x Predictors (N x M)
   * @param alpha First-layer weights (M x J)
   * @param beta Second-layer weights (J x (K - 1))
   * @return Linear predictor for output layer of RBM.
   */
  // matrix elu(matrix x, matrix beta, matrix alpha) {
  //   return expm1(x* beta + alpha);
  // }
  matrix sigmoid(matrix x, matrix beta, matrix alpha) {
    return inv_logit(x * beta + alpha);
  }
  matrix softplus(matrix x, matrix beta, matrix alpha) {
    return log1p_exp(x * beta + alpha);
  }
  matrix tanhact(matrix x, matrix beta, matrix alpha) {
    return tanh(x* beta + alpha);
  }

  matrix activate(matrix x, matrix beta, matrix alpha){
    return softplus(x,  beta,  alpha);
  }

  vector activated3(matrix x, matrix input_beta, vector input_alpha,
  matrix[] beta, vector[] alpha, vector output_beta, real output_bias,
  int L, int nodes, int N, int P) {
    matrix[N, nodes] in_alpha;
    matrix[N, nodes] hi_alpha[L];
    vector[N] out_bi_vec;

    if(L != 2) print("L is must be too");
    for(nn in 1:nodes) for(n in 1:N) in_alpha[n,nn] = input_alpha[nn];
    for(l in 1:L){
      for(nn in 1:nodes) for(n in 1:N) hi_alpha[l,n,nn] = alpha[L,nn];
    }
    for(n in 1:N) out_bi_vec[n] = output_bias;

    return softplus(softplus(softplus(x, input_beta, in_alpha), beta[1], hi_alpha[1]), beta[2], hi_alpha[2]) * output_beta + out_bi_vec;
  }

  vector activated2(matrix x, matrix input_beta, vector input_alpha,
  matrix[] beta, vector[] alpha, vector output_beta, real output_bias,
  int L, int nodes, int N, int P) {
    matrix[N, nodes] in_alpha;
    matrix[N, nodes] hi_alpha[L];
    vector[N] out_bi_vec;

    if(L != 1) print("L is must be t2o");
    for(nn in 1:nodes) for(n in 1:N) in_alpha[n,nn] = input_alpha[nn];
    for(l in 1:L){
      for(nn in 1:nodes) for(n in 1:N) hi_alpha[l,n,nn] = alpha[L,nn];
    }
    for(n in 1:N) out_bi_vec[n] = output_bias;

    return softplus(softplus(x, input_beta, in_alpha), beta[1], hi_alpha[1]) * output_beta + out_bi_vec;
  }

  matrix activated(matrix x, matrix input_beta, vector input_alpha, matrix[] beta, vector[] alpha, int L, int nodes, int N, int P);

  matrix activated(matrix x, matrix input_beta, vector input_alpha, matrix[] beta, vector[] alpha, int L, int nodes, int N, int P){

    if(L == 0) {
      matrix[N,nodes] in_alpha;
    for(nn in 1:nodes) for(n in 1:N) in_alpha[n,nn] = input_alpha[nn];
      return activate(x,  input_beta,  in_alpha);
    } else {
      matrix[N, nodes] alpha_mat;
      for(nn in 1:nodes) for(n in 1:N) alpha_mat[n,nn] = alpha[L,nn];
      return activate(activated(x, input_beta, input_alpha, beta, alpha, L-1, nodes, N, P),beta[L], alpha_mat);
    }
  }
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
  int N; //number observations
  int P; //number features
  int L; //number hidden layers
  real m0; //number non-zero param
  int nodes; // number nodes
  matrix[N,P] X; // design matrix
  int Y[N]; //outcome
}

transformed data{
  // int LL;
  // real M = P * nodes + 2 * nodes + 1;
  // real slab_scale = 3;    // Scale for large slopes
  // real slab_scale2 = square(slab_scale);
  // real slab_df = 25;      // Effective degrees of freedom for large slopes
  // real half_slab_df = 0.5 * slab_df;
  //
  // LL = L - 1;

}

parameters{
  // matrix[P,nodes] input_weights_tilde;
  // vector[nodes] output_weights_tilde;
  // vector[nodes] input_bias_tilde;
  // real output_bias_tilde;
  // //scales
  // real<lower=0> tau_tilde;
  // real<lower=0> tau_lam;
  // real<lower=0> c2_tilde;
  //
  // vector<lower=0>[nodes] input_bias_sigma;
  // vector<lower=0>[nodes] ib_lam;
  // matrix<lower=0>[P,nodes] input_weights_sigma;
  // matrix<lower=0>[P,nodes] iw_lam;
  // vector<lower=0>[nodes] output_weights_sigma;
  // vector<lower=0>[nodes] ow_lam;
  // real<lower=0> output_bias_sigma;
  // real<lower=0> ob_lam;
  //
  // matrix[nodes, nodes] hidden_weights_tilde[LL];
  // vector[nodes] hidden_bias_tilde[LL];
  //
  // vector<lower=0>[nodes] hidden_bias_sigma[LL];
  // vector<lower=0>[nodes] hb_lam[LL];
  // matrix<lower=0>[nodes, nodes] hw_lam[LL];
  // matrix<lower=0>[nodes, nodes] hidden_weights_sigma[LL];
}
transformed parameters {
  // vector[N] eta;
  // {
  //   vector[nodes] output_weights;
  //   matrix[P, nodes] input_weights;
  //   matrix[nodes,nodes] hidden_weights;
  //   real output_bias;
  //   vector[nodes] hidden_bias;
  //   vector[nodes] input_bias;
  //   real tau0 = (m0 / (M - m0)) * (1.0 / sqrt(1.0 * N));
  //   real tau = tau0 * tau_tilde; // tau ~ cauchy(0, tau0)
  //   // matrix[P,nodes] input_weights_sigma = sqrt(input_weights_sigma);
  //   // vector[nodes] output_weights_sigma = sqrt(output_weights_sigma);
  //
  //   // c2 ~ inv_gamma(half_slab_df, half_slab_df * slab_scale2)
  //   // Implies that marginally beta ~ student_t(slab_df, 0, slab_scale)
  //   // real c2 = slab_scale2 * c2_tilde;
  //   real c2 = c2_tilde;
  //
  //   // vector[nodes] lambda_output =
  //   //   sqrt( c2 * square(output_weights_sigma) ./ (c2 + square(tau) * square(output_weights_sigma)) );
  //   //
  //   // matrix[P, nodes] lambda_input = sqrt( c2 * square(input_weights_sigma) ./ (c2 + square(tau) * square(input_weights_sigma)) );
  //
  //   // vector[nodes] lambda_output =
  //   //   sqrt( c2 * output_weights_sigma ./ (c2 + tau * output_weights_sigma) );
  //   // vector[nodes] lambda_output = sqrt(output_weights_sigma);
  //   // real lambda_output_bias = sqrt( c2 * output_bias_sigma ./ (c2 + tau * output_bias_sigma) );
  //
  //   matrix[P, nodes] lambda_input = sqrt( c2 * input_weights_sigma ./ (c2 + tau * input_weights_sigma) );
  //   vector[nodes] lambda_input_bias = sqrt( c2 * input_bias_sigma ./ (c2 + tau * input_bias_sigma) );
  //
  //   matrix[nodes, nodes] lambda_hidden;
  //   vector[nodes] lambda_hidden_bias;
  //   matrix[N, nodes] activation[2];
  //   // matrix[N, nodes] activation;
  //   matrix[N, nodes] bias;
  //   vector[N] out_bi_vec;
  //
  //   // weights ~ normal(0, tau * lambda_tilde)
  //   input_weights = sqrt(tau) * lambda_input .* input_weights_tilde;
  //   input_bias = sqrt(tau) * lambda_input_bias .* input_bias_tilde;
  //   // output_weights = sqrt(tau) *lambda_output .* output_weights_tilde;
  //   // output_bias = sqrt(tau) *lambda_output_bias .* output_bias_tilde;
  //   output_weights = output_weights_tilde * 5.0;
  //   output_bias = output_bias_tilde * 5.0;
  //
  //   for(nn in 1:nodes) for(n in 1:N) bias[n,nn] = input_bias[nn];
  //   activation[2] = activate(X, input_weights, bias);
  //
  //   for(l in 1:LL) {
  //     // lambda_hidden = sqrt( c2 * square(hidden_weights_sigma[l]) ./ (c2 + square(tau) * square(hidden_weights_sigma[l])) );
  //     lambda_hidden = sqrt(c2 * hidden_weights_sigma[l] ./ (c2 + tau * hidden_weights_sigma[l]) );
  //     lambda_hidden_bias = sqrt( c2 * hidden_bias_sigma[l] ./ (c2 + tau * hidden_bias_sigma[l]) );
  //
  //     hidden_weights = sqrt(tau) * lambda_hidden .* hidden_weights_tilde[l];
  //     hidden_bias = sqrt(tau) * lambda_hidden_bias .* hidden_bias_tilde[l];
  //     activation[1] = activation[2];
  //
  //     for(nn in 1:nodes) for(n in 1:N) bias[n,nn] = hidden_bias[nn];
  //     activation[2] = activate(activation[1], hidden_weights, bias);
  //   }
  //
  //
  //
  //
  //   // NN layers
  //   // for(nn in 1:nodes) for(n in 1:N) bias[n,nn] = input_bias[nn] ;
  //
  //   //
  //   // activation = activated(X, input_weights, input_bias, hidden_weights, hidden_bias, LL, nodes, N, P);
  //
  //   for(n in 1:N) out_bi_vec[n] = output_bias;
  //
  //   eta = activation[2] * output_weights + out_bi_vec;
  //   // eta = activation * output_weights + out_bi_vec;
  //   // if(LL == 1){
  //   //   eta =  activated2(X, input_weights, input_bias, hidden_weights, hidden_bias,
  //   //                   output_weights, output_bias, LL, nodes, N, P);
  //   // } else if (LL == 2) {
  //   //   eta =  activated3(X, input_weights, input_bias, hidden_weights, hidden_bias,
  //   //                   output_weights, output_bias, LL, nodes, N, P);
  //   // }
  //   // print(input_bias);
  //   // for(ll in 1:LL) print(hidden_weights[ll]);
  //   // for(ll in 1:LL) print(hidden_bias[ll]);
  //   // print(output_weights);
  //   // print(output_bias);
  // }
}

model {
  // // prior weights
  //   output_weights_tilde ~ std_normal();
  //   // output_weights_tilde ~ std_cauchy();
  //
  //   output_bias_tilde ~ std_normal();
  //    // output_bias_tilde ~ std_cauchy();
  //
  //   // output_weights ~ cauchy(0.0,1.0);
  //   // output_bias ~ cauchy(0.0,1.0);
  //   to_vector(input_weights_tilde) ~ std_normal();
  //   input_bias_tilde ~ std_normal();
  //
  //   if(LL > 0){
  //     for(l in 1:LL){
  //       to_vector(hidden_weights_tilde[l]) ~ std_normal();
  //       hidden_bias_tilde[l] ~ std_normal();
  //       // hidden_bias_tilde[l] ~ cauchy(0.0,1.0);
  //       hb_lam[l] ~ inv_gamma(0.5, 1.0);
  //       hidden_bias_sigma[l] ~ inv_gamma(0.5, 1.0 ./hb_lam[l]);
  //
  //       // to_vector(hidden_weights_sigma[l]) ~ cauchy(0.0,1.0);
  //       to_vector(hw_lam[l]) ~ inv_gamma(0.5,1.0);
  //       to_vector(hidden_weights_sigma[l])~ inv_gamma(0.5,1.0 ./ to_vector(hw_lam[l]));
  //     }
  //   }
  //
  //   // output_weights_sigma ~ cauchy(0.0,1.0);
  //   to_vector(input_weights_sigma) ~ cauchy(0.0,1.0);
  //   ow_lam ~ inv_gamma(0.5,1.0/25.0);
  //   output_weights_sigma ~ inv_gamma(0.5,1.0 ./ ow_lam);
  //
  //   ob_lam ~ inv_gamma(0.5,1.0);
  //   output_bias_sigma ~ inv_gamma(0.5, 1 ./ ob_lam);
  // //
  //   to_vector(iw_lam) ~ inv_gamma(0.5,1.0);
  //   to_vector(input_weights_sigma) ~ inv_gamma(0.5, 1.0 ./ to_vector(iw_lam));
  //
  //   tau_tilde ~ cauchy(0.0,1.0);
  //   tau_lam ~ inv_gamma(0.5,1.0);
  //   tau_tilde ~ inv_gamma(0.5, 1/tau_lam);
  //   // c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);
  //   c2_tilde ~ inv_gamma(half_slab_df,half_slab_df);
  //
  // //likelhihood
  //   Y ~ bernoulli_logit(eta);
}

generated quantities {
  // vector[N] prob = inv_logit(eta);
  // int Y_tilde[N] = bernoulli_rng(prob);
}
