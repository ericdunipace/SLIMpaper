functions{
  real std_cauchy_lpdf(vector Y) {
    return - sum(log1p(square(Y)));
  }
  real std_cauchy_real_lpdf(real Y) {
    return (- log1p(square(Y)));
  }
  int count_steps(vector obs_t, vector t, real eps){
    int N = rows(obs_t);
    int NT = rows(t)-1;
    int count = 0;

    for(i in 1:N) {
      for(j in 1:NT) {
        if((obs_t[i] - t[j] ) > 0.0) count += 1;
      }
    }
    return (count);
  }
  vector create_times(vector times) {
    int NT = rows(times);
    vector[ NT + 1] t;
    t[1] = 0.0;
    for(i in 1:NT) t[i + 1] = times[i];
    return(t);
  }
}

data {
  int<lower=1> N; //number of data points
  int<lower=1> NT; //number of hazard times
  vector<lower=0>[N] obs_t; //observed times for data
  vector<lower=0>[NT+1] times; //times to calculate hazard function
  int<lower=0> fail[N]; //indicator for failure ==1 or censoring == 0
  int<lower=0> P; //number of predictors
  matrix[N,P] X; //matrix of predictors
  real m0; //number of non-zero parameters hypothesized
}
transformed data {
  int Y[N, NT];
  int dN[N, NT];
  real eps = .000000001;
  real inv_sqrt_N = 1/sqrt(1.0 * N);
  real m_scale = m0 / (P * 1.0 - m0);
  real slab_scale = 1.0;    // Scale for large slopes
  real slab_scale2 = square(slab_scale);
  real slab_df = 25.;      // Effective degrees of freedom for large slopes
  real half_slab_df = 0.5 * slab_df;
  vector<lower=0.>[NT + 1] t = times;
  matrix[N,NT] interval;
  vector[NT] log_interval_length;
  int count_one = count_steps(obs_t, t, eps);
  int idx_time[count_one];
  int idx_obs[count_one];
  int dNvec[count_one];
  vector[count_one] log_int_dur;
  int count = 0;
  int count_proc = 0;
  matrix[N,P] Q = qr_Q(X)[, 1:P] * sqrt(N - 1.0);
  matrix[P,P] R = qr_R(X)[1:P, ] / sqrt(N - 1.0);
  matrix[P,P] R_inv = inverse(R);

  // print(count_one);
  for(j in 1:NT) log_interval_length[j] = log(t[j+1] - t[j]);

  for(i in 1:N) {
    for(j in 1:NT) {
      //censor or event in interval
      count_proc = int_step(t[j + 1] - obs_t[i] + eps);

      //risk set. at risk if obs_t > t
      Y[i, j] = int_step(obs_t[i] - t[j]);

      //counting process jump = 1 if obs_t in ( t[j], t[j+1] ] and is a failure
      dN[i, j] = Y[i, j] * fail[i] * count_proc;

      // calculate duration in interval
      if(count_proc == 1 && Y[i,j] == 1) {
        interval[i,j] = obs_t[i] - t[j];
      } else {
        interval[i,j] = t[j+1] - t[j];
      }
    }
  }

  for(i in 1:N) {
    for(j in 1:NT) {
      if(Y[i,j] == 1) {
        count += 1;
        dNvec[count] = dN[i,j];
        idx_time[count] = j;
        idx_obs[count] = i;
        log_int_dur[count] = log(interval[i,j]);
      }
    }
  }
  // print(count);

}
parameters {
  vector[P] beta_raw;
  vector<lower = 0.>[P] beta_sd;
  real intercept_raw;
  vector[NT] log_dL0_raw;
  real<lower=0.> log_dL0_sd;
  real<lower=0.> tau_tilde;
  real<lower=0.> c2_tilde;
}

transformed parameters {
  vector[NT] log_dL0 = cumulative_sum(log_dL0_raw * log_dL0_sd);
  vector[P] beta_tilde;
  vector[N] eta;
  real intercept = intercept_raw * 5;

  {
    real tau_0 = m_scale * inv_sqrt_N;
    real tau = tau_tilde * tau_0;
    real c2 = slab_scale2 * c2_tilde;
    vector[P] lambda = sqrt(c2) * beta_sd ./ sqrt(c2 + square(tau) * square(beta_sd));
    beta_tilde = beta_raw .* lambda * tau;
  }
  eta = Q * beta_tilde;

  // log_dL0[1] = log_dL0_raw[1] + log_dL0_raw[2] * log_dL0_sd;
  // for (nt in 2:NT) {
  //   log_dL0[nt] = log_dL0[nt-1] + log_dL0_raw[nt+1] * log_dL0_sd;
  // }
  // {
  // real m_ldl0 =  mean(log_dL0);
  // for(nt in 1:NT) log_dL0[nt] -= m_ldl0;
  // }
}

model {
  vector[count_one] log_intensity = log_int_dur; // add duration to intensity

  for(i in 1:count_one) {
      log_intensity[i] += eta[idx_obs[i]] + log_dL0[idx_time[i]] + intercept;
      // log_intensity[i] += eta[idx_obs[i]] + intercept;
  }

  //priors
    intercept_raw ~ std_normal();
    tau_tilde ~ std_cauchy_real();
    c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);
    beta_sd ~ std_cauchy();
    beta_raw ~ std_normal();
    log_dL0_raw ~ std_normal();
    log_dL0_sd ~ std_normal();

  //likelihood
    dNvec ~ poisson_log(log_intensity);
}
generated quantities {
  vector[P] beta = R_inv * beta_tilde;
  vector[NT] baseline_S;
  matrix[NT, N] individ_S;
  // int dN_out[count_one] = dNvec;
  // vector[count_one] int_dur = exp(log_int_dur);
  vector[NT] log_hazard = intercept + log_dL0;

  // for (j in 1:NT) {
  //   log_hazard[j] = 0.; //delete
  // }

  for (j in 1:NT) {
    // Survivor function = exp(-Integral{l0(u)du})^exp(beta*z)
    real s = 0;
    for (i in 1:j)
      s = s + exp(log_hazard[i] + log_interval_length[i]);
      // s = s + exp(intercept + log_int_dur[i]);
    baseline_S[j] = exp(-s);
    for(n in 1:N) individ_S[j,n] = pow(baseline_S[j], exp(eta[n]));
  }

}

