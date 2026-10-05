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
        if((obs_t[i] - t[j]) > 0.0) count += 1;
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

  for(j in 1:NT) log_interval_length[j] = log(t[j+1] - t[j]);

  for(i in 1:N) {
    for(j in 1:NT) {
      //censor or event in interval
      count_proc = int_step(t[j + 1] - obs_t[i] + eps);

      //risk set. at risk if obs_t > t
      Y[i, j] = int_step(obs_t[i] - t[j]);

      //counting process jump = 1 if obs_t in ( t[j], t[j+1] ]
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

}
parameters {
  vector[P] beta_raw;
  real intercept_raw;
  vector[NT] log_dL0_raw;
  real<lower=0.> log_dL0_sd;
}

transformed parameters {
  vector[NT] log_dL0 = cumulative_sum(log_dL0_raw * log_dL0_sd);// = log_dL0_raw * log_dL0_sd; //time specific hazards
  vector[N] eta = Q * beta_raw * 5.0; //linear predictor
  real intercept = intercept_raw * 5.0; // overal mean of hazard

  // log_dL0[1] = log_dL0_raw[1] * log_dL0_sd;
  // for (nt in 2:NT) {
  //   log_dL0[nt] = log_dL0[nt-1] + log_dL0_raw[nt+1] * log_dL0_sd ;
  // }
  // {
  //   real m_ldl0 =  mean(log_dL0);
  //   for(nt in 1:NT) log_dL0[nt] -= m_ldl0;
  // }
}

model {
  // create intensity (rate) parameter
  vector[count_one] log_intensity = log_int_dur; // add duration to intensity

  for(i in 1:count_one) {
      log_intensity[i] += eta[idx_obs[i]] + log_dL0[idx_time[i]] + intercept;
  }

  //priors
    intercept_raw ~ std_normal();
    beta_raw ~ std_normal();
    log_dL0_raw ~ std_normal();
    log_dL0_sd ~ std_normal();

  //likelihood
    dNvec ~ poisson_log(log_intensity);
}
generated quantities {
  vector[P] beta = R_inv * beta_raw * 5.0;
  vector[NT] baseline_S;
  matrix[NT, N] individ_S;
  int dN_out[count_one] = dNvec;
  vector[count_one] int_dur = exp(log_int_dur);
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
