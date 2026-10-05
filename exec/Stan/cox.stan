/*
 * Leuk: Cox regression
 * URL of OpenBugs' implementation:
 *   http://www.openbugs.net/Examples/Leuk.html
 * adapted from stan file:
 *   https://github.com/stan-dev/example-models/blob/master/bugs_examples/vol1/leuk/leuk.stan
 */
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
        if((obs_t[i] - t[j] + eps) > 0.0) count += 1;
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
  vector<lower=0>[NT] times; //times to calculate hazard function
  int<lower=0> fail[N]; //indicator for failure ==1 or censoring == 0
  int<lower=0> P; //number of predictors
  matrix[N,P] X; //matrix of predictors
}
transformed data {
  int Y[N, NT];
  int dN[N, NT];
  real eps = .000000001;
  vector<lower=0.>[NT + 1] t = create_times(times);
  int count_one = count_steps(obs_t, t, eps);
  int idx_time[count_one];
  int idx_obs[count_one];
  int dNvec[count_one];
  int count = 0;
  vector[NT] gamma_vec;
  real r = 0.01;
  real c = 0.001;

  // int sums = 0;
  // t[1] = 0.0;
  // for(nt in 1:NT) {
  //   t[nt+1] = times[nt];
  // }

  for(i in 1:N) {
    for(j in 1:NT) {
      //risk set. at risk if obs_t > t
      Y[i, j] = int_step(obs_t[i] - t[j] + eps);
      //counting process jump = 1 if obs_t in [ t[j], t[j+1] )
      dN[i, j] = Y[i, j] * fail[i] * int_step(t[j + 1] - obs_t[i] - eps);
    }
  }

  // for(i in 1:N) for(j in 1:NT) if(Y[i,j] == 1) sums +=1;
  // print(sums);
  // print(count_one);
  for(i in 1:N) {
    for(j in 1:NT) {
      if(Y[i,j] == 1) {
        count += 1;
        dNvec[count] = dN[i,j];
        idx_time[count] = j;
        idx_obs[count] = i;
      }
    }
  }

  for(j in 1:NT) gamma_vec[j] = r * (t[j + 1] - t[j]) * c;

}
parameters {
  vector[P] beta;
  real<lower=0> dL0[NT];
}

transformed parameters {
  vector[N] eta;
  vector[N] rate;
  // vector[N] log_rate;

  eta = X * beta;
  rate = exp(eta);
  // for(n in 1:N) rate[n] = exp(eta[n]);
}
model {
  real log_p = 0.0;
  //priors
  beta ~ std_normal();
  dL0 ~ gamma(gamma_vec, c);

  //likelihood
  for(j in 1:NT) {
    for(i in 1:N) {
      if (Y[i, j] != 0)
        log_p += poisson_lpmf(dN[i, j] | Y[i, j] * rate[i] * dL0[j]);
    }
  }
  target += log_p;
}
generated quantities {
  vector[NT] baseline_S;
  matrix[NT, N] individ_S;
  vector[NT] log_dL0;
  real intercept;
  int dN_out[N, NT] = dN;

  for (j in 1:NT) {
    log_dL0[j] = log(dL0[j]);
  }
  intercept = mean(log_dL0);
  for (j in 1:NT) {
    // Survivor function = exp(-Integral{l0(u)du})^exp(beta*z)
    real s;
    s = 0;
    for (i in 1:j)
      s = s + dL0[i];
    baseline_S[j] = exp(-s);
    for(n in 1:N) individ_S[j,n] = pow(baseline_S[j], rate[n]);
  }

}
