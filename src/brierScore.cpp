#include "CoarsePosteriorSummary_types.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix brierScore_(vector & y,
                      vectorI & event,
                      vector & times,
                      Rcpp::List & pred,
                      Rcpp::List & cens_weight,
                      matrix & cens_prob,
                      int S,
                      int S_c)
{
  int NT = times.size();
  int N = y.size();
  double Nd = double(N);
  double S_cd = double(S_c);

  matrix bs(NT, S);
  vector int_bs(S);
  // double brier = 0;
  double temp_brier = 0;

  for (int t = 0; t < NT; t++) {
    matrix p = Rcpp::as<matrix>(pred[t]);
    matrix cens_mat = Rcpp::as<matrix>(cens_weight[t]); //cens_matrix by time
    for (int s = 0; s < S; s++) {
      double brier = 0;
      for ( int n = 0; n < N; n++) { //calculate for each data point
        double temp_brier = 0.0;

        for( int s_c = 0; s_c < S_c; s_c++) { //integrate over independent censoring probs
          double gs = cens_mat(s_c, n);
          double gi = cens_prob(s_c, n); //cens prob when event happened
          if ( y(n) <= times(t) ) {
            if(event(n) == 1 ) temp_brier += p(n, s) * p(n, s) / gi / S_cd;
          } else {
            temp_brier += ( 1. - p(n, s) ) * ( 1. - p(n, s) ) / gs / S_cd;
          }
        }
        brier += temp_brier / Nd;
      }
      bs(t, s) = brier;
    }

  }

  return Rcpp::wrap(bs);
}

/* survival probabilities */
// void pecSRC(double *pec,
//             double *Y,
//             double *D,
//             double *times,
//             double *pred,
//             double *weight,
//             double *weight_obs,
//             int *N,
//             int *NT,
//             int *cmodel,
//             int *ConstantPrediction)
// {
//   int s, i;
//   double p, brier, gs, gi;
//
//   for (s=0; s<(*NT);s++) {
//     for (i=0; i<*N;i++){
//       /* prediction */
//       if (*ConstantPrediction==0){
//         p = pred[i + s * (*N)];
//       }
//       else{
//         p = pred[s];
//       }
//       /* weights */
//       gs = weight[(i + s * (*N)) * (*cmodel) + s * (1-(*cmodel))];
//       gi = weight_obs[i];
//       if (Y[i] <= times[s])
//         brier = D[i] * p * p / gi;
//       else
//         brier = (1-p)*(1-p) / gs;
//       pec[s] += brier / (double) (*N);
//     }
//   }
// }
