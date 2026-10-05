functions{
  real std_cauchy_lpdf(vector Y){
    return - sum(log1p(square(Y)));
  }
}
