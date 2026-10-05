functions{
  real std_normal_lpdf(vector Y){
    return -0.5 * dot_self(Y);
  }
}
