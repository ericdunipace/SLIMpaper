#include "sort.h"

bool compare (double a, double b) {
  return a < b;
}

std::vector<size_t> sort_indexes( matrix::ConstColXpr & v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v(i1) < v(i2);});
  return idx;
}

std::vector<size_t> sort_indexes( const vector & v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v(i1) < v(i2);});
  return idx;
}

void sort_indexes( matrix::ConstColXpr & v, std::vector<size_t> & idx) {

  // initialize original index locations
  std::iota(idx.begin(), idx.end(), 0); //fills with increasing values

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v(i1) < v(i2);});
}

void sort_indexes( matrix::ColXpr & v, std::vector<size_t> & idx) {

  // initialize original index locations
  std::iota(idx.begin(), idx.end(), 0); //fills with increasing values

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v(i1) < v(i2);});
}

void sort_indexes( vector & v, std::vector<size_t> & idx) {

  // initialize original index locations
  std::iota(idx.begin(), idx.end(), 0); //fills with increasing values

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v(i1) < v(i2);});
}

void sort_indexes_Eigenmat(const refMatConst & v, matrixI & idx) {

  int N = v.cols();
  int P = v.rows();

  for(int n = 0; n < N;  n++) {
    // initialize original index locations
    idx.col(n) = Eigen::VectorXi::LinSpaced(Eigen::Sequential,P,0,P-1); //fills with increasing values

    // sort indexes based on comparing values in v
    std::sort(idx.col(n).data(), idx.col(n).data() + P,
              [&v,n](size_t i1, size_t i2) {return v(i1,n) < v(i2,n);});
    // for(int i =0; i < idx.size(); i++) idx(i,n) = idx_temp[i];
  }

}

void sort_matrix(refMat v) {

  int N = v.cols();
  int P = v.rows();

  // for(auto n : v.colwise()) {
  //   std::sort(n.begin(), n.end())
  // } //available in future eigen.

  for(int n = 0; n < N;  n++) {
    // sort indexes based on comparing values in v
    std::sort(v.col(n).data(), v.col(n).data() + P);
  }

}

void rel_sort(std::vector<size_t> idx, vector & y) {
  vector temp_sort = y;
  std::sort(temp_sort.data(), temp_sort.data() + temp_sort.size());
  for(int i = 0; i < y.size(); i ++) y(idx[i]) = temp_sort(i);
}

// template <typename Derived>
void rel_sorted_1(const Eigen::Ref< const vectorI>&  idx,
                  vector & y, const Eigen::Ref< const vector>& yorig) {
  for(int i = 0; i < y.size(); i ++) y(idx(i)) = yorig(i);
}
