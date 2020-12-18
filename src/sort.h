#ifndef SORT_H
#define SORT_H

#include "SLIMpaper_types.h"
#include <vector>


using namespace Rcpp;



std::vector<size_t> sort_indexes( matrix::ConstColXpr & v);

std::vector<size_t> sort_indexes( const vector & v);

void sort_indexes( matrix::ConstColXpr & v, std::vector<size_t> & idx) ;

void sort_indexes( matrix::ColXpr & v, std::vector<size_t> & idx) ;

void sort_indexes( vector & v, std::vector<size_t> & idx) ;

void sort_indexes_Eigenmat(const refMatConst & v, matrixI & idx);

void sort_matrix(refMat v) ;

void rel_sort(std::vector<size_t> idx, vector & y) ;

// template <typename Derived>
void rel_sorted_1(const Eigen::Ref< const vectorI>&  idx,
                  vector & y, const Eigen::Ref< const vector>& yorig);

#endif //SORT_H
