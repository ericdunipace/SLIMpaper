#ifndef PACKAGENAME_TYPEDEFS_H
#define PACKAGENAME_TYPEDEFS_H
#include <RcppEigen.h>

// ## // [[Rcpp::depends(RcppEigen)]]


typedef Eigen::VectorXd vector;
typedef Eigen::VectorXi vectorI;

typedef Eigen::MatrixXd matrix;
typedef Eigen::MatrixXi matrixI;
typedef Eigen::Ref<matrix> refMat;
typedef Eigen::Ref<const matrix> refMatConst;
typedef matrix::ColXpr ColXpr;
typedef matrixI::ColXpr ColXprI;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> rowMat;
typedef Eigen::LLT<matrix> llt;
typedef Eigen::LDLT<matrix> ldlt;

typedef Eigen::Map<matrix> matMap;
typedef Eigen::Map<const matrix> matMapConst;
typedef Eigen::Map<rowMat> rowMatMap;

typedef Eigen::Map<Eigen::VectorXd> vecMap;
typedef Eigen::Map<const vector> vecMapConst;
typedef Eigen::Map<Eigen::VectorXi> vecMapI;

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseVector<double> SpVec;

#endif
