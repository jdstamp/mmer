#include "boost/random.hpp"
#include <Rcpp.h>
#include <RcppEigen.h>

#include "genotype.h"

#if SSE_SUPPORT == 1
#define fastmultiply fastmultiply_sse
#define fastmultiply_pre fastmultiply_pre_sse
#else
#define fastmultiply fastmultiply_normal
#define fastmultiply_pre fastmultiply_pre_normal
#endif

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    MatrixXdr;
