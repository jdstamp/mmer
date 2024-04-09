// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "boost/random.hpp"
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

Rcpp::List fame_cpp(std::string plink_file, std::string pheno_file,
                    std::string covariate_file, int n_randvecs, int n_blocks,
                    int rand_seed, std::vector<int> gxg_indices,
                    std::string genotype_mask_file);
