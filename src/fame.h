// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <iostream>
#include <stdexcept>

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
                    std::string covariate_file, int n_randvecs,
                    int focal_snp_index, int n_blocks, int rand_seed);


void
allocate_memory(int n_randvecs, const genotype &genotype_block, int &hsegsize,
                double *&partialsums, double *&sum_op, double *&yint_e,
                double *&yint_m, double **&y_e, double **&y_m);

void deallocate_memory(int hsegsize, double *partialsums, double *sum_op,
                       double *yint_e, double *yint_m, double **y_e,
                       double **y_m, genotype &genotype_block);