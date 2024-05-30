// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <random>
#include <algorithm>

#include "boost/random.hpp"
#include "genotype.h"
#include "count_data.h"
#include "read_genotypes.h"
#ifndef FAMER_SIMULATE_A_TRAIT_H
#define FAMER_SIMULATE_A_TRAIT_H

#endif //FAMER_SIMULATE_A_TRAIT_H

MatrixXdr &scale_component(float target_variance, MatrixXdr &component);

Rcpp::List simulate_a_trait_cpp(std::string plink_file,
                               float additive_heritability,
                               float gxg_heritability,
                               std::vector<int> additive_snps,
                               std::vector<int> gxg_group_1,
                               std::vector<int> gxg_group_2);