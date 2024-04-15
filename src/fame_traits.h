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
#ifndef FAMER_SIMULATE_TRAITS_H
#define FAMER_SIMULATE_TRAITS_H

#endif //FAMER_SIMULATE_TRAITS_H

Rcpp::List
fame_traits_cpp(std::string plink_file, float additive_variance_component, float gxg_variance_component,
int focal_snp, std::vector<int> gxg_indices,
std::vector<int> additive_indices);