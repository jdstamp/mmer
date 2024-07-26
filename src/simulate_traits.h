// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include "boost/random.hpp"
#include "count_data.h"
#include "genotype.h"
#include "read_genotypes.h"
#ifndef FAMER_SIMULATE_TRAITS_H
#define FAMER_SIMULATE_TRAITS_H

#endif // FAMER_SIMULATE_TRAITS_H

std::vector<int> draw_random_ints(std::vector<int> numbers, int x);

MatrixXdr draw_normal_effects(int n);

std::vector<int> readH5File(const std::string &filename,
                            const std::string &datasetName);

void replaceH5Dataset(const std::string &filename,
                      const std::string &datasetName,
                      const std::vector<int> &newData);

Rcpp::List simulate_traits_cpp(std::string plink_file,
                               float additive_heritability,
                               float gxg_heritability, int n_additive_snps,
                               std::vector<int> gxg_group_1,
                               std::vector<int> gxg_group_2);