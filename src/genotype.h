#pragma once
#ifndef GENOTYPE_H
#define GENOTYPE_H
#include <bits/stdc++.h>
#include "fame.h"
#include <istream>

class genotype {

public:
  std::vector<int> columnsum;
  std::vector<double> columnmeans;
  MatrixXdr allelecount_means;
  MatrixXdr allelecount_stds;

  int block_size, n_snps, n_samples, n_segments_hori, segment_size_hori;
  std::vector<std::vector<int>> p;

  double get_col_mean(int snpindex) const;
  double get_col_std(int snpindex) const;
  void clear_block();
  void set_block_parameters(const int &n, const int &b);
  void compute_block_stats();
};

#endif
