#ifndef GENOTYPE_H
#define GENOTYPE_H
#include "fame.h"
#include <bits/stdc++.h>

class genotype {

public:
  std::vector<int> columnsum;
  std::vector<int> columnsum2;
  std::vector<double> columnmeans;
  std::vector<double> columnmeans2;

  int block_size;

  int Nsnp, Nindv, Nsegments_hori, segment_size_hori, segment_size_ver;
  std::vector<std::vector<int>> p;

  std::vector<std::vector<int>> not_O_j;
  std::vector<std::vector<int>> not_O_i;

  double get_col_mean(int snpindex) const;
  double get_col_std(int snpindex) const;
};

#endif
