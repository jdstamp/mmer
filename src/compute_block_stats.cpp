//
// Created by Julian Stamp on 3/29/24.
//

#include "compute_block_stats.h"

void compute_block_stats(const genotype &genotype_block,
                         MatrixXdr &allelecount_means,
                         MatrixXdr &allelecount_stds, int n_samples,
                         int block_size) {
  allelecount_stds.resize(block_size, 1);
  allelecount_means.resize(block_size, 1);
  for (int i = 0; i < block_size; i++) {
    allelecount_means(i, 0) = (double)genotype_block.columnsum[i] / n_samples;
    // TODO: handle the zero variance case
    allelecount_stds(i, 0) =
        1 /
        sqrt((allelecount_means(i, 0) * (1 - (0.5 * allelecount_means(i, 0)))));
  }
}