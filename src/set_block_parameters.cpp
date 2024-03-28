//
// Created by Julian Stamp on 3/22/24.
//
#include "set_block_parameters.h"

void set_block_parameters(genotype &genotype_block, const int &n_samples,
                          const int &block_size) {
  genotype_block.segment_size_hori =
      floor(log(n_samples) / log(3)) - 2; // object of the mailman
  genotype_block.n_segments_hori =
      ceil(block_size * 1.0 / (genotype_block.segment_size_hori * 1.0));
  genotype_block.p.resize(genotype_block.n_segments_hori,
                          std::vector<int>(n_samples));
  genotype_block.not_O_i.resize(block_size);
  genotype_block.not_O_j.resize(n_samples);
  genotype_block.block_size = 0;
  genotype_block.n_snps = block_size;
  genotype_block.n_samples = n_samples;

  genotype_block.columnsum.resize(block_size, 1);
  for (int index_temp = 0; index_temp < block_size; index_temp++) {
    genotype_block.columnsum[index_temp] = 0;
  }
}