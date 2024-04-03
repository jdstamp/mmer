//
// Created by Julian Stamp on 3/29/24.
//
#include "fame.h"

#ifndef FAMER_COMPUTE_BLOCK_STATS_H
#define FAMER_COMPUTE_BLOCK_STATS_H

#include "genotype.h"

void compute_block_stats(const genotype &genotype_block,
                         MatrixXdr &allelecount_means,
                         MatrixXdr &allelecount_stds, int n_samples,
                         int block_size);

#endif // FAMER_COMPUTE_BLOCK_STATS_H
