//
// Created by Julian Stamp on 3/22/24.
//

#ifndef FAMER_SET_BLOCK_PARAMETERS_H
#define FAMER_SET_BLOCK_PARAMETERS_H

#include "fame.h"
#include "genotype.h"
#include <vector>

void set_block_parameters(genotype &genotype_block, const int &n_samples,
                          const int &block_size);

#endif // FAMER_SET_BLOCK_PARAMETERS_H
