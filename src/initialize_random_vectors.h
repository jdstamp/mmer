//
// Created by Julian Stamp on 3/29/24.
//
#include "fame.h"

#ifndef FAMER_INITIALIZE_RANDOM_VECTORS_H
#define FAMER_INITIALIZE_RANDOM_VECTORS_H

MatrixXdr &initialize_random_vectors(int n_randvecs, int rand_seed,
                                     const MatrixXdr &pheno_mask,
                                     MatrixXdr &random_vectors, int n_samples);

#endif // FAMER_INITIALIZE_RANDOM_VECTORS_H
