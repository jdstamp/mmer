//
// Created by Julian Stamp on 3/29/24.
//

#ifndef FAMER_ALLOCATE_MEMORY_H
#define FAMER_ALLOCATE_MEMORY_H

#include "fame.h"
#include "genotype.h"

void allocate_memory(int n_randvecs, const genotype &genotype_block,
                     double *&partialsums, double *&sum_op, double *&yint_e,
                     double *&yint_m, double **&y_e, double **&y_m);

void deallocate_memory(double *partialsums, double *sum_op, double *yint_e,
                       double *yint_m, double **y_e, double **y_m,
                       const genotype &genotype_block);

#endif // FAMER_ALLOCATE_MEMORY_H
