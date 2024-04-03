//
// Created by Julian Stamp on 4/2/24.
//

#ifndef FAMER_COMPUTE_MOM_COMPONENTS_H
#define FAMER_COMPUTE_MOM_COMPONENTS_H
#include "fame.h"
#include <vector>

void compute_mom_components(int n_randvecs, int n_variance_components,
                            const MatrixXdr &pheno,
                            const MatrixXdr &random_vectors, MatrixXdr &XXz,
                            MatrixXdr &Gz, const MatrixXdr &yXXy,
                            const std::vector<int> &n_snps_variance_component,
                            int n_samples_mask, MatrixXdr &S, MatrixXdr &q);

#endif // FAMER_COMPUTE_MOM_COMPONENTS_H
