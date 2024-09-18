//
// Created by Julian Stamp on 3/28/24.
//

#ifndef MMER_FIT_COVARIATES_H
#define MMER_FIT_COVARIATES_H

#include "mme.h"

MatrixXdr &fit_covariates(const MatrixXdr &trait_mask, MatrixXdr &trait,
                          int n_samples, double y_sum, double y_mean,
                          MatrixXdr &covariate, int cov_num);

#endif // MMER_FIT_COVARIATES_H
