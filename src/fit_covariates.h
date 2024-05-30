//
// Created by Julian Stamp on 3/28/24.
//

#ifndef FAMER_FIT_COVARIATES_H
#define FAMER_FIT_COVARIATES_H

#include "fame.h"

MatrixXdr &fit_covariates(const MatrixXdr &trait_mask, MatrixXdr &trait,
                          int n_samples, double y_sum, double y_mean,
                          MatrixXdr &covariate, int cov_num);

#endif // FAMER_FIT_COVARIATES_H
