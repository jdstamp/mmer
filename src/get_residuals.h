#pragma once
#include "mme.h"
#include "count_data.h"
#include "fit_covariates.h"
#include "read_covariates.h"
#include "read_phenotypes.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <string>

Eigen::MatrixXd get_residuals(std::string pheno_file, std::string covariate_file);

