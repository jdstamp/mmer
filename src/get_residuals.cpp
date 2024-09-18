#include "get_residuals.h"

// [[Rcpp::export]]
Eigen::MatrixXd get_residuals(std::string pheno_file, std::string covariate_file) {
  MatrixXdr pheno_mask;
  MatrixXdr pheno;
  int n_samples = count_samples(pheno_file);
  read_phenotypes(n_samples, pheno_file, pheno, pheno_mask);
  MatrixXdr covariate;
  int cov_num;
  bool standardize = false;
  cov_num = read_covariates(standardize, n_samples, covariate_file, covariate);
  double y_sum = 0;
  double y_mean = 0;
  MatrixXdr residuals(pheno.rows(), pheno.cols());
    // Iterate over the columns of pheno
    for (int col = 0; col < pheno.cols(); ++col) {
        double y_sum = 0;
        double y_mean = 0;
        MatrixXdr pheno_col = pheno.col(col);
        MatrixXdr residual_col = fit_covariates(pheno_mask, pheno_col, n_samples, y_sum, y_mean, covariate, cov_num);
        residuals.col(col) = residual_col;
    }

    return residuals;
}