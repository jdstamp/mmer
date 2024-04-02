/*
 * Substantial parts of this file were published under MIT license by other
 * authors. Find the original license notice below.
 */

/* MIT License
 *
 * Copyright (c) 2024 sriramlab
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "fit_covariates.h"
#include "fame.h"
#include <Rcpp.h>
#include <RcppEigen.h>

MatrixXdr &fit_covariates(const MatrixXdr &trait_mask, MatrixXdr &trait,
                          int n_samples, double y_sum, double y_mean,
                          MatrixXdr &covariate, int cov_num) {
  MatrixXdr v1; // W^ty
  MatrixXdr v2; // QW^ty
  MatrixXdr v3; // WQW^ty
  MatrixXdr Q;
  MatrixXdr new_pheno;
  MatrixXdr mat_mask = trait_mask.replicate(1, cov_num);
  covariate = covariate.cwiseProduct(mat_mask);

  MatrixXdr WtW = covariate.transpose() * covariate;
  Q = WtW.inverse(); // Q=(W^tW)^-1

  double phen_sd = 0;
  v1 = covariate.transpose() * trait; // W^ty
  v2 = Q * v1;                        // QW^ty
  v3 = covariate * v2;                // WQW^ty
  new_pheno = trait - v3;
  trait = new_pheno.cwiseProduct(trait_mask);

  y_sum = trait.sum();
  y_mean = y_sum / trait_mask.sum();
  for (int i = 0; i < n_samples; i++) {
    phen_sd += (trait(i, 0) - y_mean) * (trait(i, 0) - y_mean);
    if (trait(i, 0) != 0)
      trait(i, 0) = trait(i, 0) - y_mean; // center phenotype
  }
  phen_sd = sqrt(phen_sd / (trait_mask.sum() - 1));
  trait = trait / phen_sd;
  return trait;
}