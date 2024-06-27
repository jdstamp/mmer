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

#include "computation.h"
#include "allocate_memory.h"
#include "mailman.h"

MatrixXdr compute_XXz(const MatrixXdr &Z_b, const MatrixXdr &means,
                      const MatrixXdr &stds, const MatrixXdr &phenotype_mask,
                      const int &Nz, const int &Nindv, genotype &genotype_block,
                      const int &sel_snp_local_index, bool exclude_sel_snp) {

  double *partialsums2 = new double[0];
  double *sum_op2;
  double *yint_e2;
  double *yint_m2;
  double **y_e2;
  double **y_m2;

  allocate_memory(Nz, genotype_block, partialsums2, sum_op2, yint_e2, yint_m2,
                  y_e2, y_m2);

  MatrixXdr res;
  res.resize(genotype_block.block_size, Nz);
  multiply_y_pre_fast(Z_b, Nz, res, false, sum_op2, genotype_block, yint_m2,
                      y_m2, genotype_block.block_size, partialsums2);
  MatrixXdr zb_sum = Z_b.colwise().sum();

  for (int j = 0; j < genotype_block.block_size; j++)
    for (int k = 0; k < Nz; k++) {
      res(j, k) = res(j, k) * stds(j, 0);
    }

  MatrixXdr resid(genotype_block.block_size, Nz);
  MatrixXdr inter = means.cwiseProduct(stds);
  resid = inter * zb_sum;
  MatrixXdr inter_zb = res - resid;

  // GxG case
  if (exclude_sel_snp == true)
    for (int k = 0; k < Nz; k++)
      inter_zb(sel_snp_local_index, k) = 0;

  // masking
  for (int j = 0; j < genotype_block.block_size; j++)
    for (int k = 0; k < Nz; k++) {
      inter_zb(j, k) = inter_zb(j, k) * stds(j, 0);
    }

  MatrixXdr new_zb = inter_zb.transpose();
  MatrixXdr new_res(Nz, Nindv);

  multiply_y_post_fast(new_zb, Nz, new_res, false, genotype_block.block_size,
                       genotype_block, yint_e2, y_e2, Nindv);

  MatrixXdr new_resid(1, Nindv);
  MatrixXdr zb_scale_sum = new_zb * means;

  new_resid = zb_scale_sum * MatrixXdr::Constant(1, Nindv, 1);

  /// new zb
  MatrixXdr temp = new_res - new_resid;

  for (int i = 0; i < Nz; i++)
    for (int j = 0; j < Nindv; j++)
      temp(i, j) = temp(i, j) * phenotype_mask(j, 0);

  return temp.transpose();
}

double compute_yXXy(const MatrixXdr &y_vec, const MatrixXdr &means,
                    const MatrixXdr &stds, const int &sel_snp_local_index,
                    genotype &genotype_block, const bool &exclude_sel_snp) {
  double *partialsums2 = new double[0];
  double *sum_op2;
  double *yint_e2;
  double *yint_m2;
  double **y_e2;
  double **y_m2;
  int Nz = 1;

  allocate_memory(Nz, genotype_block, partialsums2, sum_op2, yint_e2, yint_m2,
                  y_e2, y_m2);

  MatrixXdr res;
  res.resize(genotype_block.block_size, 1);

  multiply_y_pre_fast(y_vec, 1, res, false, sum_op2, genotype_block, yint_m2,
                      y_m2, genotype_block.block_size, partialsums2);

  res = res.cwiseProduct(stds);
  MatrixXdr resid(genotype_block.block_size, 1);
  resid = means.cwiseProduct(stds);
  resid = resid * y_vec.sum();

  MatrixXdr Xy(genotype_block.block_size, 1);
  Xy = res - resid;

  // GxG case
  if (exclude_sel_snp == true)
    Xy(sel_snp_local_index, 0) = 0;

  double yXXy = (Xy.array() * Xy.array()).sum();

  return yXXy;
}

void multiply_y_pre_fast(const MatrixXdr &op, const int &Ncol_op,
                         MatrixXdr &res, const bool &subtract_means,
                         double *&sum_op, const genotype &g, double *&yint_m,
                         double **&y_m, const int &p, double *&partialsums) {
  bool var_normalize = false;

  for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
    sum_op[k_iter] = op.col(k_iter).sum();
  }

  for (int seg_iter = 0; seg_iter < g.n_segments_hori - 1; seg_iter++) {
    mailman::fastmultiply(g.segment_size_hori, g.n_samples, Ncol_op,
                          g.p[seg_iter], op, yint_m, partialsums, y_m);
    int p_base = seg_iter * g.segment_size_hori;
    for (int p_iter = p_base;
         (p_iter < p_base + g.segment_size_hori) && (p_iter < g.n_snps);
         p_iter++) {
      for (int k_iter = 0; k_iter < Ncol_op; k_iter++)
        res(p_iter, k_iter) = y_m[p_iter - p_base][k_iter];
    }
  }

  int last_seg_size = (g.n_snps % g.segment_size_hori != 0)
                          ? g.n_snps % g.segment_size_hori
                          : g.segment_size_hori;
  mailman::fastmultiply(last_seg_size, g.n_samples, Ncol_op,
                        g.p[g.n_segments_hori - 1], op, yint_m, partialsums,
                        y_m);
  int p_base = (g.n_segments_hori - 1) * g.segment_size_hori;
  for (int p_iter = p_base;
       (p_iter < p_base + g.segment_size_hori) && (p_iter < g.n_snps);
       p_iter++) {
    for (int k_iter = 0; k_iter < Ncol_op; k_iter++)
      res(p_iter, k_iter) = y_m[p_iter - p_base][k_iter];
  }

  if (!subtract_means) {
    return;
  }

  for (int p_iter = 0; p_iter < p; p_iter++) {
    for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
      res(p_iter, k_iter) =
          res(p_iter, k_iter) - (g.get_col_mean(p_iter) * sum_op[k_iter]);
      if (var_normalize)
        res(p_iter, k_iter) = res(p_iter, k_iter) / (g.get_col_std(p_iter));
    }
  }
}

void multiply_y_post_fast(MatrixXdr &op_orig, int Nrows_op, MatrixXdr &res,
                          bool subtract_means, int p, genotype g,
                          double *&yint_e, double **&y_e, int n) {
  bool var_normalize = false;
  double *partialsums = new double[0];

  MatrixXdr op;
  op = op_orig.transpose();

  if (var_normalize && subtract_means) {
    for (int p_iter = 0; p_iter < p; p_iter++) {
      for (int k_iter = 0; k_iter < Nrows_op; k_iter++)
        op(p_iter, k_iter) = op(p_iter, k_iter) / (g.get_col_std(p_iter));
    }
  }

  int Ncol_op = Nrows_op;

  int seg_iter;
  for (seg_iter = 0; seg_iter < g.n_segments_hori - 1; seg_iter++) {
    mailman::fastmultiply_pre(g.segment_size_hori, g.n_samples, Ncol_op,
                              seg_iter * g.segment_size_hori, g.p[seg_iter], op,
                              yint_e, partialsums, y_e);
  }
  int last_seg_size = (g.n_snps % g.segment_size_hori != 0)
                          ? g.n_snps % g.segment_size_hori
                          : g.segment_size_hori;
  mailman::fastmultiply_pre(last_seg_size, g.n_samples, Ncol_op,
                            seg_iter * g.segment_size_hori, g.p[seg_iter], op,
                            yint_e, partialsums, y_e);

  for (int n_iter = 0; n_iter < n; n_iter++) {
    for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
      res(k_iter, n_iter) = y_e[n_iter][k_iter];
      y_e[n_iter][k_iter] = 0;
    }
  }

  if (!subtract_means)
    return;

  double *sums_elements = new double[Ncol_op];
  memset(sums_elements, 0, Nrows_op * sizeof(int));

  for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
    double sum_to_calc = 0.0;
    for (int p_iter = 0; p_iter < p; p_iter++)
      sum_to_calc += g.get_col_mean(p_iter) * op(p_iter, k_iter);
    sums_elements[k_iter] = sum_to_calc;
  }
  for (int k_iter = 0; k_iter < Ncol_op; k_iter++) {
    for (int n_iter = 0; n_iter < n; n_iter++)
      res(k_iter, n_iter) = res(k_iter, n_iter) - sums_elements[k_iter];
  }
}
