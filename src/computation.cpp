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
#include "mailman.h"

MatrixXdr compute_XXz(const int &num_snp, const MatrixXdr &Z_b,
                      const MatrixXdr &means, const MatrixXdr &stds,
                      const MatrixXdr &phenotype_mask,
                      const MatrixXdr &genotype_mask, const int &Nz,
                      const int &Nindv, const int &sel_snp_local_index,
                      double *&sum_op, const genotype &genotype_block,
                      double *&yint_m, double **&y_m, int p, double *&yint_e,
                      double **&y_e, double *&partialsums,
                      bool exclude_sel_snp) {
  MatrixXdr res;
  res.resize(num_snp, Nz);
  multiply_y_pre_fast(Z_b, Nz, res, false, sum_op, genotype_block, yint_m, y_m,
                      p, partialsums);

  MatrixXdr zb_sum = Z_b.colwise().sum();

  for (int j = 0; j < num_snp; j++)
    for (int k = 0; k < Nz; k++) {
      res(j, k) = res(j, k) * stds(j, 0);
      res(j, k) = res(j, k) * genotype_mask(j, 0);
    }

  ////GxG
  if (exclude_sel_snp == true)
    for (int k = 0; k < Nz; k++)
      res(sel_snp_local_index, k) = 0;

  MatrixXdr resid(num_snp, Nz);
  MatrixXdr inter = means.cwiseProduct(stds);
  resid = inter * zb_sum;
  MatrixXdr inter_zb = res - resid;

  for (int k = 0; k < Nz; k++)
    for (int j = 0; j < num_snp; j++)
      inter_zb(j, k) = inter_zb(j, k) * stds(j, 0);
  MatrixXdr new_zb = inter_zb.transpose();
  MatrixXdr new_res(Nz, Nindv);

  multiply_y_post_fast(new_zb, Nz, new_res, false, p, genotype_block, yint_e,
                       y_e, Nindv);

  MatrixXdr new_resid(Nz, num_snp);
  MatrixXdr zb_scale_sum = new_zb * means;
  new_resid = zb_scale_sum * MatrixXdr::Constant(1, Nindv, 1);

  /// new zb
  MatrixXdr temp = new_res - new_resid;

  for (int i = 0; i < Nz; i++)
    for (int j = 0; j < Nindv; j++)
      temp(i, j) = temp(i, j) * phenotype_mask(j, 0);

  return temp.transpose();
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
  double *partialsums;

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

MatrixXdr
compute_XXy(const int &num_snp, const MatrixXdr &y_vec, const MatrixXdr &means,
            const MatrixXdr &stds, const MatrixXdr &phenotype_mask,
            const MatrixXdr &genotype_mask, const int &sel_snp_local_index,
            const int &Nindv, double *&sum_op, const genotype &genotype_block,
            double *&yint_m, double **&y_m, const int &p, double *&yint_e,
            double **&y_e, double *&partialsums, const bool &exclude_sel_snp) {
  MatrixXdr res;
  res.resize(num_snp, 1);
  multiply_y_pre_fast(y_vec, 1, res, false, sum_op, genotype_block, yint_m, y_m,
                      p, partialsums);

  MatrixXdr zb_sum = y_vec.colwise().sum();

  for (int j = 0; j < num_snp; j++) {
    res(j, 0) = res(j, 0) * stds(j, 0);
    res(j, 0) = res(j, 0) * genotype_mask(j, 0);
  }

  // GxG
  if (exclude_sel_snp == true) {
    res(sel_snp_local_index, 0) = 0;
  }

  MatrixXdr resid(num_snp, 1);
  MatrixXdr inter = means.cwiseProduct(stds);
  resid = inter * zb_sum;
  MatrixXdr inter_zb = res - resid;

  for (int j = 0; j < num_snp; j++)
    inter_zb(j, 0) = inter_zb(j, 0) * stds(j, 0);

  MatrixXdr new_zb = inter_zb.transpose();
  MatrixXdr new_res(1, Nindv);
  multiply_y_post_fast(new_zb, 1, new_res, false, p, genotype_block, yint_e,
                       y_e, Nindv);
  MatrixXdr new_resid(1, num_snp);
  MatrixXdr zb_scale_sum = new_zb * means;
  new_resid = zb_scale_sum * MatrixXdr::Constant(1, Nindv, 1);

  MatrixXdr temp = new_res - new_resid;

  for (int j = 0; j < Nindv; j++)
    temp(0, j) = temp(0, j) * phenotype_mask(j, 0);

  return temp.transpose();
}

double compute_yXXy(const int &num_snp, const MatrixXdr &y_vec,
                    const MatrixXdr &means, const MatrixXdr &stds,
                    const int &sel_snp_local_index, double *&sum_op,
                    const genotype &genotype_block, double *&yint_m,
                    double **&y_m, const int &p, double *&partialsums,
                    const bool &exclude_sel_snp) {
  MatrixXdr res;          // Boyang: add declarative res
  res.resize(num_snp, 1); // Boyang: change to resize

  multiply_y_pre_fast(y_vec, 1, res, false, sum_op, genotype_block, yint_m, y_m,
                      p, partialsums);

  /// GxG
  //  bool exclude_sel_snp = false;
  if (exclude_sel_snp == true)
    res(sel_snp_local_index, 0) = 0;

  res = res.cwiseProduct(stds);
  MatrixXdr resid(num_snp, 1);
  resid = means.cwiseProduct(stds);
  resid = resid * y_vec.sum();
  MatrixXdr Xy(num_snp, 1);
  Xy = res - resid;

  double yXXy = (Xy.array() * Xy.array()).sum();

  return yXXy;
}
