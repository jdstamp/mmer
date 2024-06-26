#pragma once
#include "fame.h"

MatrixXdr compute_XXz(const int &num_snp, const MatrixXdr &Z_b,
                      const MatrixXdr &means, const MatrixXdr &stds,
                      const MatrixXdr &phenotype_mask,
                      const MatrixXdr &genotype_mask, const int &Nz,
                      const int &Nindv, genotype &genotype_block, int p,
                      const int &sel_snp_local_index, bool exclude_sel_snp);

double compute_yXXy(const int &num_snp, const MatrixXdr &y_vec,
                    const MatrixXdr &means, const MatrixXdr &stds,
                    const int &sel_snp_local_index, genotype &genotype_block,
                    const MatrixXdr &genotype_mask, const int &p,
                    const bool &exclude_sel_snp);

void multiply_y_pre_fast(const MatrixXdr &op, const int &Ncol_op,
                         MatrixXdr &res, const bool &subtract_means,
                         double *&sum_op, const genotype &g, double *&yint_m,
                         double **&y_m, const int &p, double *&partialsums);

void multiply_y_post_fast(MatrixXdr &op_orig, int Nrows_op, MatrixXdr &res,
                          bool subtract_means, int p, genotype g,
                          double *&yint_e, double **&y_e, int n);
