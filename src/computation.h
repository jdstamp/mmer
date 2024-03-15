#pragma once
#include "fame.h"

MatrixXdr compute_XXz(int num_snp, MatrixXdr Z_b, MatrixXdr means,
                      MatrixXdr stds, MatrixXdr mask, int Nz, int Nindv,
                      int sel_snp_local_index, double *&sum_op, genotype g,
                      double *&yint_m, double **&y_m, int p, double *&yint_e,
                      double **&y_e, double *&partialsums);

MatrixXdr compute_XXUz(int num_snp, int Nz, int Nindv, MatrixXdr means,
                       MatrixXdr stds, MatrixXdr mask, double *&sum_op,
                       genotype g, double *&yint_m, double **&y_m, int p,
                       double *&yint_e, double **&y_e, double *&partialsums);

MatrixXdr compute_Xz(int num_snp, int Nz, int Nindv, MatrixXdr means,
                     MatrixXdr mask, int p, genotype g, double &yint_e,
                     double &y_e);

void multiply_y_pre_fast(MatrixXdr &op, int Ncol_op, MatrixXdr &res,
                         bool subtract_means, double *&sum_op, genotype g,
                         double *&yint_m, double **&y_m, int p, double *&partialsums);

void multiply_y_post_fast(MatrixXdr &op_orig, int Nrows_op, MatrixXdr &res,
                          bool subtract_means, int p, genotype g,
                          double *&yint_e, double **&y_e, int n);

double compute_yXXy(int num_snp, MatrixXdr y_vec, MatrixXdr Z_b,
                    MatrixXdr means, MatrixXdr stds, int Nz, double *&sum_op,
                    genotype g, double *&yint_m, double **&y_m, int p,
                    int sel_snp_local_index, double *&partialsums);

double compute_yVXXVy(int num_snp, MatrixXdr new_pheno, MatrixXdr means,
                      MatrixXdr stds, int Nz, double *&sum_op, genotype g,
                      double *&yint_m, double **&y_m, int p, double *&partialsums);

MatrixXdr compute_XXy(int num_snp, MatrixXdr y_vec, MatrixXdr means,
                      MatrixXdr stds, MatrixXdr mask, int sel_snp_local_index,
                      int Nindv, double *&sum_op, genotype g, double *&yint_m,
                      double **&y_m, int p, double *&yint_e, double **&y_e, double *&partialsums);

double compute_yXXy(int num_snp, MatrixXdr y_vec, MatrixXdr means,
                    MatrixXdr stds, int sel_snp_local_index, double *&sum_op,
                    genotype g, double *&yint_m, double **&y_m, int p, double *&partialsums);
