#pragma once
#include "fame.h"

MatrixXdr compute_XXz(const int &num_snp, const MatrixXdr &Z_b,
                      const MatrixXdr &means, const MatrixXdr &stds,
                      const MatrixXdr &mask, const int &Nz, const int &Nindv,
                      const int &sel_snp_local_index, double *&sum_op,
                      const genotype &g, double *&yint_m, double **&y_m, int p,
                      double *&yint_e, double **&y_e, double *&partialsums,
                      bool exclude_sel_snp);

MatrixXdr compute_XXUz(const int &num_snp, const int &Nz, const int &Nindv,
                       const MatrixXdr &means, const MatrixXdr &stds,
                       const MatrixXdr &mask, double *&sum_op,
                       const genotype &g, double *&yint_m, double **&y_m,
                       const int &p, double *&yint_e, double **&y_e,
                       double *&partialsums);

MatrixXdr compute_Xz(int num_snp, int Nz, int Nindv, MatrixXdr means,
                     MatrixXdr mask, int p, genotype g, double &yint_e,
                     double &y_e);

void multiply_y_pre_fast(const MatrixXdr &op, const int &Ncol_op,
                         MatrixXdr &res, const bool &subtract_means,
                         double *&sum_op, const genotype &g, double *&yint_m,
                         double **&y_m, const int &p, double *&partialsums);

void multiply_y_post_fast(MatrixXdr &op_orig, int Nrows_op, MatrixXdr &res,
                          bool subtract_means, int p, genotype g,
                          double *&yint_e, double **&y_e, int n);

double compute_yXXy(int num_snp, MatrixXdr y_vec, MatrixXdr Z_b,
                    MatrixXdr means, MatrixXdr stds, int Nz, double *&sum_op,
                    genotype g, double *&yint_m, double **&y_m, int p,
                    int sel_snp_local_index, double *&partialsums);

double compute_yVXXVy(const int &num_snp, const MatrixXdr &new_pheno,
                      const MatrixXdr &means, const MatrixXdr &stds,
                      const int &Nz, double *&sum_op, const genotype &g,
                      double *&yint_m, double **&y_m, const int &p,
                      double *&partialsums);

MatrixXdr compute_XXy(const int &num_snp, const MatrixXdr &y_vec,
                      const MatrixXdr &means, const MatrixXdr &stds,
                      const MatrixXdr &mask, const int &sel_snp_local_index,
                      const int &Nindv, double *&sum_op, const genotype &g,
                      double *&yint_m, double **&y_m, const int &p,
                      double *&yint_e, double **&y_e, double *&partialsums,
                      const bool &exclude_sel_snp);

double compute_yXXy(const int &num_snp, const MatrixXdr &y_vec,
                    const MatrixXdr &means, const MatrixXdr &stds,
                    const int &sel_snp_local_index, double *&sum_op,
                    const genotype &g, double *&yint_m, double **&y_m,
                    const int &p, double *&partialsums,
                    const bool &exclude_sel_snp);
