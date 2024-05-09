/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include "allocate_memory.h"
#include "computation.h"
#include "compute_block_stats.h"
#include "compute_mom_components.h"
#include "fame.h"
#include "genotype.h"
#include "initialize_random_vectors.h"
#include "read_genotypes.h"
#include "read_phenotypes.h"
#include "set_block_parameters.h"
#include "set_metadata.h"
#include <testthat.h>

#include <fstream>
#include <iostream>
#include <unistd.h>

// this should come from one configuration file
std::string testdata_dir_mqs = "../../inst/testdata/";
std::string test_bed_mqs = testdata_dir_mqs + "test.bed";
std::string test_csv_mqs = testdata_dir_mqs + "test.csv";
std::string test_pheno_mqs = testdata_dir_mqs + "test_h2_0.5.pheno";

double tolerance_mqs = 1e-6;
int n_samples_mqs = 200;
int block_size_mqs = 10;
metaData metadata_mqs = set_metadata(n_samples_mqs, block_size_mqs);

MatrixXdr readCSVToMatrixXdr_mqs(const std::string &filename) {
  std::ifstream data(filename);
  std::string line;
  std::vector<double> values;
  int rows = 0;
  int cols = 0;

  // Skip the header line
  std::getline(data, line);

  // Read data, line by line
  while (std::getline(data, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    int temp_cols = 0;

    // Read each cell
    while (std::getline(lineStream, cell, ',')) {
      values.push_back(std::stod(cell));
      temp_cols++;
    }

    // Update the number of columns
    if (cols == 0) {
      cols = temp_cols;
    }

    rows++;
  }

  // Convert the vector of values into a MatrixXdr
  MatrixXdr mat(rows, cols);
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      mat(i, j) = values[i * cols + j];
  return mat;
}

void centerAndStandardize(MatrixXdr &mat) {
  // Compute column-wise mean
  MatrixXdr mean = mat.colwise().mean();

  // iterate over elements in mean and subtract it from corresponding mat
  // columns
  for (int i = 0; i < mean.cols(); i++) {
    mat.col(i) = mat.col(i).array() - mean(0, i);
  }

  // Compute column-wise standard deviation
  MatrixXdr std_dev =
      ((mat.array().square().colwise().sum()) / (mat.rows() - 1)).sqrt();

  for (int i = 0; i < mean.cols(); i++) {
    mat.col(i) = mat.col(i).array() / std_dev(0, i);
  }
}

void centerAndStandardize2(MatrixXdr &mat) {
  // Compute column-wise mean
  MatrixXdr mean = mat.colwise().mean();

  // iterate over elements in mean and subtract it from corresponding mat
  // columns
  for (int i = 0; i < mean.cols(); i++) {
    mat.col(i) = mat.col(i).array() - mean(0, i);
  }

  MatrixXdr stds = MatrixXdr::Zero(mean.cols(), 1);
  for (int i = 0; i < mean.cols(); i++) {
    stds(i, 0) = 1 / sqrt((mean(0, i) * (1 - (0.5 * mean(0, i)))));
  }

  for (int i = 0; i < mean.cols(); i++) {
    mat.col(i) = mat.col(i).array() * stds(i, 0);
  }
}

context("C++ test mqs method is implemented correctly") {
  test_that("test compute_mom_components") {
    // given
    MatrixXdr data = readCSVToMatrixXdr_mqs(test_csv_mqs);
    MatrixXdr matrix_block = data.block(0, 0, n_samples_mqs, block_size_mqs);
    // MatrixXdr TWO = MatrixXdr::Ones(n_samples_mqs, block_size_mqs) * 2;
    // matrix_block = TWO - matrix_block; // because encoding is inverted
    centerAndStandardize2(
        matrix_block); // TODO: check if we can switch to proper
                       // standardizing (version 1 above)

    MatrixXdr pheno_mask;
    MatrixXdr pheno;
    read_phenotypes(n_samples_mqs, test_pheno_mqs, pheno, pheno_mask);

    // select focal SNP as the first column of the matrix block
    MatrixXdr focal_snp = matrix_block.col(0);

    // compute K, G
    MatrixXdr K = matrix_block * matrix_block.transpose() / block_size_mqs;
    MatrixXdr G = (K * block_size_mqs - focal_snp * focal_snp.transpose()) /
                  (block_size_mqs - 1);

    for (int i = 0; i < G.cols(); ++i) {
      G.col(i) = G.col(i).array() * focal_snp.array();
      G.row(i) = G.row(i).array() * focal_snp.transpose().array();
    }

    MatrixXdr q_expected(3, 1);
    q_expected(0, 0) = ((K * pheno).array() * pheno.array()).sum();
    q_expected(1, 0) = ((G * pheno).array() * pheno.array()).sum();
    q_expected(2, 0) = (pheno.array() * pheno.array()).sum();

    MatrixXdr S_expected(3, 3);
    S_expected(0, 0) = (K * K).trace();
    S_expected(1, 0) = (K * G).trace();
    S_expected(0, 1) = (G * K).trace();
    S_expected(1, 1) = (G * G).trace();
    S_expected(2, 0) = K.trace();
    S_expected(0, 2) = K.trace();
    S_expected(1, 2) = G.trace();
    S_expected(2, 1) = G.trace();
    S_expected(2, 2) = n_samples_mqs;

    MatrixXdr sigma_expected = S_expected.inverse() * q_expected;

    // when
    int n_randvecs = 1000;
    int rand_seed = 1234;

    int n_variance_components = 2;
    int focal_snp_local_index = 0;
    genotype genotype_block;

    MatrixXdr genotype_mask_matrix = MatrixXdr::Ones(block_size_mqs, 1);

    set_block_parameters(genotype_block, n_samples_mqs, block_size_mqs);
    std::ifstream bed_ifs(test_bed_mqs.c_str(), ios::in | ios::binary);
    int global_snp_index = -1;
    read_genotype_block(bed_ifs, block_size_mqs, genotype_block, n_samples_mqs,
                        global_snp_index, metadata_mqs);
    MatrixXdr fame_means(block_size_mqs, 1);
    MatrixXdr fame_stds(block_size_mqs, 1);
    compute_block_stats(genotype_block, fame_means, fame_stds, n_samples_mqs,
                        block_size_mqs);

    // memory setup for mailman
    double *partialsums = new double[0];
    double *sum_op;
    double *yint_e;
    double *yint_m;
    double **y_e;
    double **y_m;
    allocate_memory(n_randvecs, genotype_block, partialsums, sum_op, yint_e,
                    yint_m, y_e, y_m);

    bool exclude_sel_snp = false;

    MatrixXdr yXXy;
    yXXy = MatrixXdr::Zero(n_variance_components, 1);
    yXXy(0, 0) += compute_yXXy(block_size_mqs, pheno, fame_means, fame_stds,
                               focal_snp_local_index, sum_op, genotype_block,
                               genotype_mask_matrix, yint_m, y_m,
                               block_size_mqs, partialsums, exclude_sel_snp);

    MatrixXdr gxg_pheno;
    gxg_pheno = pheno.array() * focal_snp.col(0).array();
    yXXy(1, 0) += compute_yXXy(block_size_mqs, gxg_pheno, fame_means, fame_stds,
                               focal_snp_local_index, sum_op, genotype_block,
                               genotype_mask_matrix, yint_m, y_m,
                               block_size_mqs, partialsums, true);

    MatrixXdr random_vectors;
    MatrixXdr gxg_random_vectors;

    random_vectors = initialize_random_vectors(
        n_randvecs, rand_seed, pheno_mask, random_vectors, n_samples_mqs);
    gxg_random_vectors =
        random_vectors.array().colwise() * focal_snp.col(0).array();

    MatrixXdr GxGz;
    MatrixXdr XXz;
    MatrixXdr collect_XXy;
    MatrixXdr collect_XXUy;
    MatrixXdr temp_grm;
    MatrixXdr temp_gxg;

    XXz = MatrixXdr::Zero(n_samples_mqs, n_randvecs);
    GxGz = MatrixXdr::Zero(n_samples_mqs, n_randvecs);
    collect_XXy = MatrixXdr::Zero(n_samples_mqs, n_variance_components + 1);
    collect_XXUy =
        MatrixXdr::Zero(n_samples_mqs, (n_variance_components + 1) *
                                           (n_variance_components + 1));

    temp_grm = compute_XXz(
        block_size_mqs, random_vectors, fame_means, fame_stds, pheno_mask,
        genotype_mask_matrix, n_randvecs, n_samples_mqs, sum_op, genotype_block,
        yint_m, y_m, block_size_mqs, yint_e, y_e, partialsums, 0, false);

    for (int z_index = 0; z_index < n_randvecs; z_index++) {
      XXz.col(z_index) += temp_grm.col(z_index);
    }
    bool in_gxg_block = true;
    temp_gxg =
        compute_XXz(block_size_mqs, gxg_random_vectors, fame_means, fame_stds,
                    pheno_mask, genotype_mask_matrix, n_randvecs, n_samples_mqs,
                    sum_op, genotype_block, yint_m, y_m, block_size_mqs, yint_e,
                    y_e, partialsums, focal_snp_local_index, in_gxg_block);
    temp_gxg = temp_gxg.array().colwise() * focal_snp.col(0).array();

    for (int z_index = 0; z_index < n_randvecs; z_index++) {
      GxGz.col(z_index) += temp_gxg.col(z_index);
    }
    collect_XXy.col(0) +=
        compute_XXz(block_size_mqs, pheno, fame_means, fame_stds, pheno_mask,
                    genotype_mask_matrix, 1, n_samples_mqs, sum_op,
                    genotype_block, yint_m, y_m, block_size_mqs, yint_e, y_e,
                    partialsums, focal_snp_local_index, false);

    MatrixXdr temp_Gy =
        compute_XXz(block_size_mqs, gxg_pheno, fame_means, fame_stds,
                    pheno_mask, genotype_mask_matrix, 1, n_samples_mqs, sum_op,
                    genotype_block, yint_m, y_m, block_size_mqs, yint_e, y_e,
                    partialsums, focal_snp_local_index, true);
    temp_Gy = temp_Gy.array() * focal_snp.col(0).array();
    collect_XXy.col(1) += temp_Gy;

    collect_XXy.col(0) = collect_XXy.col(0) / block_size_mqs;
    collect_XXy.col(1) = collect_XXy.col(1) / (block_size_mqs - 1);
    collect_XXy.col(2) = pheno;
    //

    MatrixXdr scaled_vec;
    for (int i = 0; i < (n_variance_components + 1); i++) {
      MatrixXdr temp_XXUy =
          compute_XXz(block_size_mqs, collect_XXy.col(i), fame_means, fame_stds,
                      pheno_mask, genotype_mask_matrix, 1, n_samples_mqs,
                      sum_op, genotype_block, yint_m, y_m, block_size_mqs,
                      yint_e, y_e, partialsums, focal_snp_local_index, false);

      scaled_vec = collect_XXy.col(i).array() * focal_snp.col(0).array();
      MatrixXdr temp_GUy =
          compute_XXz(block_size_mqs, scaled_vec, fame_means, fame_stds,
                      pheno_mask, genotype_mask_matrix, 1, n_samples_mqs,
                      sum_op, genotype_block, yint_m, y_m, block_size_mqs,
                      yint_e, y_e, partialsums, focal_snp_local_index, false);
      temp_GUy = temp_GUy.array() * focal_snp.col(0).array();

      collect_XXUy.col(i) += temp_XXUy / block_size_mqs;
      collect_XXUy.col(((n_variance_components + 1)) + i) +=
          temp_GUy / (block_size_mqs - 1);
    }

    for (int i = 0; i < (n_variance_components + 1); i++) {
      collect_XXUy.col((n_variance_components * (n_variance_components + 1)) +
                       i) = collect_XXy.col(i);
    }
    vector<int> n_snps_variance_component = {block_size_mqs,
                                             block_size_mqs - 1};
    int n_samples_mask = pheno_mask.sum();
    MatrixXdr S_observed(n_variance_components + 1, n_variance_components + 1);
    MatrixXdr q_observed(n_variance_components + 1, 1);
    compute_mom_components(n_randvecs, n_variance_components, pheno,
                           random_vectors, XXz, GxGz, yXXy,
                           n_snps_variance_component, n_samples_mask,
                           S_observed, q_observed);

    MatrixXdr sigma_observed = S_observed.inverse() * q_observed;

    //    // then
    float S_error =
        std::sqrt((S_observed.array() - S_expected.array()).square().sum());
    float q_error = (q_observed.array() - q_expected.array()).square().sum();

    float sigma_error =
        (sigma_observed.array() - sigma_expected.array()).square().sum();

    expect_true(S_error < 100);
    expect_true(q_error < tolerance_mqs);
    expect_true(sigma_error < tolerance_mqs);
  }
}
