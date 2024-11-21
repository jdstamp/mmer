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

#include <testthat.h>

#include "testing_utils.h"

namespace local {
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
} // namespace local

context("C++ test covariance is computed correctly") {
  test_that("test second bed file") {
    correctTestFiles(test_csv, test_bed_2, test_pheno, test_h5, test_ld_h5);
    correctTestFiles(test_csv, test_bed, test_pheno, test_h5, test_ld_h5);
    // given
    int n_grm_snps = 200;
    int n_samples = 200;
    int n_gxg_snps = 100;
    int focal_snp_local_index = 96;

    MatrixXdr pheno_mask;
    MatrixXdr pheno;
    read_phenotypes(n_samples, test_pheno, pheno, pheno_mask);

    MatrixXdr grm_matrix_block = MatrixXdr::Zero(n_samples, n_grm_snps);
    MatrixXdr gxg_matrix_block = MatrixXdr::Zero(n_samples, n_gxg_snps);
    MatrixXdr snp_matrix = MatrixXdr::Zero(n_samples, 1);
    std::ifstream bed_ifs(test_bed_2.c_str(), ios::in | ios::binary);
    int global_snp_index = -1;
    for (int i = 0; i < n_grm_snps; i++) {
      read_snp(bed_ifs, global_snp_index, snp_matrix);
      grm_matrix_block.col(i) = snp_matrix;
    }

    std::ifstream gxg_bed_ifs(test_bed.c_str(), ios::in | ios::binary);
    global_snp_index = -1;
    for (int i = 0; i < n_gxg_snps; i++) {
      read_snp(gxg_bed_ifs, global_snp_index, snp_matrix);
      gxg_matrix_block.col(i) = snp_matrix;
    }

    local::centerAndStandardize2(grm_matrix_block);
    local::centerAndStandardize2(gxg_matrix_block);

    // compute K from grm_matrix_block
    MatrixXdr K = grm_matrix_block * grm_matrix_block.transpose() / n_grm_snps;
    // compute G from gxg_matrix_block
    MatrixXdr focal_snp = gxg_matrix_block.col(focal_snp_local_index);
    MatrixXdr G = gxg_matrix_block * gxg_matrix_block.transpose();
    G = (G - focal_snp * focal_snp.transpose()) / (n_gxg_snps - 1);

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
    S_expected(2, 2) = n_samples;

    // compute yKy and yGy
    MatrixXdr Sinv = S_expected.inverse();
    MatrixXdr sigma_expected = Sinv * q_expected;

    MatrixXdr V =
        sigma_expected(0, 0) * K + sigma_expected(1, 0) * G +
        sigma_expected(2, 0) * MatrixXdr::Identity(n_samples, n_samples);
    MatrixXdr H = Sinv(0, 1) * K + Sinv(1, 1) * G +
                  Sinv(2, 1) * MatrixXdr::Identity(n_samples, n_samples);

    MatrixXdr sigma_variance = 2 * pheno.transpose() * H * V * H * pheno;

    // when
    int n_randvecs = 100;
    int rand_seed = 1234;

    int n_variance_components = 2;
    genotype grm_genotype_block;
    genotype gxg_genotype_block;

    grm_genotype_block.set_block_parameters(n_samples, n_grm_snps);
    gxg_genotype_block.set_block_parameters(n_samples, n_gxg_snps - 1);

    bed_ifs.seekg(0, std::ios::beg);     // reset file pointer to beginning
    gxg_bed_ifs.seekg(0, std::ios::beg); // reset file pointer to beginning
    global_snp_index = -1;

    for (int i = 0; i < n_grm_snps; i++) {
      read_snp(bed_ifs, global_snp_index, snp_matrix);
      grm_genotype_block.encode_snp(snp_matrix);
    }
    global_snp_index = -1;
    for (int i = 0; i < n_gxg_snps; i++) {
      read_snp(gxg_bed_ifs, global_snp_index, snp_matrix);
      if (i != focal_snp_local_index)
        gxg_genotype_block.encode_snp(snp_matrix);
    }
    //
    grm_genotype_block.compute_block_stats();
    gxg_genotype_block.compute_block_stats();

    MatrixXdr yXXy;
    yXXy = MatrixXdr::Zero(n_variance_components, 1);
    yXXy(0, 0) += compute_yXXy(grm_genotype_block, pheno);

    MatrixXdr gxg_pheno;
    gxg_pheno = pheno.array() * focal_snp.col(0).array();
    yXXy(1, 0) += compute_yXXy(gxg_genotype_block, gxg_pheno);

    MatrixXdr random_vectors;
    MatrixXdr gxg_random_vectors;

    random_vectors = initialize_random_vectors(
        n_randvecs, rand_seed, pheno_mask, random_vectors, n_samples);
    gxg_random_vectors =
        random_vectors.array().colwise() * focal_snp.col(0).array();

    MatrixXdr GxGz;
    MatrixXdr XXz;
    MatrixXdr collect_XXy;
    MatrixXdr collect_XXUy;
    MatrixXdr temp_grm;
    MatrixXdr temp_gxg;

    XXz = MatrixXdr::Zero(n_samples, n_randvecs);
    GxGz = MatrixXdr::Zero(n_samples, n_randvecs);
    collect_XXy = MatrixXdr::Zero(n_samples, n_variance_components + 1);
    collect_XXUy = MatrixXdr::Zero(n_samples, (n_variance_components + 1) *
                                                  (n_variance_components + 1));

    temp_grm =
        compute_XXz(random_vectors, pheno_mask, n_randvecs, grm_genotype_block);

    for (int z_index = 0; z_index < n_randvecs; z_index++) {
      XXz.col(z_index) += temp_grm.col(z_index);
    }
    temp_gxg = compute_XXz(gxg_random_vectors, pheno_mask, n_randvecs,
                           gxg_genotype_block);
    temp_gxg = temp_gxg.array().colwise() * focal_snp.col(0).array();

    for (int z_index = 0; z_index < n_randvecs; z_index++) {
      GxGz.col(z_index) += temp_gxg.col(z_index);
    }
    collect_XXy.col(0) += compute_XXz(pheno, pheno_mask, 1, grm_genotype_block);

    MatrixXdr temp_Gy =
        compute_XXz(gxg_pheno, pheno_mask, 1, gxg_genotype_block);
    temp_Gy = temp_Gy.array() * focal_snp.col(0).array();
    collect_XXy.col(1) += temp_Gy;

    collect_XXy.col(0) = collect_XXy.col(0) / n_grm_snps;
    collect_XXy.col(1) = collect_XXy.col(1) / (n_gxg_snps - 1);
    collect_XXy.col(2) = pheno;

    MatrixXdr scaled_vec;
    for (int i = 0; i < (n_variance_components + 1); i++) {
      MatrixXdr temp_XXUy =
          compute_XXz(collect_XXy.col(i), pheno_mask, 1, grm_genotype_block);

      scaled_vec = collect_XXy.col(i).array() * focal_snp.col(0).array();
      MatrixXdr temp_GUy =
          compute_XXz(scaled_vec, pheno_mask, 1, gxg_genotype_block);
      temp_GUy = temp_GUy.array() * focal_snp.col(0).array();

      collect_XXUy.col(i) += temp_XXUy / block_size;
      collect_XXUy.col(((n_variance_components + 1)) + i) +=
          temp_GUy / (n_gxg_snps - 1);
    }

    for (int i = 0; i < (n_variance_components + 1); i++) {
      collect_XXUy.col((n_variance_components * (n_variance_components + 1)) +
                       i) = collect_XXy.col(i);
    }
    vector<int> n_snps_variance_component = {n_grm_snps, n_gxg_snps - 1};
    int n_samples_mask = pheno_mask.sum();
    MatrixXdr S_observed(n_variance_components + 1, n_variance_components + 1);
    MatrixXdr q_observed(n_variance_components + 1, 1);
    compute_mom_components(n_randvecs, n_variance_components, pheno,
                           random_vectors, XXz, GxGz, yXXy(0, 0), yXXy(1, 0),
                           n_snps_variance_component, n_samples_mask,
                           S_observed, q_observed);

    MatrixXdr sigma_observed = S_observed.inverse() * q_observed;

    MatrixXdr cov_q_observed;
    compute_covariance_q(n_variance_components, collect_XXUy, sigma_observed,
                         cov_q_observed);

    MatrixXdr invS = S_observed.inverse();

    MatrixXdr cov_sigma = invS * cov_q_observed * invS;

    //    // then
    float S_error =
        std::sqrt((S_observed.array() - S_expected.array()).square().sum());
    float q_error = (q_observed.array() - q_expected.array()).square().sum();
    float sigma_error =
        (sigma_observed.array() - sigma_expected.array()).square().sum();

    expect_true(S_error < 1e2);
    expect_true(q_error < 1e-5);
    expect_true(sigma_error < 1e-4);
  }
}
