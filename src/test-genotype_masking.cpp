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

context("C++ test yXXy with masking") {
  test_that("mailman algorithm reproduces naive implementation for yXXy") {
    // given
    int n_variance_components = 1;
    int focal_snp_local_index = 0;
    genotype genotype_block;

    std::vector<int> genotype_mask = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int mask_size = genotype_mask.size();
    MatrixXdr genotype_mask_matrix = MatrixXdr::Zero(block_size, 1);
    for (int i = 0; i < block_size; i++) {
      genotype_mask_matrix(i, 0) = genotype_mask[i];
    }

    set_block_parameters(genotype_block, n_samples, block_size);
    std::ifstream bed_ifs(test_bed.c_str(), ios::in | ios::binary);
    int global_snp_index = -1;
    read_genotype_block(bed_ifs, block_size, genotype_block, n_samples,
                        global_snp_index, metadata);

    MatrixXdr data = readCSVToMatrixXdr(test_csv);

    MatrixXdr fame_means(block_size, 1);
    MatrixXdr fame_stds(block_size, 1);
    compute_block_stats(genotype_block, fame_means, fame_stds, n_samples,
                        block_size);

    MatrixXdr matrix_block = data.block(0, 0, n_samples, block_size);
    MatrixXdr TWO = MatrixXdr::Ones(n_samples, block_size) * 2;
    matrix_block = TWO - matrix_block;
    MatrixXdr means(block_size, 1);
    MatrixXdr stds(block_size, 1);
    means = matrix_block.colwise().sum().transpose() / n_samples;

    // compute stds_minor by iterating over the columns
    for (int i = 0; i < block_size; i++) {
      stds(i, 0) = 1 / sqrt((means(i, 0) * (1 - (0.5 * means(i, 0)))));
    }
    for (int i = 0; i < block_size; i++) {
      matrix_block.col(i) = matrix_block.col(i).array() - means(i, 0);
      matrix_block.col(i) = matrix_block.col(i).array() * stds(i, 0);
    }

    MatrixXdr pheno_mask;
    MatrixXdr pheno;
    read_phenotypes(n_samples, test_pheno, pheno, pheno_mask);

    MatrixXdr Xy(block_size, n_samples);
    MatrixXdr yXXy_expected(1, 1);
    Xy = matrix_block.transpose() * pheno;
    yXXy_expected(0, 0) = (Xy.array() * Xy.array()).sum();

    // when
    bool exclude_sel_snp = false;

    MatrixXdr yXXy_observed;
    yXXy_observed = MatrixXdr::Zero(n_variance_components, 1);
    yXXy_observed(0, 0) += compute_yXXy(
        block_size, pheno, fame_means, fame_stds, focal_snp_local_index,
        genotype_block, block_size, exclude_sel_snp);

    // then
    expect_true(std::abs(yXXy_expected(0, 0) - yXXy_observed(0, 0)) <
                tolerance);
  }

  test_that("masking gives the expected result for yXXy") {
    // given
    int n_variance_components = 1;
    int focal_snp_local_index = 0;
    genotype genotype_block;

    std::vector<int> genotype_mask = {0, 1, 1, 0, 1, 0, 1, 1, 0, 0};
    // convert the genotype mask to a MatrixXdr
    MatrixXdr genotype_mask_matrix = MatrixXdr::Zero(block_size, 1);
    for (int i = 0; i < block_size; i++) {
      genotype_mask_matrix(i, 0) = genotype_mask[i];
    }
    int mask_size = genotype_mask_matrix.sum();

    set_block_parameters(genotype_block, n_samples, mask_size);
    std::ifstream bed_ifs(test_bed.c_str(), ios::in | ios::binary);
    int global_snp_index = -1;

    MatrixXdr snp_matrix = MatrixXdr::Zero(n_samples, 1);
for (int i = 0; i < block_size; i++) {
read_snp(bed_ifs, n_samples, global_snp_index, metadata, snp_matrix);
if (genotype_mask_matrix(i, 0) == 1) {
encode_snp(genotype_block, snp_matrix);
}
}

    MatrixXdr data = readCSVToMatrixXdr(test_csv);
    MatrixXdr matrix_block = data.block(0, 0, n_samples, block_size);
    MatrixXdr masked_matrix_block(n_samples, mask_size);
    int count = 0;
    for (int i = 0; i < block_size; i++) {
      if (genotype_mask[i] == 1) {
        masked_matrix_block.col(count) = matrix_block.col(i);
        count++;
      }
    }

    MatrixXdr fame_means(mask_size, 1);
    MatrixXdr fame_stds(mask_size, 1);
    compute_block_stats(genotype_block, fame_means, fame_stds, n_samples,
            mask_size);

    MatrixXdr TWO = MatrixXdr::Ones(n_samples, mask_size) * 2;
    masked_matrix_block = TWO - masked_matrix_block;
    MatrixXdr means(mask_size, 1);
    MatrixXdr stds(mask_size, 1);
    means = masked_matrix_block.colwise().mean().transpose();
    for (int i = 0; i < mask_size; i++) {
      stds(i, 0) = 1 / sqrt((means(i, 0) * (1 - (0.5 * means(i, 0)))));
    }
    for (int i = 0; i < mask_size; i++) {
      masked_matrix_block.col(i) =
          masked_matrix_block.col(i).array() - means(i, 0);
      masked_matrix_block.col(i) =
          masked_matrix_block.col(i).array() * stds(i, 0);
    }

    MatrixXdr pheno_mask;
    MatrixXdr pheno;
    read_phenotypes(n_samples, test_pheno, pheno, pheno_mask);

    MatrixXdr Xy(mask_size, n_samples);
    MatrixXdr yXXy_expected(1, 1);
    Xy = masked_matrix_block.transpose() * pheno;
    yXXy_expected(0, 0) = (Xy.array() * Xy.array()).sum();

    // when
    bool exclude_sel_snp = false;

    MatrixXdr yXXy_observed;
    yXXy_observed = MatrixXdr::Zero(n_variance_components, 1);
    yXXy_observed(0, 0) += compute_yXXy(
        mask_size, pheno, fame_means, fame_stds, focal_snp_local_index,
        genotype_block, mask_size, exclude_sel_snp);

    // then
    expect_true(std::abs(yXXy_expected(0, 0) - yXXy_observed(0, 0)) <
                tolerance);
  }
}

context("C++ test major or minor allele count encoding gives same result") {
  test_that("Standard deviation and mean are equivalent for major or minor "
            "allele count encoding") {
    // given
    genotype genotype_block;

    MatrixXdr means_minor(block_size, 1);
    MatrixXdr stds_minor(block_size, 1);

    MatrixXdr means_major(block_size, 1);
    MatrixXdr stds_major(block_size, 1);

    MatrixXdr fame_means(block_size, 1);
    MatrixXdr fame_stds(block_size, 1);

    set_block_parameters(genotype_block, n_samples, block_size);
    std::ifstream bed_ifs(test_bed.c_str(), ios::in | ios::binary);
    int global_snp_index = -1;
    read_genotype_block(bed_ifs, block_size, genotype_block, n_samples,
                        global_snp_index, metadata);

    MatrixXdr data = readCSVToMatrixXdr(test_csv);
    MatrixXdr matrix_block_minor = data.block(0, 0, n_samples, block_size);
    MatrixXdr matrix_block_major(n_samples, block_size);
    MatrixXdr TWO = MatrixXdr::Ones(n_samples, block_size) * 2;

    // when
    compute_block_stats(genotype_block, fame_means, fame_stds, n_samples,
                        block_size);

    matrix_block_major = TWO - matrix_block_minor;
    means_minor = matrix_block_minor.colwise().mean().transpose();
    means_major = matrix_block_major.colwise().mean().transpose();

    for (int i = 0; i < block_size; i++) {
      stds_minor(i, 0) =
          1 / sqrt((means_minor(i, 0) * (1 - (0.5 * means_minor(i, 0)))));
      stds_major(i, 0) =
          1 / sqrt((means_major(i, 0) * (1 - (0.5 * means_major(i, 0)))));
    }

    double err_stds_naive = (stds_minor - stds_major).array().abs().sum();
    double err_means_fame = (means_major - fame_means).array().abs().sum();
    double err_stds_fame = (stds_major - fame_stds).array().abs().sum();

    // then

    expect_true(err_stds_naive < tolerance);
    expect_true(err_means_fame < tolerance);
    expect_true(err_stds_fame < tolerance);
  }
}

context("C++ test XXz with masking") {
  test_that("mailman algorithm reproduces naive implementation for XXz") {
    // given
    int n_variance_components = 1;
    int focal_snp_local_index = 0;
    int n_randvecs = 1;
    genotype genotype_block;

    std::vector<int> genotype_mask = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int mask_size = genotype_mask.size();
    MatrixXdr genotype_mask_matrix = MatrixXdr::Zero(block_size, 1);
    for (int i = 0; i < block_size; i++) {
      genotype_mask_matrix(i, 0) = genotype_mask[i];
    }

    set_block_parameters(genotype_block, n_samples, block_size);
    std::ifstream bed_ifs(test_bed.c_str(), ios::in | ios::binary);
    int global_snp_index = -1;
    read_genotype_block(bed_ifs, block_size, genotype_block, n_samples,
                        global_snp_index, metadata);

    MatrixXdr data = readCSVToMatrixXdr(test_csv);

    MatrixXdr fame_means(block_size, 1);
    MatrixXdr fame_stds(block_size, 1);
    compute_block_stats(genotype_block, fame_means, fame_stds, n_samples,
                        block_size);

    MatrixXdr matrix_block = data.block(0, 0, n_samples, block_size);
    MatrixXdr TWO = MatrixXdr::Ones(n_samples, block_size) * 2;
    matrix_block = TWO - matrix_block;
    MatrixXdr means(block_size, 1);
    MatrixXdr stds(block_size, 1);
    means = matrix_block.colwise().sum().transpose() / n_samples;
    // compute stds_minor by iterating over the columns
    for (int i = 0; i < block_size; i++) {
      stds(i, 0) = 1 / sqrt((means(i, 0) * (1 - (0.5 * means(i, 0)))));
    }
    for (int i = 0; i < block_size; i++) {
      matrix_block.col(i) = matrix_block.col(i).array() - means(i, 0);
      matrix_block.col(i) = matrix_block.col(i).array() * stds(i, 0);
    }

    MatrixXdr pheno_mask;
    MatrixXdr pheno;
    read_phenotypes(n_samples, test_pheno, pheno, pheno_mask);

    MatrixXdr Xz(block_size, n_samples);
    MatrixXdr XXz_expected(n_samples, 1);
    Xz = matrix_block.transpose() * pheno;
    XXz_expected = matrix_block * Xz;

    // when
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

    MatrixXdr XXz_observed;
    XXz_observed = compute_XXz(block_size, pheno, fame_means, fame_stds,
                               pheno_mask, n_randvecs,
                               n_samples, genotype_block, block_size, 0, false);

    double abs_error = (XXz_expected - XXz_observed).array().abs().sum();

    // then
    expect_true(abs_error < tolerance);
  }

  test_that("masking gives the expected result for XXz") {
    // given
    int n_variance_components = 1;
    int focal_snp_local_index = 0;
    int n_randvecs = 1;
    genotype genotype_block;

    std::vector<int> genotype_mask = {0, 1, 1, 0, 1, 0, 1, 1, 0, 0};
    MatrixXdr genotype_mask_matrix = MatrixXdr::Zero(block_size, 1);
    for (int i = 0; i < block_size; i++) {
      genotype_mask_matrix(i, 0) = genotype_mask[i];
    }
    int mask_size = genotype_mask_matrix.sum();

    set_block_parameters(genotype_block, n_samples, mask_size);
    std::ifstream bed_ifs(test_bed.c_str(), ios::in | ios::binary);
    int global_snp_index = -1;

MatrixXdr snp_matrix = MatrixXdr::Zero(n_samples, 1);
for (int i = 0; i < block_size; i++) {
read_snp(bed_ifs, n_samples, global_snp_index, metadata, snp_matrix);
if (genotype_mask_matrix(i, 0) == 1) {
encode_snp(genotype_block, snp_matrix);
}
}

    MatrixXdr data = readCSVToMatrixXdr(test_csv);
    MatrixXdr matrix_block = data.block(0, 0, n_samples, block_size);
    MatrixXdr masked_matrix_block(n_samples, mask_size);
    int count = 0;
    for (int i = 0; i < block_size; i++) {
      if (genotype_mask[i] == 1) {
        masked_matrix_block.col(count) = matrix_block.col(i);
        count++;
      }
    }

    MatrixXdr fame_means(mask_size, 1);
    MatrixXdr fame_stds(mask_size, 1);
    compute_block_stats(genotype_block, fame_means, fame_stds, n_samples,
            mask_size);

    MatrixXdr TWO = MatrixXdr::Ones(n_samples, mask_size) * 2;
    masked_matrix_block = TWO - masked_matrix_block;
    MatrixXdr means(mask_size, 1);
    MatrixXdr stds(mask_size, 1);
    means = masked_matrix_block.colwise().mean().transpose();
    for (int i = 0; i < mask_size; i++) {
      stds(i, 0) = 1 / sqrt((means(i, 0) * (1 - (0.5 * means(i, 0)))));
    }
    for (int i = 0; i < mask_size; i++) {
      masked_matrix_block.col(i) =
          masked_matrix_block.col(i).array() - means(i, 0);
      masked_matrix_block.col(i) =
          masked_matrix_block.col(i).array() * stds(i, 0);
    }

    MatrixXdr pheno_mask;
    MatrixXdr pheno;
    read_phenotypes(n_samples, test_pheno, pheno, pheno_mask);

    MatrixXdr Xz(block_size, n_samples);
    MatrixXdr XXz_expected(n_samples, 1);
    Xz = masked_matrix_block.transpose() * pheno;
    XXz_expected = masked_matrix_block * Xz;

    // when

    bool exclude_sel_snp = false;

    MatrixXdr XXz_observed = MatrixXdr::Zero(n_samples, 1);
    XXz_observed = compute_XXz(mask_size, pheno, fame_means, fame_stds,
                               pheno_mask, n_randvecs,
                               n_samples, genotype_block, mask_size, 0, false);

    double abs_error = (XXz_expected - XXz_observed).array().abs().sum();

    // then
    expect_true(abs_error < tolerance);
  }
}

context("C++ test XXy with masking") {
  test_that("mailman algorithm reproduces naive implementation for XXy") {
    // given
    int n_variance_components = 1;
    int focal_snp_local_index = 0;
    int n_randvecs = 1;
    genotype genotype_block;

    std::vector<int> genotype_mask = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int mask_size = genotype_mask.size();
    MatrixXdr genotype_mask_matrix = MatrixXdr::Zero(block_size, 1);
    for (int i = 0; i < block_size; i++) {
      genotype_mask_matrix(i, 0) = genotype_mask[i];
    }

    set_block_parameters(genotype_block, n_samples, block_size);
    std::ifstream bed_ifs(test_bed.c_str(), ios::in | ios::binary);
    int global_snp_index = -1;
    read_genotype_block(bed_ifs, block_size, genotype_block, n_samples,
                        global_snp_index, metadata);

    MatrixXdr data = readCSVToMatrixXdr(test_csv);

    MatrixXdr fame_means(block_size, 1);
    MatrixXdr fame_stds(block_size, 1);
    compute_block_stats(genotype_block, fame_means, fame_stds, n_samples,
                        block_size);

    MatrixXdr matrix_block = data.block(0, 0, n_samples, block_size);
    MatrixXdr TWO = MatrixXdr::Ones(n_samples, block_size) * 2;
    matrix_block = TWO - matrix_block;
    MatrixXdr means(block_size, 1);
    MatrixXdr stds(block_size, 1);
    means = matrix_block.colwise().sum().transpose() / n_samples;
    // compute stds_minor by iterating over the columns
    for (int i = 0; i < block_size; i++) {
      stds(i, 0) = 1 / sqrt((means(i, 0) * (1 - (0.5 * means(i, 0)))));
    }
    for (int i = 0; i < block_size; i++) {
      matrix_block.col(i) = matrix_block.col(i).array() - means(i, 0);
      matrix_block.col(i) = matrix_block.col(i).array() * stds(i, 0);
    }

    MatrixXdr pheno_mask;
    MatrixXdr pheno;
    read_phenotypes(n_samples, test_pheno, pheno, pheno_mask);

    MatrixXdr Xy(block_size, n_samples);
    MatrixXdr XXy_expected(n_samples, 1);
    Xy = matrix_block.transpose() * pheno;
    XXy_expected = matrix_block * Xy;

    // when
    // memory setup for mailman
    double *partialsums = new double[0];
    double *sum_op;
    double *yint_e;
    double *yint_m;
    double **y_e;
    double **y_m;
    allocate_memory(n_randvecs, genotype_block, partialsums, sum_op, yint_e,
                    yint_m, y_e, y_m);

    MatrixXdr XXy_observed = MatrixXdr::Zero(n_samples, 1);

    XXy_observed.col(0) +=
        compute_XXz(block_size, pheno, fame_means, fame_stds, pheno_mask,
                    n_randvecs, n_samples, genotype_block,
                    block_size, focal_snp_local_index, false);

    double abs_error = (XXy_expected - XXy_observed).array().abs().sum();

    // then
    expect_true(abs_error < tolerance);
  }

  test_that("masking gives the expected result for XXy") {
    // given
    int n_variance_components = 1;
    int focal_snp_local_index = 0;
    int n_randvecs = 1;
    genotype genotype_block;

    std::vector<int> genotype_mask = {0, 1, 1, 0, 1, 0, 1, 1, 0, 0};
    MatrixXdr genotype_mask_matrix = MatrixXdr::Zero(block_size, 1);
    for (int i = 0; i < block_size; i++) {
      genotype_mask_matrix(i, 0) = genotype_mask[i];
    }
    int mask_size = genotype_mask_matrix.sum();

    set_block_parameters(genotype_block, n_samples, mask_size);
    std::ifstream bed_ifs(test_bed.c_str(), ios::in | ios::binary);
    int global_snp_index = -1;
//    read_genotype_block(bed_ifs, block_size, genotype_block, n_samples,
//                        global_snp_index, metadata);

MatrixXdr snp_matrix = MatrixXdr::Zero(n_samples, 1);
for (int i = 0; i < block_size; i++) {
read_snp(bed_ifs, n_samples, global_snp_index, metadata, snp_matrix);
if (genotype_mask_matrix(i, 0) == 1) {
encode_snp(genotype_block, snp_matrix);
}
}


    MatrixXdr data = readCSVToMatrixXdr(test_csv);
    MatrixXdr matrix_block = data.block(0, 0, n_samples, block_size);
    MatrixXdr masked_matrix_block(n_samples, mask_size);
    int count = 0;
    for (int i = 0; i < block_size; i++) {
      if (genotype_mask[i] == 1) {
        masked_matrix_block.col(count) = matrix_block.col(i);
        count++;
      }
    }

    MatrixXdr fame_means(block_size, 1);
    MatrixXdr fame_stds(block_size, 1);
    compute_block_stats(genotype_block, fame_means, fame_stds, n_samples,
            mask_size);

    MatrixXdr TWO = MatrixXdr::Ones(n_samples, mask_size) * 2;
    masked_matrix_block = TWO - masked_matrix_block;
    MatrixXdr means(mask_size, 1);
    MatrixXdr stds(mask_size, 1);
    means = masked_matrix_block.colwise().mean().transpose();
    for (int i = 0; i < mask_size; i++) {
      stds(i, 0) = 1 / sqrt((means(i, 0) * (1 - (0.5 * means(i, 0)))));
    }
    for (int i = 0; i < mask_size; i++) {
      masked_matrix_block.col(i) =
          masked_matrix_block.col(i).array() - means(i, 0);
      masked_matrix_block.col(i) =
          masked_matrix_block.col(i).array() * stds(i, 0);
    }

    MatrixXdr pheno_mask;
    MatrixXdr pheno;
    read_phenotypes(n_samples, test_pheno, pheno, pheno_mask);

    MatrixXdr Xy(mask_size, n_samples);
    MatrixXdr XXy_expected(n_samples, 1);
    Xy = masked_matrix_block.transpose() * pheno;
    XXy_expected = masked_matrix_block * Xy;

    // when

    MatrixXdr XXy_observed = MatrixXdr::Zero(n_samples, 1);

    XXy_observed.col(0) +=
        compute_XXz(mask_size, pheno, fame_means, fame_stds, pheno_mask,
                    n_randvecs, n_samples, genotype_block,
                mask_size, focal_snp_local_index, false);

    double abs_error = (XXy_expected - XXy_observed).array().abs().sum();

    // then
    expect_true(abs_error < tolerance);
  }
}

context("C++ test focal SNP exclusion") {
  test_that("excluding focal SNP gives the expected result for XXz") {
    // given
    int n_variance_components = 1;
    int focal_snp_local_index = 0;
    int n_randvecs = 1;
    genotype genotype_block;

    std::vector<int> genotype_mask = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int mask_size = genotype_mask.size();
    MatrixXdr genotype_mask_matrix = MatrixXdr::Zero(block_size, 1);
    for (int i = 0; i < block_size; i++) {
      genotype_mask_matrix(i, 0) = genotype_mask[i];
    }

    set_block_parameters(genotype_block, n_samples, block_size);
    std::ifstream bed_ifs(test_bed.c_str(), ios::in | ios::binary);
    int global_snp_index = -1;
    read_genotype_block(bed_ifs, block_size, genotype_block, n_samples,
                        global_snp_index, metadata);

    MatrixXdr data = readCSVToMatrixXdr(test_csv);

    MatrixXdr fame_means(block_size, 1);
    MatrixXdr fame_stds(block_size, 1);
    compute_block_stats(genotype_block, fame_means, fame_stds, n_samples,
                        block_size);

    MatrixXdr matrix_block = data.block(0, 1, n_samples, block_size - 1);
    MatrixXdr TWO = MatrixXdr::Ones(n_samples, block_size - 1) * 2;
    matrix_block = TWO - matrix_block;
    MatrixXdr means(block_size - 1, 1);
    MatrixXdr stds(block_size - 1, 1);
    means = matrix_block.colwise().sum().transpose() / n_samples;
    // compute stds_minor by iterating over the columns
    for (int i = 0; i < block_size - 1; i++) {
      stds(i, 0) = 1 / sqrt((means(i, 0) * (1 - (0.5 * means(i, 0)))));
    }
    for (int i = 0; i < block_size - 1; i++) {
      matrix_block.col(i) = matrix_block.col(i).array() - means(i, 0);
      matrix_block.col(i) = matrix_block.col(i).array() * stds(i, 0);
    }

    MatrixXdr pheno_mask;
    MatrixXdr pheno;
    read_phenotypes(n_samples, test_pheno, pheno, pheno_mask);

    MatrixXdr Xy(block_size - 1, n_samples);
    MatrixXdr XXy_expected(n_samples, 1);
    Xy = matrix_block.transpose() * pheno;
    XXy_expected = matrix_block * Xy;

    // when
    // memory setup for mailman
    double *partialsums = new double[0];
    double *sum_op;
    double *yint_e;
    double *yint_m;
    double **y_e;
    double **y_m;
    allocate_memory(n_randvecs, genotype_block, partialsums, sum_op, yint_e,
                    yint_m, y_e, y_m);

    MatrixXdr XXy_observed = MatrixXdr::Zero(n_samples, 1);
    //    XXy_observed.col(0) +=
    //        compute_XXy(block_size, pheno, fame_means, fame_stds, pheno_mask,
    //                    genotype_mask_matrix, n_samples,
    //                    sum_op, genotype_block, yint_m, y_m, block_size,
    //                    yint_e, y_e, partialsums, focal_snp_local_index,
    //                    false);
    bool exclude_sel_snp = true;
    XXy_observed.col(0) +=
        compute_XXz(block_size, pheno, fame_means, fame_stds, pheno_mask,
                    n_randvecs, n_samples, genotype_block,
                    block_size, focal_snp_local_index, exclude_sel_snp);

    double abs_error = (XXy_expected - XXy_observed).array().abs().sum();

    // then
    expect_true(abs_error < tolerance);
  }

  test_that("excluding focal SNP gives the expected result for yXXy") {
    // given
    int n_variance_components = 1;
    int focal_snp_local_index = 0;
    int n_randvecs = 1;
    genotype genotype_block;

    std::vector<int> genotype_mask = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int mask_size = genotype_mask.size();
    MatrixXdr genotype_mask_matrix = MatrixXdr::Zero(block_size, 1);
    for (int i = 0; i < block_size; i++) {
      genotype_mask_matrix(i, 0) = genotype_mask[i];
    }

    set_block_parameters(genotype_block, n_samples, block_size);
    std::ifstream bed_ifs(test_bed.c_str(), ios::in | ios::binary);
    int global_snp_index = -1;
    read_genotype_block(bed_ifs, block_size, genotype_block, n_samples,
                        global_snp_index, metadata);

    MatrixXdr data = readCSVToMatrixXdr(test_csv);

    MatrixXdr fame_means(block_size, 1);
    MatrixXdr fame_stds(block_size, 1);
    compute_block_stats(genotype_block, fame_means, fame_stds, n_samples,
                        block_size);

    MatrixXdr matrix_block = data.block(0, 1, n_samples, block_size - 1);
    MatrixXdr TWO = MatrixXdr::Ones(n_samples, block_size - 1) * 2;
    matrix_block = TWO - matrix_block;
    MatrixXdr means(block_size - 1, 1);
    MatrixXdr stds(block_size - 1, 1);
    means = matrix_block.colwise().sum().transpose() / n_samples;
    // compute stds_minor by iterating over the columns
    for (int i = 0; i < block_size - 1; i++) {
      stds(i, 0) = 1 / sqrt((means(i, 0) * (1 - (0.5 * means(i, 0)))));
    }
    for (int i = 0; i < block_size - 1; i++) {
      matrix_block.col(i) = matrix_block.col(i).array() - means(i, 0);
      matrix_block.col(i) = matrix_block.col(i).array() * stds(i, 0);
    }

    MatrixXdr pheno_mask;
    MatrixXdr pheno;
    read_phenotypes(n_samples, test_pheno, pheno, pheno_mask);

    // when
    // memory setup for mailman
    double *partialsums = new double[0];
    double *sum_op;
    double *yint_e;
    double *yint_m;
    double **y_e;
    double **y_m;
    allocate_memory(n_randvecs, genotype_block, partialsums, sum_op, yint_e,
                    yint_m, y_e, y_m);

    MatrixXdr Xy(block_size - 1, n_samples);
    MatrixXdr yXXy_expected(1, 1);
    Xy = matrix_block.transpose() * pheno;
    yXXy_expected(0, 0) = (Xy.array() * Xy.array()).sum();

    bool exclude_sel_snp = true;
    MatrixXdr yXXy_observed;
    yXXy_observed = MatrixXdr::Zero(n_variance_components, 1);
    yXXy_observed(0, 0) += compute_yXXy(
        block_size, pheno, fame_means, fame_stds, focal_snp_local_index,
        genotype_block, block_size, exclude_sel_snp);

    double abs_error = (yXXy_expected - yXXy_observed).array().abs().sum();

    // then
    expect_true(abs_error < tolerance);
  }
}