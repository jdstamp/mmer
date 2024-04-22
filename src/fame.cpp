/*  famer: An R Package implementation of the Fast Marginal Epistasis test
 *  Copyright (C) 2024  Julian Stamp
 *  This code is licensed under MIT license (see LICENSE.md for details)
 */

#include "fame.h"
#include "allocate_memory.h"
#include "computation.h"
#include "compute_block_stats.h"
#include "compute_covariance_q.h"
#include "compute_mom_components.h"
#include "count_data.h"
#include "fit_covariates.h"
#include "initialize_random_vectors.h"
#include "read_covariates.h"
#include "read_genotypes.h"
#include "read_phenotypes.h"
#include <algorithm>
#include <highfive/H5Easy.hpp>
#include <random>

#include "set_block_parameters.h"

void read_genotype_mask(const std::string &genotype_mask_file, int n_snps,
                        int gxg_i, const string &gxg_h5_dataset,
                        MatrixXdr &genotype_mask, int &n_gxg_snps);

MatrixXdr indices_to_binary(const std::vector<int> &indices, int length) {
  MatrixXdr binaryVec = MatrixXdr::Zero(length, 1);
  for (int index : indices) {
    if (index >= 0 && index < length) {
      binaryVec(index, 0) = 1;
    }
  }
  return binaryVec;
}

// [[Rcpp::export]]
Rcpp::List fame_cpp(std::string plink_file, std::string pheno_file,
                    std::string covariate_file, int n_randvecs, int n_blocks,
                    int rand_seed, std::vector<int> gxg_indices,
                    std::string genotype_mask_file) {

  // Mailman algo variables.
  // TODO: Give meaningful names or move into mailman algo file
  double *partialsums = new double[0];
  double *sum_op;
  double *yint_e;
  double *yint_m;
  double **y_e;
  double **y_m;

  auto start_grm = std::chrono::high_resolution_clock::now();
  int n_variance_components = 2; // For now limited to GRM and GXG
  // initialize object to collect point_est and cov_sigma
  MatrixXdr VC;
  MatrixXdr SE;

  std::stringstream bim_stream;
  std::stringstream fam_stream;
  std::stringstream bed_stream;
  bim_stream << plink_file << ".bim";
  fam_stream << plink_file << ".fam";
  bed_stream << plink_file << ".bed";
  std::string bim_file = bim_stream.str();
  std::string fam_file = fam_stream.str();
  std::string bed_file = bed_stream.str();

  MatrixXdr pheno_mask;
  MatrixXdr pheno;
  MatrixXdr XXz; // Matrix for GRM related estimates
  genotype genotype_block;

  MatrixXdr allelecount_means;
  MatrixXdr allelecount_stds;

  int n_snps = count_snps_bim(bim_file);
  int n_samples = count_samples(pheno_file);
  int n_fam_lines = count_fam(fam_file);

  int step_size = n_snps / n_blocks;
  int step_size_remainder = n_snps % n_blocks;
  std::vector<int> block_sizes(n_blocks, step_size);
  block_sizes.back() += step_size_remainder; // add remainder to last block

  if (n_fam_lines != n_samples) {
    throw std::runtime_error(
        "Number of samples in fam file and pheno file do not match.");
  } // TODO: this is already done in the R part. do we need it here?

  VC.resize(gxg_indices.size(), n_variance_components + 1);
  SE.resize(gxg_indices.size(), n_variance_components + 1);

  // read data into matrices pheno and pheno_mask
  read_phenotypes(n_samples, pheno_file, pheno, pheno_mask);

  // Covariate handling - needs cleanup
  double y_sum = 0;
  double y_mean = 0;
  if (covariate_file != "") {
    bool snp_fix_ef = false;
    MatrixXdr covariate;
    int cov_num;
    cov_num = read_covariates(false, n_samples, covariate_file, covariate,
                              snp_fix_ef);
    //        if (snp_fix_ef == true) {
    //            covariate.col(cov_num - 1) = focal_snp_gtype.col(0);
    //        }

    pheno = fit_covariates(pheno_mask, pheno, n_samples, y_sum, y_mean,
                           covariate, cov_num);

  } else {
    y_sum = pheno.sum();
    y_mean = y_sum / pheno_mask.sum();
    for (int i = 0; i < n_samples; i++) {
      if (pheno(i, 0) != 0) {
        pheno(i, 0) = pheno(i, 0) - y_mean;
      } // center phenotype
    }
  }

  // inserted below
  //    MatrixXdr XXz_observed;
  //    MatrixXdr random_vectors;
  //
  //    XXz = MatrixXdr::Zero(n_samples, n_randvecs);
  //    random_vectors = initialize_random_vectors(n_randvecs, rand_seed,
  //    pheno_mask,
  //                                             random_vectors, n_samples);

  metaData metadata = set_metadata(n_samples, n_snps);
  ifstream bed_ifs(bed_file.c_str(), ios::in | ios::binary);
  int global_snp_index = -1;
  // inserted below
  //  for (int block_index = 0; block_index < n_blocks; block_index++) {
  //
  //    int block_size = block_sizes[block_index];
  //
  //    set_block_parameters(genotype_block, n_samples, block_size);
  //
  //      read_genotype_block(bed_ifs, block_size, genotype_block, n_samples,
  //                          global_snp_index, metadata);
  //
  //    if (block_size != 0) {
  //      compute_block_stats(genotype_block, allelecount_means,
  //      allelecount_stds,
  //                          n_samples, block_size);
  //
  //      allocate_memory(n_randvecs, genotype_block, partialsums, sum_op,
  //      yint_e,
  //                      yint_m, y_e, y_m);
  //
  //      XXz_observed = compute_XXz(block_size, random_vectors,
  //      allelecount_means,
  //                             allelecount_stds, pheno_mask, n_randvecs,
  //                             n_samples, 0, sum_op, genotype_block, yint_m,
  //                             y_m, block_size, yint_e, y_e, partialsums,
  //                             false);
  //
  //      for (int z_index = 0; z_index < n_randvecs; z_index++) {
  //        XXz.col(z_index) += XXz_observed.col(z_index);
  //      }
  //      deallocate_memory(partialsums, sum_op, yint_e, yint_m, y_e, y_m,
  //                        genotype_block);
  //    }
  //  }
  //  auto end_grm = std::chrono::high_resolution_clock::now();
  //  std::chrono::duration<double> elapsed_grm = end_grm - start_grm;
  //  std::cout << "Execution time of grm: " << elapsed_grm.count() << "
  //  seconds."
  //            << std::endl;
  auto start_gxg = std::chrono::high_resolution_clock::now();

  // TODO: parallel loop starts here?
  int process_count = -1;
  for (int gxg_i : gxg_indices) {
    process_count++;
    //    std::cout << "gxg_i: " << gxg_i << std::endl;
    MatrixXdr focal_snp_gtype;
    // insert declaration GRM
    MatrixXdr temp_grm;
    MatrixXdr random_vectors;

    XXz = MatrixXdr::Zero(n_samples, n_randvecs);
    random_vectors = initialize_random_vectors(
        n_randvecs, rand_seed, pheno_mask, random_vectors, n_samples);
    // end insert

    // experimental masking
    string gxg_h5_dataset = "gxg";
    MatrixXdr genotype_mask;
    int n_gxg_snps;
    // read mask file into genotype_mask and n_gxg_snps
    read_genotype_mask(genotype_mask_file, n_snps, gxg_i, gxg_h5_dataset,
                       genotype_mask, n_gxg_snps);

    // end experimental masking
    vector<int> n_snps_variance_component = {n_snps, n_gxg_snps};

    MatrixXdr gxg_random_vectors;
    MatrixXdr GxGz;
    MatrixXdr yXXy;

    MatrixXdr collect_XXy;
    MatrixXdr collect_XXUy;

    focal_snp_gtype.resize(n_samples, 1);

    global_snp_index = -1;
    read_focal_snp(bed_file, focal_snp_gtype, gxg_i, n_samples, n_snps,
                   global_snp_index);
    normalize_genotype(focal_snp_gtype, n_samples);

    int focal_snp_block = std::min(gxg_i / step_size, n_blocks - 1);

    int focal_snp_local_index = gxg_i - (step_size * focal_snp_block);

    random_vectors = initialize_random_vectors(
        n_randvecs, rand_seed, pheno_mask, random_vectors, n_samples);

    gxg_random_vectors =
        random_vectors.array().colwise() * focal_snp_gtype.col(0).array();
    MatrixXdr temp_gxg;
    GxGz = MatrixXdr::Zero(n_samples, n_randvecs);
    yXXy = MatrixXdr::Zero(n_variance_components, 1);

    collect_XXy = MatrixXdr::Zero(n_samples, n_variance_components + 1);
    collect_XXUy = MatrixXdr::Zero(n_samples, (n_variance_components + 1) *
                                                  (n_variance_components + 1));

    bed_ifs.seekg(0, std::ios::beg); // reset file pointer to beginning
    global_snp_index = -1;

    for (int block_index = 0; block_index < n_blocks; block_index++) {

      int block_size = block_sizes[block_index];
      // initialize gxg_mask
      MatrixXdr gxg_mask =
          genotype_mask.block(block_index * block_size, 0, block_size, 1);
      MatrixXdr grm_mask = MatrixXdr::Ones(block_size, 1);

      set_block_parameters(genotype_block, n_samples, block_size);

      read_genotype_block(bed_ifs, block_size, genotype_block, n_samples,
                          global_snp_index, metadata);

      if (block_size != 0) {
        compute_block_stats(genotype_block, allelecount_means, allelecount_stds,
                            n_samples, block_size);

        allocate_memory(n_randvecs, genotype_block, partialsums, sum_op, yint_e,
                        yint_m, y_e, y_m);
        // insert GRM
        temp_grm = compute_XXz(
            block_size, random_vectors, allelecount_means, allelecount_stds,
            pheno_mask, grm_mask, n_randvecs, n_samples, sum_op, genotype_block,
            yint_m, y_m, block_size, yint_e, y_e, partialsums, 0, false);

        for (int z_index = 0; z_index < n_randvecs; z_index++) {
          XXz.col(z_index) += temp_grm.col(z_index);
        }
        // end insert GRM
        bool in_gxg_block = (focal_snp_block == block_index);
        temp_gxg = compute_XXz(
            block_size, gxg_random_vectors, allelecount_means, allelecount_stds,
            pheno_mask, gxg_mask, n_randvecs, n_samples, sum_op, genotype_block,
            yint_m, y_m, block_size, yint_e, y_e, partialsums,
            focal_snp_local_index, in_gxg_block);
        temp_gxg = temp_gxg.array().colwise() * focal_snp_gtype.col(0).array();

        for (int z_index = 0; z_index < n_randvecs; z_index++) {
          GxGz.col(z_index) += temp_gxg.col(z_index);
        }
        collect_XXy.col(0) += compute_XXz(
            block_size, pheno, allelecount_means, allelecount_stds, pheno_mask,
            grm_mask, 1, n_samples, sum_op, genotype_block, yint_m, y_m,
            block_size, yint_e, y_e, partialsums, focal_snp_local_index, false);

        MatrixXdr gxg_pheno;
        gxg_pheno = pheno.array() * focal_snp_gtype.col(0).array();
        MatrixXdr temp_Gy =
            compute_XXz(block_size, gxg_pheno, allelecount_means,
                        allelecount_stds, pheno_mask, gxg_mask, 1, n_samples,
                        sum_op, genotype_block, yint_m, y_m, block_size, yint_e,
                        y_e, partialsums, focal_snp_local_index, in_gxg_block);
        temp_Gy = temp_Gy.array() * focal_snp_gtype.col(0).array();
        collect_XXy.col(1) += temp_Gy;

        yXXy(0, 0) +=
            compute_yXXy(block_size, pheno, allelecount_means, allelecount_stds,
                         focal_snp_local_index, sum_op, genotype_block,
                         grm_mask, yint_m, y_m, block_size, partialsums, false);

        yXXy(1, 0) += compute_yXXy(block_size, gxg_pheno, allelecount_means,
                                   allelecount_stds, focal_snp_local_index,
                                   sum_op, genotype_block, gxg_mask, yint_m,
                                   y_m, block_size, partialsums, in_gxg_block);
        deallocate_memory(partialsums, sum_op, yint_e, yint_m, y_e, y_m,
                          genotype_block);
      }
    }

    for (int i = 0; i < n_variance_components; i++) {
      collect_XXy.col(i) = collect_XXy.col(i) / n_snps_variance_component[i];
    }
    collect_XXy.col(n_variance_components) = pheno;

    bed_ifs.seekg(0, std::ios::beg); // reset file pointer to beginning
    global_snp_index = -1;

    for (int block_index = 0; block_index < n_blocks; block_index++) {

      int block_size = block_sizes[block_index];
      // initialize gxg_mask
      MatrixXdr gxg_mask =
          genotype_mask.block(block_index * block_size, 0, block_size, 1);
      MatrixXdr grm_mask = MatrixXdr::Ones(block_size, 1);

      set_block_parameters(genotype_block, n_samples, block_size);

      read_genotype_block(bed_ifs, block_size, genotype_block, n_samples,
                          global_snp_index, metadata);

      if (block_size != 0) {
        compute_block_stats(genotype_block, allelecount_means, allelecount_stds,
                            n_samples, block_size);

        allocate_memory(n_randvecs, genotype_block, partialsums, sum_op, yint_e,
                        yint_m, y_e, y_m);

        MatrixXdr scaled_vec;
        for (int i = 0; i < (n_variance_components + 1); i++) {
          MatrixXdr temp_XXUy = compute_XXz(
              block_size, collect_XXy.col(i), allelecount_means,
              allelecount_stds, pheno_mask, grm_mask, 1, n_samples, sum_op,
              genotype_block, yint_m, y_m, block_size, yint_e, y_e, partialsums,
              focal_snp_local_index, false);

          scaled_vec =
              collect_XXy.col(i).array() * focal_snp_gtype.col(0).array();
          MatrixXdr temp_GUy = compute_XXz(
              block_size, scaled_vec, allelecount_means, allelecount_stds,
              pheno_mask, gxg_mask, 1, n_samples, sum_op, genotype_block,
              yint_m, y_m, block_size, yint_e, y_e, partialsums,
              focal_snp_local_index, false);
          temp_GUy = temp_GUy.array() * focal_snp_gtype.col(0).array();

          collect_XXUy.col(i) += temp_XXUy / n_snps_variance_component[0];
          collect_XXUy.col(((n_variance_components + 1)) + i) +=
              temp_GUy / n_snps_variance_component[1];
        }

        deallocate_memory(partialsums, sum_op, yint_e, yint_m, y_e, y_m,
                          genotype_block);
      }
    }

    for (int i = 0; i < (n_variance_components + 1); i++) {
      collect_XXUy.col((n_variance_components * (n_variance_components + 1)) +
                       i) = collect_XXy.col(i);
    }

    int n_samples_mask = pheno_mask.sum();

    // naming from eq. (5) in mvMAPIT paper
    MatrixXdr S(n_variance_components + 1, n_variance_components + 1);
    MatrixXdr q(n_variance_components + 1, 1);
    compute_mom_components(n_randvecs, n_variance_components, pheno,
                           random_vectors, XXz, GxGz, yXXy,
                           n_snps_variance_component, n_samples_mask, S, q);

    MatrixXdr point_est = S.colPivHouseholderQr().solve(q);

    MatrixXdr cov_q;
    compute_covariance_q(n_variance_components, collect_XXUy, point_est, cov_q);

    MatrixXdr invS = S.inverse();

    MatrixXdr cov_sigma = invS * cov_q * invS;

    VC.row(process_count) = point_est.col(0).transpose();
    SE.row(process_count) = cov_sigma.diagonal().array().sqrt().transpose();
  }

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_gxg = end - start_gxg;
  //  std::cout << "Execution time of gxg: " << elapsed_gxg.count() << "
  //  seconds."
  //            << std::endl;

  return Rcpp::List::create(Rcpp::Named("vc_estimate") = VC,
                            Rcpp::Named("vc_se") = SE);
}

void read_genotype_mask(const std::string &genotype_mask_file, int n_snps,
                        int gxg_i, const string &gxg_h5_dataset,
                        MatrixXdr &genotype_mask, int &n_gxg_snps) {
  std::vector<int> genotype_mask_indices;
  if (genotype_mask_file != "") {
    H5Easy::File file(genotype_mask_file, H5Easy::File::ReadOnly);
    genotype_mask_indices = H5Easy::load<std::vector<int>>(
        file, gxg_h5_dataset + "/" + std::to_string(gxg_i));
    genotype_mask = indices_to_binary(genotype_mask_indices, n_snps);
    n_gxg_snps = genotype_mask_indices.size();
  } else {
    genotype_mask = MatrixXdr::Ones(n_snps, 1);
    n_gxg_snps = n_snps - 1;
  }
}
