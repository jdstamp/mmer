/*  famer: An R Package implementation of the Fast Marginal Epistasis test
 *  Copyright (C) 2024  Julian Stamp
 *  This code is licensed under MIT license (see LICENSE.md for details)
 */

// [[Rcpp::plugins(openmp)]]
#include <omp.h>

#include "computation.h"
#include "compute_block_stats.h"
#include "compute_covariance_q.h"
#include "compute_mom_components.h"
#include "count_data.h"
#include "fame.h"
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
                    std::string genotype_mask_file, int n_threads) {

  auto start = std::chrono::high_resolution_clock::now();

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif
  //  #pragma omp parallel for schedule(dynamic)
  //  1. try oscar because linux
  //  a. maybe get stack trace of crash?
  // 2. provide bash script for parallelization accross jobs

  // TODO: make "gxg" configurable
  string gxg_h5_dataset = "gxg";

  int n_variance_components = 2; // For now limited to GRM and GXG
  // initialize object to collect point_est and cov_sigma
  MatrixXdr VC;
  MatrixXdr SE;

  std::stringstream bim_stream;
  std::stringstream bed_stream;
  bim_stream << plink_file << ".bim";
  bed_stream << plink_file << ".bed";
  std::string bim_file = bim_stream.str();
  std::string bed_file = bed_stream.str();
  int n_gxg_idx = gxg_indices.size();

  MatrixXdr pheno_mask;
  MatrixXdr pheno;
  MatrixXdr XXz;
  MatrixXdr GxGz;
  MatrixXdr yXXy;
  MatrixXdr yGxGy;
  MatrixXdr collect_XXy;
  MatrixXdr collect_Gy;
  MatrixXdr collect_XXUy;

  MatrixXdr random_vectors;

  genotype grm_genotype_block;
  std::vector<genotype> gxg_genotype_blocks(n_gxg_idx);

  MatrixXdr allelecount_means;
  MatrixXdr allelecount_stds;

  int n_snps = count_snps_bim(bim_file);
  int n_samples = count_samples(pheno_file);

  int step_size = n_snps / n_blocks;
  int step_size_remainder = n_snps % n_blocks;
  std::vector<int> block_sizes(n_blocks, step_size);
  block_sizes.back() += step_size_remainder; // add remainder to last block

  VC.resize(n_gxg_idx, n_variance_components + 1);
  SE.resize(n_gxg_idx, n_variance_components + 1);

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
    // TODO: make this a default true and remove the option ? this might have
    //  to move to the SNP specific portion if we do that
    //        if (snp_fix_ef == true) {
    //            covariate.col(cov_num - 1) = focal_snp_gtype.col(0);
    //        }

    pheno = fit_covariates(pheno_mask, pheno, n_samples, y_sum, y_mean,
                           covariate, cov_num);

  } else {
    y_sum = pheno.sum();
    y_mean = y_sum / pheno_mask.sum();
    // TODO: should we scale to 1 as well?
    for (int i = 0; i < n_samples; i++) {
      if (pheno(i, 0) != 0) {
        pheno(i, 0) = pheno(i, 0) - y_mean;
      } // center phenotype
    }
  }

  vector<int> n_mask_gxg(n_gxg_idx);
  XXz = MatrixXdr::Zero(n_samples, n_randvecs);
  GxGz = MatrixXdr::Zero(n_samples, n_randvecs * n_gxg_idx);
  yXXy = MatrixXdr::Zero(1, 1);
  yGxGy = MatrixXdr::Zero(n_gxg_idx, 1);
  random_vectors = initialize_random_vectors(n_randvecs, rand_seed, pheno_mask,
                                             random_vectors, n_samples);

  collect_Gy = MatrixXdr::Zero(n_samples, n_gxg_idx);
  collect_XXUy =
      MatrixXdr::Zero(n_samples, (n_variance_components + 1) *
                                     (n_variance_components + 1) * n_gxg_idx);

  metaData metadata = set_metadata(n_samples, n_snps);
  ifstream bed_ifs(bed_file.c_str(), ios::in | ios::binary);
  int global_snp_index = -1;

  XXz = MatrixXdr::Zero(n_samples, n_randvecs);
  yXXy = MatrixXdr::Zero(1, 1);
  collect_XXy = MatrixXdr::Zero(n_samples, 1);
  MatrixXdr temp_grm;

  std::vector<int> n_gxg_snps_list(n_gxg_idx);
  MatrixXdr binary_gxg_mask = MatrixXdr::Ones(n_snps, n_gxg_idx);

  random_vectors = initialize_random_vectors(n_randvecs, rand_seed, pheno_mask,
                                             random_vectors, n_samples);

  bed_ifs.seekg(0, std::ios::beg); // reset file pointer to beginning
  global_snp_index = -1;

  for (int block_index = 0; block_index < n_blocks; block_index++) {
    Rcpp::checkUserInterrupt();
    int block_size = block_sizes[block_index];
    MatrixXdr snp_matrix = MatrixXdr::Zero(n_samples, 1);
    MatrixXdr grm_mask = MatrixXdr::Ones(block_size, 1);
    set_block_parameters(grm_genotype_block, n_samples, block_size);

    for (int parallel_idx = 0; parallel_idx < n_gxg_idx; parallel_idx++) {
      int n_gxg_snps;
      MatrixXdr genotype_mask;
      int gxg_i = gxg_indices[parallel_idx];
      read_genotype_mask(genotype_mask_file, n_snps, gxg_i, gxg_h5_dataset,
                         genotype_mask, n_gxg_snps);
      MatrixXdr gxg_mask =
          genotype_mask.block(block_index * block_size, 0, block_size, 1);
      int gxg_snps_in_block = gxg_mask.sum();
      set_block_parameters(gxg_genotype_blocks[parallel_idx], n_samples,
                           gxg_snps_in_block);
      binary_gxg_mask.col(parallel_idx) = genotype_mask;
      n_gxg_snps_list[parallel_idx] = n_gxg_snps;
    }

    for (int i = 0; i < block_size; i++) {
      read_snp(bed_ifs, n_samples, global_snp_index, metadata, snp_matrix);
      encode_snp(grm_genotype_block, snp_matrix);

#pragma omp parallel for schedule(dynamic)
      for (int parallel_idx = 0; parallel_idx < n_gxg_idx;
           parallel_idx++) { // parallel loop 1

        MatrixXdr gxg_mask =
            binary_gxg_mask.col(parallel_idx)
                .block(block_index * block_size, 0, block_size, 1);
        if (gxg_mask(i, 0) == 1) {
          encode_snp(gxg_genotype_blocks[parallel_idx], snp_matrix);
        }
      } // end of parallel loop 1
    }

    if (block_size != 0) {
      compute_block_stats(grm_genotype_block, allelecount_means,
                          allelecount_stds, n_samples, block_size);

      temp_grm = compute_XXz(random_vectors, allelecount_means,
                             allelecount_stds, pheno_mask, n_randvecs,
                             n_samples, grm_genotype_block, 0, false);

      for (int z_index = 0; z_index < n_randvecs; z_index++) {
        XXz.col(z_index) += temp_grm.col(z_index);
      }

      yXXy(0, 0) += compute_yXXy(pheno, allelecount_means, allelecount_stds, 0,
                                 grm_genotype_block, false);
      collect_XXy.col(0) +=
          compute_XXz(pheno, allelecount_means, allelecount_stds, pheno_mask, 1,
                      n_samples, grm_genotype_block, 0, false);

      grm_genotype_block.clear_block();

#pragma omp parallel for schedule(dynamic)
      for (int parallel_idx = 0; parallel_idx < n_gxg_idx; parallel_idx++) {
        // parallel loop 2

        MatrixXdr focal_snp_gtype;
        MatrixXdr gxg_allelecount_means;
        MatrixXdr gxg_allelecount_stds;
        int n_gxg_snps = n_gxg_snps_list[parallel_idx];
        MatrixXdr genotype_mask;
        int gxg_i = gxg_indices[parallel_idx];
        MatrixXdr gxg_mask =
            binary_gxg_mask.col(parallel_idx)
                .block(block_index * block_size, 0, block_size, 1);
        int gxg_snps_in_block = gxg_mask.sum();
        n_mask_gxg[parallel_idx] = n_gxg_snps;
        vector<int> n_snps_variance_component = {n_snps, n_gxg_snps};

        focal_snp_gtype.resize(n_samples, 1);
        int gindex = -1;
        read_focal_snp(bed_file, focal_snp_gtype, gxg_i, n_samples, n_snps,
                       gindex);
        normalize_genotype(focal_snp_gtype, n_samples);

        int focal_snp_block = std::min(gxg_i / step_size, n_blocks - 1);
        int focal_snp_local_index = gxg_i - (step_size * focal_snp_block);

        MatrixXdr gxg_random_vectors =
            random_vectors.array().colwise() * focal_snp_gtype.col(0).array();
        MatrixXdr temp_gxg;
        compute_block_stats(gxg_genotype_blocks[parallel_idx],
                            gxg_allelecount_means, gxg_allelecount_stds,
                            n_samples, gxg_snps_in_block);

        bool in_gxg_block = (focal_snp_block == block_index);
        temp_gxg = compute_XXz(gxg_random_vectors, gxg_allelecount_means,
                               gxg_allelecount_stds, pheno_mask, n_randvecs,
                               n_samples, gxg_genotype_blocks[parallel_idx],
                               focal_snp_local_index, in_gxg_block);

        temp_gxg = temp_gxg.array().colwise() * focal_snp_gtype.col(0).array();

        for (int z_index = 0; z_index < n_randvecs; z_index++) {
          GxGz.col(z_index + parallel_idx * n_randvecs) +=
              temp_gxg.col(z_index);
        }

        MatrixXdr gxg_pheno;
        gxg_pheno = pheno.array() * focal_snp_gtype.col(0).array();
        MatrixXdr temp_Gy = compute_XXz(
            gxg_pheno, gxg_allelecount_means, gxg_allelecount_stds, pheno_mask,
            1, n_samples, gxg_genotype_blocks[parallel_idx],
            focal_snp_local_index, in_gxg_block);
        temp_Gy = temp_Gy.array() * focal_snp_gtype.col(0).array();
        collect_Gy.col(parallel_idx) += temp_Gy;

        yGxGy(parallel_idx, 0) +=
            compute_yXXy(gxg_pheno, gxg_allelecount_means, gxg_allelecount_stds,
                         focal_snp_local_index,
                         gxg_genotype_blocks[parallel_idx], in_gxg_block);
        gxg_genotype_blocks[parallel_idx].clear_block();
      } // end of parallel loop 2
    }
  }

  collect_XXy = collect_XXy / n_snps;

  // TODO: enable this when loop is in block
  for (int i = 0; i < n_gxg_idx; i++) {
    collect_Gy.col(i) = collect_Gy.col(i) / n_mask_gxg[i];
  }

  bed_ifs.seekg(0, std::ios::beg); // reset file pointer to beginning
  global_snp_index = -1;

  for (int block_index = 0; block_index < n_blocks; block_index++) {
    Rcpp::checkUserInterrupt();
    int block_size = block_sizes[block_index];
    MatrixXdr snp_matrix = MatrixXdr::Zero(n_samples, 1);
    MatrixXdr grm_mask = MatrixXdr::Ones(block_size, 1);
    set_block_parameters(grm_genotype_block, n_samples, block_size);

#pragma omp parallel for schedule(dynamic)
    for (int parallel_idx = 0; parallel_idx < n_gxg_idx; parallel_idx++) {
      // parallel loop 3

      // initialize gxg_mask
      int n_gxg_snps;
      MatrixXdr genotype_mask;
      int gxg_i = gxg_indices[parallel_idx];

      MatrixXdr gxg_mask =
          binary_gxg_mask.col(parallel_idx)
              .block(block_index * block_size, 0, block_size, 1);
      int gxg_snps_in_block = gxg_mask.sum();
      set_block_parameters(gxg_genotype_blocks[parallel_idx], n_samples,
                           gxg_snps_in_block);
    } // end loop 3

    for (int i = 0; i < block_size; i++) {
      read_snp(bed_ifs, n_samples, global_snp_index, metadata, snp_matrix);
      encode_snp(grm_genotype_block, snp_matrix);

#pragma omp parallel for schedule(dynamic)
      for (int parallel_idx = 0; parallel_idx < n_gxg_idx;
           parallel_idx++) { // parallel loop 4
        // initialize gxg_mask
        int n_gxg_snps;
        MatrixXdr genotype_mask;
        int gxg_i = gxg_indices[parallel_idx];
        MatrixXdr gxg_mask =
            binary_gxg_mask.col(parallel_idx)
                .block(block_index * block_size, 0, block_size, 1);
        if (gxg_mask(i, 0) == 1) {
          encode_snp(gxg_genotype_blocks[parallel_idx], snp_matrix);
        }
      } // end of parallel loop 4
    }

    if (block_size != 0) {
      compute_block_stats(grm_genotype_block, allelecount_means,
                          allelecount_stds, n_samples, block_size);

      MatrixXdr temp_XXXXy =
          compute_XXz(collect_XXy.col(0), allelecount_means, allelecount_stds,
                      pheno_mask, 1, n_samples, grm_genotype_block, 0, false) /
          n_snps;

#pragma omp parallel for schedule(dynamic)
      for (int parallel_idx = 0; parallel_idx < n_gxg_idx;
           parallel_idx++) { // parallel loop 5
        MatrixXdr temp_XXUy = compute_XXz(
            collect_Gy.col(parallel_idx), allelecount_means, allelecount_stds,
            pheno_mask, 1, n_samples, grm_genotype_block, 0, false);

        collect_XXUy.col(0 + (n_variance_components + 1) *
                                 (n_variance_components + 1) * parallel_idx) +=
            temp_XXXXy;
        collect_XXUy.col(1 + (n_variance_components + 1) *
                                 (n_variance_components + 1) * parallel_idx) +=
            temp_XXUy / n_snps; // n_snps_variance_component[0];

      } // end of parallel loop 5

      grm_genotype_block.clear_block();

#pragma omp parallel for schedule(dynamic)
      for (int parallel_idx = 0; parallel_idx < n_gxg_idx;
           parallel_idx++) { // parallel loop 6
        collect_XXUy.col(2 + (n_variance_components + 1) *
                                 (n_variance_components + 1) * parallel_idx) =
            collect_XXy.col(0);
      } // end loop 6

#pragma omp parallel for schedule(dynamic)
      for (int parallel_idx = 0; parallel_idx < n_gxg_idx;
           parallel_idx++) { // parallel loop 7
        // insert GXG
        MatrixXdr genotype_mask;
        int n_gxg_snps = n_gxg_snps_list[parallel_idx];
        MatrixXdr focal_snp_gtype;
        MatrixXdr gxg_allelecount_means;
        MatrixXdr gxg_allelecount_stds;
        int gxg_i = gxg_indices[parallel_idx];

        MatrixXdr gxg_mask =
            binary_gxg_mask.col(parallel_idx)
                .block(block_index * block_size, 0, block_size, 1);
        int gxg_snps_in_block = gxg_mask.sum();
        n_mask_gxg[parallel_idx] = n_gxg_snps;
        vector<int> n_snps_variance_component = {n_snps, n_gxg_snps};

        focal_snp_gtype.resize(n_samples, 1);
        global_snp_index = -1;
        read_focal_snp(bed_file, focal_snp_gtype, gxg_i, n_samples, n_snps,
                       global_snp_index);
        normalize_genotype(focal_snp_gtype, n_samples);

        int focal_snp_block = std::min(gxg_i / step_size, n_blocks - 1);
        int focal_snp_local_index = gxg_i - (step_size * focal_snp_block);

        MatrixXdr temp_gxg;

        compute_block_stats(gxg_genotype_blocks[parallel_idx],
                            gxg_allelecount_means, gxg_allelecount_stds,
                            n_samples, gxg_snps_in_block);

        MatrixXdr temp_XXy =
            MatrixXdr::Zero(n_samples, n_variance_components + 1);
        temp_XXy.col(0) = collect_XXy.col(0);
        temp_XXy.col(1) = collect_Gy.col(parallel_idx);
        temp_XXy.col(2) = pheno;

        MatrixXdr scaled_vec;
        for (int i = 0; i < (n_variance_components + 1); i++) {

          scaled_vec = temp_XXy.col(i).array() * focal_snp_gtype.col(0).array();
          MatrixXdr temp_GUy = compute_XXz(
              scaled_vec, gxg_allelecount_means, gxg_allelecount_stds,
              pheno_mask, 1, n_samples, gxg_genotype_blocks[parallel_idx],
              focal_snp_local_index, false);
          temp_GUy = temp_GUy.array() * focal_snp_gtype.col(0).array();
          collect_XXUy.col(((n_variance_components + 1)) + i +
                           (n_variance_components + 1) *
                               (n_variance_components + 1) * parallel_idx) +=
              temp_GUy / n_snps_variance_component[1];
        }

        gxg_genotype_blocks[parallel_idx].clear_block();
      } // end of parallel loop 7
    }
  }

#pragma omp parallel for schedule(dynamic)
  for (int parallel_idx = 0; parallel_idx < n_gxg_idx;
       parallel_idx++) { // parallel loop 8
    collect_XXUy.col((n_variance_components * (n_variance_components + 1)) + 0 +
                     (n_variance_components + 1) * (n_variance_components + 1) *
                         parallel_idx) = collect_XXy.col(0);
    collect_XXUy.col((n_variance_components * (n_variance_components + 1)) + 1 +
                     (n_variance_components + 1) * (n_variance_components + 1) *
                         parallel_idx) = collect_Gy.col(parallel_idx);
    collect_XXUy.col((n_variance_components * (n_variance_components + 1)) + 2 +
                     (n_variance_components + 1) * (n_variance_components + 1) *
                         parallel_idx) = pheno;
  } // end of parallel loop 8

#pragma omp parallel for schedule(dynamic)
  for (int parallel_idx = 0; parallel_idx < n_gxg_idx; parallel_idx++) {
    // parallel loop 9
    vector<int> n_snps_variance_component = {n_snps, n_mask_gxg[parallel_idx]};

    int n_samples_mask = pheno_mask.sum();

    MatrixXdr GxG =
        GxGz.block(0, parallel_idx * n_randvecs, n_samples, n_randvecs);
    MatrixXdr XXUy = collect_XXUy.block(
        0,
        +(n_variance_components + 1) * (n_variance_components + 1) *
            parallel_idx,
        n_samples, +(n_variance_components + 1) * (n_variance_components + 1));

    // naming from eq. (5) in mvMAPIT paper
    MatrixXdr S(n_variance_components + 1, n_variance_components + 1);
    MatrixXdr q(n_variance_components + 1, 1);
    compute_mom_components(n_randvecs, n_variance_components, pheno,
                           random_vectors, XXz, GxG, yXXy(0, 0),
                           yGxGy(parallel_idx, 0), n_snps_variance_component,
                           n_samples_mask, S, q);

    MatrixXdr point_est = S.colPivHouseholderQr().solve(q);
    MatrixXdr cov_q;
    compute_covariance_q(n_variance_components, XXUy, point_est, cov_q);

    MatrixXdr invS = S.inverse();

    MatrixXdr cov_sigma = invS * cov_q * invS;

    VC.row(parallel_idx) = point_est.col(0).transpose();
    SE.row(parallel_idx) = cov_sigma.diagonal().array().sqrt().transpose();
  } // end of parallel loop 9

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  return Rcpp::List::create(Rcpp::Named("vc_estimate") = VC,
                            Rcpp::Named("vc_se") = SE,
                            Rcpp::Named("duration") = elapsed.count());
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
