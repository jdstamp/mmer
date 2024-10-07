/*  mmer: An R Package implementation of the Multimodal Marginal Epistasis test
 * with Rcpp Copyright (C) 2024  Julian Stamp This code is licensed under MIT
 * license (see LICENSE.md for details)
 */

// [[Rcpp::plugins(openmp)]]
#include <algorithm>
#include <omp.h>
#include <random>

#include "computation.h"
#include "compute_covariance_q.h"
#include "compute_mom_components.h"
#include "count_data.h"
#include "initialize_random_vectors.h"
#include "mme.h"
#include "read_genotype_mask.h"
#include "read_genotypes.h"
#include "read_phenotypes.h"

// [[Rcpp::export]]
Rcpp::List mme_cpp(std::string plink_file, std::string pheno_file,
                   std::string genotype_mask_file, int n_randvecs, int n_blocks,
                   int rand_seed, std::vector<int> gxg_indices, int n_threads,
                   std::string gxg_h5_dataset, std::string ld_h5_dataset) {

  auto start = std::chrono::high_resolution_clock::now();

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

  int n_variance_components = 2; // For now limited to GRM and GXG
  // initialize object to collect point_est and cov_sigma
  MatrixXdr VC;
  MatrixXdr SE;

  std::string bim_file = plink_file + ".bim";
  std::string bed_file = plink_file + ".bed";
  int n_gxg_idx = gxg_indices.size();

  MatrixXdr pheno_mask;
  MatrixXdr pheno;
  MatrixXdr focal_snps_matrix;
  MatrixXdr snp_matrix;
  MatrixXdr XXz;
  MatrixXdr GxGz;
  MatrixXdr yXXy;
  MatrixXdr yGxGy;
  MatrixXdr collect_XXy;
  MatrixXdr collect_Gy;
  MatrixXdr collect_XXUy;

  MatrixXdr random_vectors;

  genotype grm_genotype_block;

  int n_snps = count_snps_bim(bim_file);
  int n_samples = count_samples(pheno_file);

  int step_size = n_snps / n_blocks;
  int step_size_remainder = n_snps % n_blocks;
  std::vector<int> block_sizes(n_blocks, step_size);
  block_sizes.back() += step_size_remainder; // add remainder to last block

  VC.resize(n_gxg_idx, n_variance_components + 1);
  SE.resize(n_gxg_idx, n_variance_components + 1);

  read_phenotypes(n_samples, pheno_file, pheno, pheno_mask);

  double y_sum = 0;
  double y_mean = 0;
  y_sum = pheno.sum();
  y_mean = y_sum / pheno_mask.sum();
  // TODO: should we scale to 1 as well?
  for (int i = 0; i < n_samples; i++) {
    if (pheno(i, 0) != 0) {
      pheno(i, 0) = pheno(i, 0) - y_mean;
    } // center phenotype
  }
  pheno = pheno.col(0);

  focal_snps_matrix = MatrixXdr::Zero(n_samples, n_gxg_idx);
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

  XXz = MatrixXdr::Zero(n_samples, n_randvecs);
  yXXy = MatrixXdr::Zero(1, 1);
  collect_XXy = MatrixXdr::Zero(n_samples, 1);
  MatrixXdr temp_grm;

  std::vector<int> n_gxg_snps_list(n_gxg_idx);
  MatrixXdr binary_gxg_mask = MatrixXdr::Ones(n_snps, n_gxg_idx);

  ifstream bed_ifs(bed_file.c_str(), ios::in | ios::binary);
  int global_snp_index = -1;
  ifstream focal_bed_ifs(bed_file.c_str(), ios::in | ios::binary);
  int focal_global_snp_index = -1;

  // read focal snps
  std::unordered_map<int, int>
      index_to_column; // Maps line number to position in selected_lines
  for (int i = 0; i < gxg_indices.size(); ++i) {
    index_to_column[gxg_indices[i]] = i;
    // print mapping
    std::cout << "Mapping: " << gxg_indices[i] << " -> " << i << std::endl;
  }

  snp_matrix = MatrixXdr::Zero(n_samples, 1);
  // Read bed file snp by snp
  while (focal_global_snp_index < n_snps) {
    if (index_to_column.count(focal_global_snp_index + 1)) { // + 1 to know
      // if the next snp is in the list
      read_snp(focal_bed_ifs, focal_global_snp_index, snp_matrix);
      normalize_genotype(snp_matrix, n_samples);
      int column = index_to_column[focal_global_snp_index];
      focal_snps_matrix.col(column) = snp_matrix;
      std::cout << "Reading focal snp " << focal_global_snp_index << " -> " << "Column index " << column << std::endl;
    } else {
      skip_snp(focal_bed_ifs, focal_global_snp_index, n_samples);
    }
  }

  //    bed_ifs.seekg(0, std::ios::beg); // reset file pointer to beginning
  //    global_snp_index = -1;

  for (int block_index = 0; block_index < n_blocks; block_index++) {
    Rcpp::checkUserInterrupt();
    std::vector<genotype> gxg_genotype_blocks(n_gxg_idx);
    int block_size = block_sizes[block_index];
    snp_matrix = MatrixXdr::Zero(n_samples, 1);
    grm_genotype_block.set_block_parameters(n_samples, block_size);

    for (int parallel_idx = 0; parallel_idx < n_gxg_idx; parallel_idx++) {
      int n_gxg_snps;
      MatrixXdr genotype_mask = MatrixXdr::Zero(n_snps, 1);
      int gxg_i = gxg_indices[parallel_idx];
      read_genotype_mask(genotype_mask_file, n_snps, gxg_i, gxg_h5_dataset,
                         ld_h5_dataset, genotype_mask, n_gxg_snps);
      genotype_mask(gxg_i, 0) = 0;
      MatrixXdr gxg_mask =
          genotype_mask.block(block_index * step_size, 0, block_size, 1);
      int gxg_snps_in_block = gxg_mask.sum();
      gxg_genotype_blocks[parallel_idx].set_block_parameters(n_samples,
                                                             gxg_snps_in_block);
      binary_gxg_mask.col(parallel_idx) = genotype_mask;
      n_gxg_snps_list[parallel_idx] = n_gxg_snps;
    }

    auto start_encoding = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < block_size; i++) {
      read_snp(bed_ifs, global_snp_index, snp_matrix);
      grm_genotype_block.encode_snp(snp_matrix);

#pragma omp parallel for schedule(dynamic)
      for (int parallel_idx = 0; parallel_idx < n_gxg_idx;
           parallel_idx++) { // parallel loop 1
        MatrixXdr gxg_mask =
            binary_gxg_mask.col(parallel_idx)
                .block(block_index * step_size, 0, block_size, 1);
        if (gxg_mask(i, 0) == 1) {
          gxg_genotype_blocks[parallel_idx].encode_snp(snp_matrix);
        }
      } // end of parallel loop 1
    }
    auto end_encoding = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_encoding =
        end_encoding - start_encoding;
    // print duration of encoding in seconds
    std::cout << "Duration of encoding: " << elapsed_encoding.count()
              << " seconds" << std::endl;

    auto start_parallel = std::chrono::high_resolution_clock::now();
    if (block_size != 0) {
#pragma omp parallel for schedule(dynamic)
      for (int parallel_idx = -1; parallel_idx < n_gxg_idx; parallel_idx++) {
        // parallel loop 2
        if (parallel_idx < 0) {
          auto start_grm = std::chrono::high_resolution_clock::now();
          grm_genotype_block.compute_block_stats();

          temp_grm = compute_XXz(random_vectors, pheno_mask, n_randvecs,
                                 grm_genotype_block);

          for (int z_index = 0; z_index < n_randvecs; z_index++) {
            XXz.col(z_index) += temp_grm.col(z_index);
          }

          yXXy(0, 0) += compute_yXXy(grm_genotype_block, pheno);
          collect_XXy.col(0) +=
              compute_XXz(pheno, pheno_mask, 1, grm_genotype_block);
          grm_genotype_block.clear_block();
          auto end_grm = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double> elapsed_grm = end_grm - start_grm;
          // print duration of encoding in seconds
          std::cout << "Duration of GRM part: " << elapsed_grm.count()
                    << "seconds" << std::endl;
        } else {
          if (gxg_genotype_blocks[parallel_idx].n_encoded == 0) {
            continue;
          }
          auto start_gxg = std::chrono::high_resolution_clock::now();

          MatrixXdr focal_snp_gtype = focal_snps_matrix.col(parallel_idx);
          int n_gxg_snps = n_gxg_snps_list[parallel_idx];

          vector<int> n_snps_variance_component = {n_snps, n_gxg_snps};

          MatrixXdr gxg_random_vectors =
              random_vectors.array().colwise() * focal_snp_gtype.col(0).array();
          MatrixXdr temp_gxg;
          gxg_genotype_blocks[parallel_idx].compute_block_stats();

          temp_gxg = compute_XXz(gxg_random_vectors, pheno_mask, n_randvecs,
                                 gxg_genotype_blocks[parallel_idx]);
          temp_gxg =
              temp_gxg.array().colwise() * focal_snp_gtype.col(0).array();

          for (int z_index = 0; z_index < n_randvecs; z_index++) {
            GxGz.col(z_index + parallel_idx * n_randvecs) +=
                temp_gxg.col(z_index);
          }

          MatrixXdr gxg_pheno;
          gxg_pheno = pheno.array() * focal_snp_gtype.col(0).array();
          MatrixXdr temp_Gy = compute_XXz(gxg_pheno, pheno_mask, 1,
                                          gxg_genotype_blocks[parallel_idx]);
          temp_Gy = temp_Gy.array() * focal_snp_gtype.col(0).array();
          collect_Gy.col(parallel_idx) += temp_Gy;

          yGxGy(parallel_idx, 0) +=
              compute_yXXy(gxg_genotype_blocks[parallel_idx], gxg_pheno);
          //        gxg_genotype_blocks[parallel_idx].clear_block();
          auto end_gxg = std::chrono::high_resolution_clock::now();
          std::chrono::duration<double> elapsed_gxg = end_gxg - start_gxg;
          // print duration of encoding in seconds
          std::cout << "Duration of GxG part: " << elapsed_gxg.count()
                    << "seconds" << std::endl;
        }
      } // end of parallel loop 2
    }
    auto end_parallel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_parallel =
        end_parallel - start_parallel;
    // print duration of encoding in seconds
    std::cout << "Duration of parallel loop 2: " << elapsed_parallel.count()
              << " seconds" << std::endl;
  }

  collect_XXy = collect_XXy / n_snps;

  // TODO: enable this when loop is in block
  for (int i = 0; i < n_gxg_idx; i++) {
    collect_Gy.col(i) = collect_Gy.col(i) / n_gxg_snps_list[i];
  }

  bed_ifs.seekg(0, std::ios::beg); // reset file pointer to beginning
  global_snp_index = -1;

  for (int block_index = 0; block_index < n_blocks; block_index++) {
    auto start_block2 = std::chrono::high_resolution_clock::now();
    Rcpp::checkUserInterrupt();
    std::vector<genotype> gxg_genotype_blocks(n_gxg_idx);
    int block_size = block_sizes[block_index];
    MatrixXdr snp_matrix = MatrixXdr::Zero(n_samples, 1);
    grm_genotype_block.set_block_parameters(n_samples, block_size);

#pragma omp parallel for schedule(dynamic)
    for (int parallel_idx = 0; parallel_idx < n_gxg_idx; parallel_idx++) {
      // parallel loop 3
      MatrixXdr gxg_mask =
          binary_gxg_mask.col(parallel_idx)
              .block(block_index * step_size, 0, block_size, 1);
      int gxg_snps_in_block = gxg_mask.sum();
      gxg_genotype_blocks[parallel_idx].set_block_parameters(n_samples,
                                                             gxg_snps_in_block);
    } // end loop 3

    for (int i = 0; i < block_size; i++) {
      read_snp(bed_ifs, global_snp_index, snp_matrix);
      grm_genotype_block.encode_snp(snp_matrix);

#pragma omp parallel for schedule(dynamic)
      for (int parallel_idx = 0; parallel_idx < n_gxg_idx;
           parallel_idx++) { // parallel loop 4
        // initialize gxg_mask
        MatrixXdr gxg_mask =
            binary_gxg_mask.col(parallel_idx)
                .block(block_index * step_size, 0, block_size, 1);
        if (gxg_mask(i, 0) == 1) {
          gxg_genotype_blocks[parallel_idx].encode_snp(snp_matrix);
        }
      } // end of parallel loop 4
    }

    if (block_size != 0) {
      grm_genotype_block.compute_block_stats();

      MatrixXdr temp_XXXXy =
          compute_XXz(collect_XXy.col(0), pheno_mask, 1, grm_genotype_block);

#pragma omp parallel for schedule(dynamic)
      for (int parallel_idx = 0; parallel_idx < n_gxg_idx;
           parallel_idx++) { // parallel loop 5
        MatrixXdr temp_XXUy = compute_XXz(collect_Gy.col(parallel_idx),
                                          pheno_mask, 1, grm_genotype_block);

        collect_XXUy.col(0 + (n_variance_components + 1) *
                                 (n_variance_components + 1) * parallel_idx) +=
            temp_XXXXy / n_snps;
        collect_XXUy.col(1 + (n_variance_components + 1) *
                                 (n_variance_components + 1) * parallel_idx) +=
            temp_XXUy / n_snps;

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

        if (gxg_genotype_blocks[parallel_idx].n_encoded == 0) {
          continue;
        }

        int n_gxg_snps = n_gxg_snps_list[parallel_idx];
        MatrixXdr focal_snp_gtype = focal_snps_matrix.col(parallel_idx);

        gxg_genotype_blocks[parallel_idx].compute_block_stats();

        MatrixXdr temp_XXy =
            MatrixXdr::Zero(n_samples, n_variance_components + 1);
        temp_XXy.col(0) = collect_XXy.col(0);
        temp_XXy.col(1) = collect_Gy.col(parallel_idx);
        temp_XXy.col(2) = pheno;

        MatrixXdr scaled_vec;
        for (int i = 0; i < (n_variance_components + 1); i++) {

          scaled_vec = temp_XXy.col(i).array() * focal_snp_gtype.col(0).array();
          MatrixXdr temp_GUy = compute_XXz(scaled_vec, pheno_mask, 1,
                                           gxg_genotype_blocks[parallel_idx]);
          temp_GUy = temp_GUy.array() * focal_snp_gtype.col(0).array();
          collect_XXUy.col(((n_variance_components + 1)) + i +
                           (n_variance_components + 1) *
                               (n_variance_components + 1) * parallel_idx) +=
              temp_GUy / n_gxg_snps;
        }
      } // end of parallel loop 7
    }

    auto end_block2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_block2 = end_block2 - start_block2;
    // print duration of encoding in seconds
    std::cout << "Duration of block " << block_index + 1
              << elapsed_block2.count() << "seconds" << std::endl;
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
    vector<int> n_snps_variance_component = {n_snps,
                                             n_gxg_snps_list[parallel_idx]};

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
