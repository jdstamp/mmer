/*  famer: An R Package implementation of the Fast Marginal Epistasis test
 *  Copyright (C) 2024  Julian Stamp
 *  This code is licensed under MIT license (see LICENSE.md for details)
 */

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

#include "computation.h"
#include "count_data.h"
#include "error_handling.h"
#include "fame.h"
#include "genotype.h"
#include "read_annotation_file.h"
#include "read_covariates.h"
#include "read_genotypes.h"
#include "read_phenotypes.h"
#include "set_block_parameters.h"
#include <iostream>
#include <chrono>

// [[Rcpp::export]]
Rcpp::List fame_cpp(std::string plink_file, std::string pheno_file,
                    std::string annotation_file, std::string covariate_file,
                    int n_randvecs, int focal_snp_index, int n_blocks) {

  auto start = std::chrono::high_resolution_clock::now();
  int gxgbin = 0;
  // Initialize all variables like in FAME
  MatrixXdr focal_snp_gtype;
  int hsegsize; // = log_3(n)
  double *partialsums;
  double *sum_op;
  double *yint_e;
  double *yint_m;
  double **y_e;
  double **y_m;

  MatrixXdr mask;
  MatrixXdr pheno;
  MatrixXdr covariate;
  MatrixXdr v1; // W^ty
  MatrixXdr v2; // QW^ty
  MatrixXdr v3; // WQW^ty
  MatrixXdr Q;
  MatrixXdr new_pheno;

  genotype g;
  MatrixXdr geno_matrix;      //(p,n)
  bool allow_missing = false; // bolean for reading genotype data

  double y_sum;
  double y_mean;

  MatrixXdr c;     //(p,k)
  MatrixXdr x;     //(k,n)
  MatrixXdr v;     //(p,k)
  MatrixXdr means; //(p,1)
  MatrixXdr stds;  //(p,1)
  MatrixXdr sum2;
  MatrixXdr sum;

  MatrixXdr all_zb;
  MatrixXdr all_Uzb;
  MatrixXdr res;
  MatrixXdr XXz;
  MatrixXdr Xy;
  MatrixXdr UXXz;
  MatrixXdr XXUz;
  MatrixXdr yXXy;
  vector<genotype> allgen_mail;

  MatrixXdr wt;
  MatrixXdr vt;
  //

  std::stringstream bim_stream;
  std::stringstream fam_stream;
  std::stringstream bed_stream;
  bim_stream << plink_file << ".bim";
  fam_stream << plink_file << ".fam";
  bed_stream << plink_file << ".bed";
  std::string bim_file = bim_stream.str();
  std::string fam_file = fam_stream.str();
  std::string bed_file = bed_stream.str();

  int n_snps = count_snps_bim(bim_file);
  int n_samples = count_samples(pheno_file);
  int fam_lines = count_fam(fam_file);
  annotationStruct annotation =
      read_annotation_file(annotation_file, n_snps, n_blocks);

  if (fam_lines != n_samples) {
    exitWithError("# samples in fam file and pheno file does not match ");
  }

  // read data into matrices pheno and mask
  read_phenotypes(n_samples, pheno_file, pheno, mask);

  focal_snp_gtype.resize(n_samples, 1);

  int global_snp_index = -1; // keep track of how much was read for
                             // block wise reading of genotypes
  read_focal_snp(bed_file, focal_snp_gtype, focal_snp_index, n_samples, n_snps,
                 global_snp_index);
  double mean_sel_snp = focal_snp_gtype.array().sum() / n_samples;
  double sd_sel_snp = sqrt((mean_sel_snp * (1 - (0.5 * mean_sel_snp))));
  focal_snp_gtype.array() = focal_snp_gtype.array() - mean_sel_snp;
  focal_snp_gtype.array() = focal_snp_gtype.array() / sd_sel_snp;

  // GxG : to remove selected snp from X
  int focal_snp_block = (focal_snp_index - 1) / annotation.step_size;
  // step_size is the number of snps per block. step_size * n_blocks = #snps
  if (focal_snp_block >= n_blocks) {
    // in case above math did not work - fix
    focal_snp_block = n_blocks - 1;
  }
  int sel_snp_local_index = // block_size of focal SNP inside it's block
      focal_snp_index - (annotation.step_size * focal_snp_block) - 1;

  // Covariate handling - needs cleanup
  std::string covname = "";
  bool use_cov = false;
  bool both_side_cov = false;
  bool snp_fix_ef = false;
  int cov_num;
  if (covariate_file != "") {
    use_cov = true;
    cov_num = read_covariates(false, n_samples, covariate_file, covname,
                              covariate, snp_fix_ef);
    if (snp_fix_ef == true)
      covariate.col(cov_num - 1) = focal_snp_gtype.col(0);
  } else if (covariate_file == "") {
    both_side_cov = false;
  }

  /// regress out cov from phenotypes - needs cleanup
  if (use_cov == true) {
    MatrixXdr mat_mask = mask.replicate(1, cov_num);
    covariate = covariate.cwiseProduct(mat_mask);

    MatrixXdr WtW = covariate.transpose() * covariate;
    Q = WtW.inverse(); // Q=(W^tW)^-1

    double phen_sd = 0;
    if (both_side_cov == false) {
      MatrixXdr v1 = covariate.transpose() * pheno; // W^ty
      MatrixXdr v2 = Q * v1;                        // QW^ty
      MatrixXdr v3 = covariate * v2;                // WQW^ty
      new_pheno = pheno - v3;
      pheno = new_pheno.cwiseProduct(mask);

      y_sum = pheno.sum();
      y_mean = y_sum / mask.sum();
      for (int i = 0; i < n_samples; i++) {
        phen_sd += (pheno(i, 0) - y_mean) * (pheno(i, 0) - y_mean);
        if (pheno(i, 0) != 0)
          pheno(i, 0) = pheno(i, 0) - y_mean; // center phenotype
      }
      phen_sd = sqrt(phen_sd / (mask.sum() - 1));
      pheno = pheno / phen_sd;
      y_sum = pheno.sum();
    }

    if (both_side_cov == true) {
      y_sum = pheno.sum();
      y_mean = y_sum / mask.sum();
      for (int i = 0; i < n_samples; i++) {
        phen_sd += (pheno(i, 0) - y_mean) * (pheno(i, 0) - y_mean);
        if (pheno(i, 0) != 0)
          pheno(i, 0) = pheno(i, 0) - y_mean; // center phenotype
      }
      phen_sd = sqrt(phen_sd / (mask.sum() - 1));
      pheno = pheno / phen_sd;
      y_sum = pheno.sum();

      v1 = covariate.transpose() * pheno; // W^ty
      v2 = Q * v1;                        // QW^ty
      v3 = covariate * v2;                // WQW^ty
      new_pheno = pheno - v3;
      new_pheno = new_pheno.cwiseProduct(mask);
    }
  }
  if (use_cov == false) {
    y_sum = pheno.sum();
    y_mean = y_sum / mask.sum();
    for (int i = 0; i < n_samples; i++) {
      if (pheno(i, 0) != 0)
        pheno(i, 0) = pheno(i, 0) - y_mean; // center phenotype
    }
    y_sum = pheno.sum();
  }

  // random vector stuff
  // define random vector z's
  all_zb = MatrixXdr::Random(n_samples, n_randvecs);
  all_zb = all_zb * sqrt(3);

  boost::mt19937 seedr;
  seedr.seed(std::time(0));
  seedr.seed(123); // TODO: remove seed for actual algorithm
  boost::normal_distribution<> dist(0, 1);
  boost::variate_generator<boost::mt19937 &, boost::normal_distribution<>>
      z_vec(seedr, dist);

  for (int i = 0; i < n_randvecs; i++)
    for (int j = 0; j < n_samples; j++)
      all_zb(j, i) = z_vec();

  for (int i = 0; i < n_randvecs; i++)
    for (int j = 0; j < n_samples; j++)
      all_zb(j, i) = all_zb(j, i) * mask(j, 0);

  if (both_side_cov == true) {
    all_Uzb.resize(n_samples, n_randvecs);
    for (int j = 0; j < n_randvecs; j++) {
      MatrixXdr w1 = covariate.transpose() * all_zb.col(j);
      MatrixXdr w2 = Q * w1;
      MatrixXdr w3 = covariate * w2;
      all_Uzb.col(j) = w3;
    }
  }

  // Computation?
  MatrixXdr output;
  MatrixXdr output_env;

  bool snp_in_annot = false;
  int selected_snp_bin;

  for (int j = 0; j < annotation.n_bin; j++) {
    if (annotation.annot_bool[focal_snp_index - 1][j] == 1) {
      selected_snp_bin = j;
      snp_in_annot = true;
    }
  }

  for (int i = 0; i < annotation.n_bin;
       i++) { // Boyang: here extend len to deal with gxg; change Nenv to
              // annotation.n_bin
    if (i == gxgbin) {
      if (snp_in_annot == true &&
          i == selected_snp_bin) // Boyang: double change here, only remove
                                 // length if selected snp in current set
        annotation.len.push_back(
            annotation.selected_snps_vec[i] -
            1); /// GxG ; remove selected snp; Boyang: replace
                /// selected_snps with selected_snps_vec
      else
        annotation.len.push_back(annotation.selected_snps_vec[i]);
    }
  }

  XXz = MatrixXdr::Zero(n_samples,
                        (annotation.n_bin + 1) *
                            n_randvecs); // Boyang: v3 change nongen_Nbin to 1;
                                         // we only want 1 gxg component
  if (both_side_cov == true) {
    UXXz = MatrixXdr::Zero(n_samples,
                           (annotation.n_bin + 1) *
                               n_randvecs); // Boyang: v3 change nongen_Nbin to
                                            // 1;  we only want 1 gxg component
    XXUz = MatrixXdr::Zero(n_samples,
                           (annotation.n_bin + 1) *
                               n_randvecs); // Boyang: v3 change nongen_Nbin to
                                            // 1;  we only want 1 gxg component
  }
  yXXy = MatrixXdr::Zero(
      annotation.n_bin + 1,
      1); // Boyang: v3 change nongen_Nbin to 1;  we only want 1 gxg component

  metaData metadata = set_metadata(n_samples, n_snps);
  allgen_mail.resize(annotation.n_bin);
  ifstream bed_ifs(bed_file.c_str(), ios::in | ios::binary);
  global_snp_index = -1;

  MatrixXdr vec1;
  MatrixXdr w1;
  MatrixXdr w2;
  MatrixXdr w3;
  //////////// analytic se
  int total_bin_num = annotation.n_bin + 1;
  wt = MatrixXdr::Zero(n_samples, total_bin_num + 1);
  vt = MatrixXdr::Zero(n_samples, (total_bin_num + 1) * (total_bin_num + 1));

  for (int block_index = 0; block_index < n_blocks;
       block_index++) { // Boyang: stream block

    int read_Nsnp = (block_index < (n_blocks - 1))
                        ? (annotation.step_size)
                        : (annotation.step_size + annotation.step_size_rem);
    // could this be changed to have the remainder as the last block only?

    set_block_parameters(allgen_mail, n_samples, annotation, block_index);

      read_genotype_block(bed_ifs, read_Nsnp, allgen_mail, n_samples, n_snps,
                          global_snp_index, annotation, metadata);

    for (int bin_index = 0; bin_index < annotation.n_bin; bin_index++) {
      int block_size = allgen_mail[bin_index].block_size;

      if (block_size != 0) {
        stds.resize(block_size, 1);
        means.resize(block_size, 1);
        for (int i = 0; i < block_size; i++) {
          means(i, 0) = (double)allgen_mail[bin_index].columnsum[i] / n_samples;
          if (means(i, 0) == 2 | means(i, 0) == 0)
            cout << "mean :" << means(i, 0) << endl;
        }

        for (int i = 0; i < block_size; i++) {
          stds(i, 0) = 1 / sqrt((means(i, 0) * (1 - (0.5 * means(i, 0)))));
        }

        g = allgen_mail[bin_index];
        g.segment_size_hori = floor(log(n_samples) / log(3)) - 2;
        g.Nsegments_hori = ceil(annotation.jack_bin[block_index][bin_index] *
                                1.0 / (g.segment_size_hori * 1.0));
        g.p.resize(g.Nsegments_hori, std::vector<int>(n_samples));
        g.not_O_i.resize(annotation.jack_bin[block_index][bin_index]);
        g.not_O_j.resize(n_samples);

        int p = g.Nsnp;
        int n = g.Nindv;
        c.resize(p, n_randvecs);
        x.resize(n_randvecs, n);
        v.resize(p, n_randvecs);
        sum2.resize(p, 1);
        sum.resize(p, 1);

        // TODO: Initialization of c with gaussian distribution
        c = MatrixXdr::Random(p, n_randvecs);

        // Initial intermediate data structures
        hsegsize = g.segment_size_hori; // = log_3(n)
        int hsize = pow(3, hsegsize);
        int vsegsize = g.segment_size_ver; // = log_3(p)
        int vsize = pow(3, vsegsize);

        partialsums = new double[n_randvecs];
        sum_op = new double[n_randvecs];
        yint_e = new double[hsize * n_randvecs];
        yint_m = new double[hsize * n_randvecs];
        memset(yint_m, 0, hsize * n_randvecs * sizeof(double));
        memset(yint_e, 0, hsize * n_randvecs * sizeof(double));

        y_e = new double *[g.Nindv];
        for (int i = 0; i < g.Nindv; i++) {
          y_e[i] = new double[n_randvecs];
          memset(y_e[i], 0, n_randvecs * sizeof(double));
        }

        y_m = new double *[hsegsize];
        for (int i = 0; i < hsegsize; i++)
          y_m[i] = new double[n_randvecs];

        output = compute_XXz(block_size, all_zb, means, stds, mask, n_randvecs,
                             n_samples, sel_snp_local_index, sum_op, g, yint_m,
                             y_m, p, yint_e, y_e, partialsums, false);

        MatrixXdr scaled_pheno;

        int env_index = 0;
        if (bin_index == gxgbin) {
          bool in_gxg_block = (focal_snp_block == block_index);
            MatrixXdr env_all_zb =
              all_zb.array().colwise() * focal_snp_gtype.col(env_index).array();
          output_env = compute_XXz(block_size, env_all_zb, means, stds, mask,
                                   n_randvecs, n_samples, sel_snp_local_index,
                                   sum_op, g, yint_m, y_m, p, yint_e, y_e,
                                   partialsums, in_gxg_block); // z
                                   // is the random vector

          output_env = output_env.array().colwise() *
                       focal_snp_gtype.col(env_index).array();
          for (int z_index = 0; z_index < n_randvecs; z_index++) {
            XXz.col(((annotation.n_bin) * n_randvecs) + z_index) +=
                output_env.col(z_index);
          }
        }

        ///
        if (bin_index == gxgbin) {
            bool in_gxg_block = (focal_snp_block == block_index);
          scaled_pheno = pheno.array() * focal_snp_gtype.col(env_index).array();
          MatrixXdr temp = compute_XXy(block_size, scaled_pheno, means, stds,
                                       mask, sel_snp_local_index, n_samples,
                                       sum_op, g, yint_m, y_m, p, yint_e, y_e,
                                       partialsums, in_gxg_block); // Here is the memory error
          temp = temp.array() * focal_snp_gtype.col(env_index).array();
          wt.col(annotation.n_bin) += temp;
          if (both_side_cov == false)
            yXXy(annotation.n_bin, 0) += compute_yXXy(
                    block_size, scaled_pheno, means, stds, sel_snp_local_index,
                    sum_op, g, yint_m, y_m, p, partialsums, in_gxg_block);
        }

        ////// wt
        wt.col(bin_index) += compute_XXy(
                block_size, pheno, means, stds, mask, sel_snp_local_index,
                n_samples, sum_op, g, yint_m, y_m, p, yint_e, y_e, partialsums,
                false);

        for (int z_index = 0; z_index < n_randvecs; z_index++) {
          XXz.col((bin_index * n_randvecs) + z_index) +=
              output.col(z_index); /// save whole sample

          if (both_side_cov == true) {
            vec1 = output.col(z_index);
            w1 = covariate.transpose() * vec1;
            w2 = Q * w1;
            w3 = covariate * w2;
            UXXz.col((bin_index * n_randvecs) + z_index) += w3;
          }
        }

        if (both_side_cov == true) {
          output =
              compute_XXUz(block_size, n_randvecs, n_samples, means, stds, mask,
                           sum_op, g, yint_m, y_m, p, yint_e, y_e, partialsums);

          for (int z_index = 0; z_index < n_randvecs; z_index++) {
            XXUz.col((bin_index * n_randvecs) + z_index) +=
                output.col(z_index); /// save whole sample
          }
        }

        if (both_side_cov == false)
          yXXy(bin_index, 0) +=
                  compute_yXXy(block_size, pheno, means, stds,
                               sel_snp_local_index,
                               sum_op, g, yint_m, y_m, p, partialsums, false);
        else
          yXXy(bin_index, 0) +=
              compute_yVXXVy(block_size, pheno, means, stds, n_randvecs, sum_op,
                             g, yint_m, y_m, p, partialsums);
        delete[] sum_op;
        delete[] partialsums;
        delete[] yint_e;
        delete[] yint_m;
        for (int i = 0; i < hsegsize; i++)
          delete[] y_m[i];
        delete[] y_m;
        for (int i = 0; i < g.Nindv; i++)
          delete[] y_e[i];
        delete[] y_e;
        std::vector<std::vector<int>>().swap(g.p);
        std::vector<std::vector<int>>().swap(g.not_O_j);
        std::vector<std::vector<int>>().swap(g.not_O_i);
        std::vector<std::vector<int>>().swap(allgen_mail[bin_index].p);
        std::vector<std::vector<int>>().swap(allgen_mail[bin_index].not_O_j);
        std::vector<std::vector<int>>().swap(allgen_mail[bin_index].not_O_i);
        g.columnsum.clear();
        g.columnsum2.clear();
        g.columnmeans.clear();
        g.columnmeans2.clear();
        allgen_mail[bin_index].columnsum.clear();
        allgen_mail[bin_index].columnsum2.clear();
        allgen_mail[bin_index].columnmeans.clear();
        allgen_mail[bin_index].columnmeans2.clear();
      }
    }
  }

  /// analytic se
  for (int i = 0; i < total_bin_num; i++) {
    wt.col(i) = wt.col(i) / annotation.len[i];
    double dx = (wt.col(i).array() * pheno.array()).sum();
  }
  wt.col(total_bin_num) = pheno;

  ifstream ifs_2(bed_file.c_str(), ios::in | ios::binary);
  global_snp_index = -1;

  for (int block_index = 0; block_index < n_blocks; block_index++) {

    int read_Nsnp = (block_index < (n_blocks - 1))
                        ? (annotation.step_size)
                        : (annotation.step_size + annotation.step_size_rem);

    set_block_parameters(allgen_mail, n_samples, annotation, block_index);

    read_genotype_block(ifs_2, read_Nsnp, allgen_mail, n_samples, n_snps,
                          global_snp_index, annotation, metadata);
    MatrixXdr means; //(p,1)
    MatrixXdr stds;  //(p,1)
    genotype g;
    for (int bin_index = 0; bin_index < annotation.n_bin; bin_index++) {
      int num_snp = allgen_mail[bin_index].block_size;

      if (num_snp != 0) {
        stds.resize(num_snp, 1);
        means.resize(num_snp, 1);

        for (int i = 0; i < num_snp; i++)
          means(i, 0) = (double)allgen_mail[bin_index].columnsum[i] / n_samples;

        for (int i = 0; i < num_snp; i++)
          stds(i, 0) = 1 / sqrt((means(i, 0) * (1 - (0.5 * means(i, 0)))));
        g = allgen_mail[bin_index];
        g.segment_size_hori = floor(log(n_samples) / log(3)) - 2;
        g.Nsegments_hori = ceil(annotation.jack_bin[block_index][bin_index] *
                                1.0 / (g.segment_size_hori * 1.0));
        g.p.resize(g.Nsegments_hori, std::vector<int>(n_samples));
        g.not_O_i.resize(annotation.jack_bin[block_index][bin_index]);
        g.not_O_j.resize(n_samples);

        /// duplicate
        int p = g.Nsnp;
        int n = g.Nindv;
        MatrixXdr c; //(p,n_randvecs)
        MatrixXdr x; //(n_randvecs,n)
        MatrixXdr v; //(p,n_randvecs)
        MatrixXdr sum2;
        MatrixXdr sum;
        MatrixXdr geno_matrix; //(p,n)
        c.resize(p, n_randvecs);
        x.resize(n_randvecs, n);
        v.resize(p, n_randvecs);
        sum2.resize(p, 1);
        sum.resize(p, 1);

        // TODO: Initialization of c with gaussian distribution
        c = MatrixXdr::Random(p, n_randvecs);

        // Initial intermediate data structures
        int hsize = pow(3, hsegsize);
        int vsegsize = g.segment_size_ver; // = log_3(p)
        int vsize = pow(3, vsegsize);

        partialsums = new double[n_randvecs];
        sum_op = new double[n_randvecs];
        yint_e = new double[hsize * n_randvecs];
        yint_m = new double[hsize * n_randvecs];
        memset(yint_m, 0, hsize * n_randvecs * sizeof(double));
        memset(yint_e, 0, hsize * n_randvecs * sizeof(double));

        y_e = new double *[g.Nindv];
        for (int i = 0; i < g.Nindv; i++) {
          y_e[i] = new double[n_randvecs];
          memset(y_e[i], 0, n_randvecs * sizeof(double));
        }

        y_m = new double *[hsegsize];
        for (int i = 0; i < hsegsize; i++)
          y_m[i] = new double[n_randvecs];
        //// end duplicate

        // total_bin_num=annotation.n_bin+nongen_Nbin;
        MatrixXdr val_temp;
        for (int i = 0; i < (total_bin_num + 1); i++) {
          MatrixXdr val_temp = compute_XXy(
                  num_snp, wt.col(i), means, stds, mask, sel_snp_local_index,
                  n_samples, sum_op, g, yint_m, y_m, p, yint_e, y_e,
                  partialsums, false);
          vt.col((bin_index * (total_bin_num + 1)) + i) +=
              val_temp / annotation.len[bin_index]; // Boyang: what is vt??
        }

        MatrixXdr scaled_vec;
        if (bin_index == gxgbin) {
          for (int i = 0; i < (total_bin_num + 1); i++) {
            scaled_vec =
                wt.col(i).array() *
                focal_snp_gtype.col(0).array(); // change env_index to 0
            MatrixXdr temp = compute_XXy(
                    num_snp, scaled_vec, means, stds, mask, sel_snp_local_index,
                    n_samples, sum_op, g, yint_m, y_m, p, yint_e, y_e,
                    partialsums, false);
            temp = temp.array() *
                   focal_snp_gtype.col(0).array(); // change env_index to 0

            vt.col(((annotation.n_bin) * (total_bin_num + 1)) + i) +=
                temp /
                annotation.len[annotation.n_bin]; // Boyang: could be right; v3:
                                                  // change env_index to 0
          }
        }

        delete[] sum_op;
        delete[] partialsums;
        delete[] yint_e;
        delete[] yint_m;
        for (int i = 0; i < hsegsize; i++)
          delete[] y_m[i];
        delete[] y_m;

        for (int i = 0; i < g.Nindv; i++)
          delete[] y_e[i];
        delete[] y_e;

        std::vector<std::vector<int>>().swap(g.p);
        std::vector<std::vector<int>>().swap(g.not_O_j);
        std::vector<std::vector<int>>().swap(g.not_O_i);
        std::vector<std::vector<int>>().swap(allgen_mail[bin_index].p);
        std::vector<std::vector<int>>().swap(allgen_mail[bin_index].not_O_j);
        std::vector<std::vector<int>>().swap(allgen_mail[bin_index].not_O_i);

        g.columnsum.clear();
        g.columnsum2.clear();
        g.columnmeans.clear();
        g.columnmeans2.clear();
        allgen_mail[bin_index].columnsum.clear();
        allgen_mail[bin_index].columnsum2.clear();
        allgen_mail[bin_index].columnmeans.clear();
        allgen_mail[bin_index].columnmeans2.clear();
      }
    }
  }

  for (int i = 0; i < (total_bin_num + 1); i++) {
    vt.col((total_bin_num * (total_bin_num + 1)) + i) = wt.col(i);
  }

  annotation.n_bin = annotation.n_bin + 1;
  // compute the elements of normal equation:

  /// normal equations LHS
  MatrixXdr A_trs(annotation.n_bin, annotation.n_bin);
  MatrixXdr b_trk(annotation.n_bin, 1);
  MatrixXdr c_yky(annotation.n_bin, 1);

  MatrixXdr X_l(annotation.n_bin + 1, annotation.n_bin + 1);
  MatrixXdr Y_r(annotation.n_bin + 1, 1);
  int jack_index = n_blocks;
  MatrixXdr B1;
  MatrixXdr B2;
  MatrixXdr C1;
  MatrixXdr C2;
  double trkij;
  double yy = (pheno.array() * pheno.array()).sum();

  if (both_side_cov == true) {
    MatrixXdr Wty = covariate.transpose() * pheno; // W^ty
    MatrixXdr QWty = Q * Wty;
    double temp = (Wty.array() * QWty.array()).sum();
    yy = yy - temp;
  }

  int Nindv_mask = mask.sum();
  int NC;
  if (both_side_cov == true)
    NC = Nindv_mask - cov_num;
  else
    NC = Nindv_mask;

  MatrixXdr point_est;
  MatrixXdr herit_est;

  point_est.resize(annotation.n_bin + 1, 1);
  herit_est.resize(annotation.n_bin + 1, 1);

  MatrixXdr h1;
  MatrixXdr h2;
  MatrixXdr h3;

  double trkij_res1;
  double trkij_res2;
  double trkij_res3;
  double tk_res;

  for (int i = 0; i < annotation.n_bin; i++) {

    if (both_side_cov == false)
      b_trk(i, 0) = Nindv_mask;

    if (i >= (annotation.n_bin - 1)) { // change nongen_Nbin to 1

      B1 = XXz.block(0, i * n_randvecs, n_samples, n_randvecs);
      B1 = all_zb.array() * B1.array();
      b_trk(i, 0) = B1.sum() / annotation.len[i] / n_randvecs;
    }

    c_yky(i, 0) = yXXy(i, 0) / annotation.len[i];

    if (both_side_cov == true) {
      B1 = XXz.block(0, i * n_randvecs, n_samples, n_randvecs);
      C1 = B1.array() * all_Uzb.array();
      C2 = C1.colwise().sum();
      tk_res = C2.sum();
      tk_res = tk_res / annotation.len[i] / n_randvecs;
      b_trk(i, 0) = Nindv_mask - tk_res;
    }
    for (int j = i; j < annotation.n_bin; j++) {
      B1 = XXz.block(0, i * n_randvecs, n_samples, n_randvecs);
      B2 = XXz.block(0, j * n_randvecs, n_samples, n_randvecs);
      C1 = B1.array() * B2.array();
      C2 = C1.colwise().sum();
      trkij = C2.sum();

      if (both_side_cov == true) {

        h1 = covariate.transpose() * B1;
        h2 = Q * h1;
        h3 = covariate * h2;
        C1 = h3.array() * B2.array();
        C2 = C1.colwise().sum();
        trkij_res1 = C2.sum();

        B1 = XXUz.block(0, (i * n_randvecs), n_samples, n_randvecs);
        B2 = UXXz.block(0, (j * n_randvecs), n_samples, n_randvecs);
        C1 = B1.array() * B2.array();
        C2 = C1.colwise().sum();
        trkij_res3 = C2.sum();

        trkij += trkij_res3 - trkij_res1 - trkij_res1;
      }

      trkij = trkij / annotation.len[i] / annotation.len[j] / n_randvecs;
      A_trs(i, j) = trkij;
      A_trs(j, i) = trkij;
    }
  }

  ////solve normal equation as :

  X_l << A_trs, b_trk, b_trk.transpose(), NC;
  Y_r << c_yky, yy;

  MatrixXdr herit = X_l.colPivHouseholderQr().solve(Y_r);

  for (int i = 0; i < (annotation.n_bin + 1); i++)
    point_est(i, 0) = herit(i, 0);

  double temp_sum = point_est.sum();
  double temp_sig = 0;
  for (int j = 0; j < annotation.n_bin; j++) {
    herit_est(j, 0) = point_est(j, 0) / temp_sum;
    temp_sig += herit_est(j, 0);
  }
  herit_est(annotation.n_bin, 0) = temp_sig;

  /////compute SE analytic version
  MatrixXdr e;
  e = MatrixXdr::Zero(n_samples, 1);
  for (int i = 0; i < (annotation.n_bin + 1); i++)
    e += herit(i, 0) * wt.col(i);

  double d = 0;
  for (int i = 0; i < annotation.n_bin; i++)
    d += herit(i, 0) * c_yky(i, 0);

  d += herit(annotation.n_bin, 0) * yy;

  MatrixXdr cov_q;
  cov_q.resize(annotation.n_bin + 1, annotation.n_bin + 1);
  for (int k = 0; k < (annotation.n_bin + 1); k++) {
    for (int l = 0; l < (annotation.n_bin + 1); l++) {
      double sum1 = 0;
      for (int t = 0; t < (annotation.n_bin + 1); t++) {

        MatrixXdr z1 = vt.col((t * (annotation.n_bin + 1)) + l);
        MatrixXdr z2 = wt.col(k);
        double temp = (z1.array() * z2.array()).sum();
        sum1 += herit(t, 0) * temp;
      }
      double all = 2 * sum1;
      cov_q(k, l) = all;
    }
  }

  MatrixXdr inver_X = X_l.inverse();

  MatrixXdr cov_sigma = inver_X * cov_q * inver_X;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution time of main function: " << elapsed.count() << " seconds" << std::endl;
    return Rcpp::List::create(Rcpp::Named("Est") = point_est.col(0),
                            Rcpp::Named("SE") =
                                cov_sigma.diagonal().array().sqrt());
}

