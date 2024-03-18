/*  famer: An R Package implementation of the Fast Marginal Epistasis test
 *  Copyright (C) 2024  Julian Stamp
 *  This code is licensed under MIT license (see LICENSE.md for details)
 */

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

#include "computation.h"
#include "count_samples.h"
#include "count_snps_bim.h"
#include "error_handling.h"
#include "fame.h"
#include "genotype.h"
#include "read_annotation_file.h"
#include "read_covariates.h"
#include "read_genotypes.h"
#include "read_phenotypes.h"
#include <iostream>

// [[Rcpp::export]]
Rcpp::List fame_cpp(std::string plink_file, std::string pheno_file,
                    std::string annotation_file, std::string covariate_file,
                    int gxgbin, int num_evec, int snp_index, int jack_number) {
  // Initialize all variables like in FAME
  int Nenv = 1;
  MatrixXdr Enviro;
  int blocksize;
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
  MatrixXdr geno_matrix; //(p,n)

  double y_sum;
  double y_mean;

  MatrixXdr c; //(p,k)
  MatrixXdr x; //(k,n)
  MatrixXdr v; //(p,k)
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

  std::stringstream bim_file;
  std::stringstream fam_file;
  std::stringstream bed_file;
  bim_file << plink_file << ".bim";
  fam_file << plink_file << ".fam";
  bed_file << plink_file << ".bed";

  int n_snps = count_snps_bim(bim_file.str());
  annotationStruct annotation =
      read_annotation_file(annotation_file, n_snps, jack_number);

  int n_samples = count_samples(pheno_file);
  std::cout << "Number of samples: " << n_samples << std::endl;
  int fam_lines = count_fam(fam_file.str());
  if (fam_lines != n_samples) {
    exitWithError("# samples in fam file and pheno file does not match ");
  }
  std::cout << "Number of fam lines: " << fam_lines << std::endl;
  pheno = read_phenotypes(n_samples, pheno_file, mask);

  y_sum = pheno.sum();
  std::cout << "Phenotype sum: " << y_sum << std::endl;
  Enviro.resize(n_samples, Nenv);
  cout << "selected snp index: " << snp_index << endl;
  int global_snp_index = -1;
  bool missing = false;
  read_genotypes(bed_file.str(), Enviro, missing, snp_index, n_samples, n_snps,
                 global_snp_index);
  double mean_sel_snp = Enviro.array().sum() / n_samples;
  double sd_sel_snp = sqrt((mean_sel_snp * (1 - (0.5 * mean_sel_snp))));
  cout << "mean selected snp " << mean_sel_snp << endl;
  cout << "sd seletected snp " << sd_sel_snp << endl;
  Enviro.array() = Enviro.array() - mean_sel_snp;
  Enviro.array() = Enviro.array() / sd_sel_snp;

  /////GxG : to remove selected snp from X
  bool remove_self_inter = true;
  int sel_snp_jack;
  int sel_snp_local_index;
  if (remove_self_inter == true) {
    // step_size=Nsnp/jack_number;
    sel_snp_jack =
        (snp_index - 1) /
        annotation.step_size; // Boyang: step_size is the number of snps per
                              // block. step_size * block number = all snps
    if (sel_snp_jack >= jack_number) {
      std::cout << "sel_snp_jack out of bound: " << sel_snp_jack << std::endl;
      sel_snp_jack = jack_number - 1;
    }

    sel_snp_local_index =
        snp_index - (annotation.step_size *
                     sel_snp_jack); // Boyang: get the local index !!!
    std::cout << "selected snp is in " << sel_snp_jack << " jackknife block"
              << std::endl;
    std::cout << "selected snp is " << sel_snp_local_index << " th snp of block"
              << std::endl;
  }

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
    std::cout << "Total number of fixed effects: " << cov_num << std::endl;
    if (snp_fix_ef == true)
      covariate.col(cov_num - 1) = Enviro.col(0);
  } else if (covariate_file == "") {
    std::cout << "No Covariate File Specified" << std::endl;
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
  all_zb = MatrixXdr::Random(n_samples, num_evec);
  all_zb = all_zb * sqrt(3);

  boost::mt19937 seedr;
  seedr.seed(std::time(0));
  boost::normal_distribution<> dist(0, 1);
  boost::variate_generator<boost::mt19937 &, boost::normal_distribution<>>
      z_vec(seedr, dist);
  for (int i = 0; i < num_evec; i++)
    for (int j = 0; j < n_samples; j++)
      all_zb(j, i) = z_vec();

  for (int i = 0; i < num_evec; i++)
    for (int j = 0; j < n_samples; j++)
      all_zb(j, i) = all_zb(j, i) * mask(j, 0);

  if (both_side_cov == true) {
    all_Uzb.resize(n_samples, num_evec);
    for (int j = 0; j < num_evec; j++) {
      MatrixXdr w1 = covariate.transpose() * all_zb.col(j);
      MatrixXdr w2 = Q * w1;
      MatrixXdr w3 = covariate * w2;
      all_Uzb.col(j) = w3;
    }
  }

  // Computation?
  MatrixXdr output;
  MatrixXdr output_env;

  int nongen_Nbin =
      annotation.n_bin; // Boyang change here -- to make each linear component
                        // corresponding to a nonlinear component

  bool snp_in_annot = false;
  int selected_snp_bin;

  for (int j = 0; j < annotation.n_bin; j++) {
    if (annotation.annot_bool[snp_index - 1][j] == 1) {
      selected_snp_bin = j;
      snp_in_annot = true;
    }
  }

  for (int i = 0; i < nongen_Nbin;
       i++) { // Boyang: here extend len to deal with gxg; change Nenv to
              // nongen_Nbin
    if (i == gxgbin) {
      if (remove_self_inter == true && snp_in_annot == true &&
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
                            num_evec); // Boyang: v3 change nongen_Nbin to 1;
                                       // we only want 1 gxg component
  std::cout << "Matrix dimensions XXz: " << XXz.rows() << "x" << XXz.cols() << std::endl;
  if (both_side_cov == true) {
    UXXz = MatrixXdr::Zero(n_samples,
                           (annotation.n_bin + 1) *
                               num_evec); // Boyang: v3 change nongen_Nbin to
                                          // 1;  we only want 1 gxg component
    XXUz = MatrixXdr::Zero(n_samples,
                           (annotation.n_bin + 1) *
                               num_evec); // Boyang: v3 change nongen_Nbin to
                                          // 1;  we only want 1 gxg component
  }
  yXXy = MatrixXdr::Zero(
      annotation.n_bin + 1,
      1); // Boyang: v3 change nongen_Nbin to 1;  we only want 1 gxg component

  allgen_mail.resize(annotation.n_bin);

  int bin_index = 0;
  ///// code for handeling overlapping annotations
  string name = bed_file.str();
  cout << name << endl;
  ifstream ifs(name.c_str(), ios::in | ios::binary);
  bool read_header = true;
  global_snp_index = -1;
  //////////////////////////
  MatrixXdr vec1;
  MatrixXdr w1;
  MatrixXdr w2;
  MatrixXdr w3;
  //////////// analytic se
  int total_bin_num = annotation.n_bin + 1;
  wt = MatrixXdr::Zero(n_samples, total_bin_num + 1);
  vt = MatrixXdr::Zero(n_samples, (total_bin_num + 1) * (total_bin_num + 1));
  /////

  cout << "reading genotype matrix :" << endl;
  for (int jack_index = 0; jack_index < jack_number;
       jack_index++) { // Boyang: stream block

    cout << "reading " << jack_index << " block" << endl;

    int read_Nsnp = (jack_index < (jack_number - 1))
                        ? (annotation.step_size)
                        : (annotation.step_size + annotation.step_size_rem);

    for (int i = 0; i < annotation.n_bin; i++) {

      allgen_mail[i].segment_size_hori =
          floor(log(n_samples) / log(3)) - 2; // object of the mailman
      allgen_mail[i].Nsegments_hori =
          ceil(annotation.jack_bin[jack_index][i] * 1.0 /
               (allgen_mail[i].segment_size_hori * 1.0));
      allgen_mail[i].p.resize(allgen_mail[i].Nsegments_hori,
                              std::vector<int>(n_samples));
      allgen_mail[i].not_O_i.resize(annotation.jack_bin[jack_index][i]);
      allgen_mail[i].not_O_j.resize(n_samples);
      allgen_mail[i].index = 0;
      allgen_mail[i].Nsnp = annotation.jack_bin[jack_index][i];
      allgen_mail[i].Nindv = n_samples;

      allgen_mail[i].columnsum.resize(
          annotation.jack_bin[jack_index][i],
          1); // Boyang: calculate how many columns belongs to bin i
      for (int index_temp = 0; index_temp < annotation.jack_bin[jack_index][i];
           index_temp++) {
        allgen_mail[i].columnsum[index_temp] = 0;
      }
    }

    read_bed2(ifs, missing, read_Nsnp, allgen_mail, n_samples, n_snps, global_snp_index,
              annotation, snp_index);
    read_header = false;


    for (int bin_index = 0; bin_index < annotation.n_bin; bin_index++) {
      int num_snp;
      num_snp = allgen_mail[bin_index].index;

      cout << "bin_index: " << bin_index << " num_snps: " << num_snp << endl;

      if (num_snp != 0) {
        stds.resize(num_snp, 1);
        means.resize(num_snp, 1);
        for (int i = 0; i < num_snp; i++) {
          means(i, 0) = (double)allgen_mail[bin_index].columnsum[i] / n_samples;
          if (means(i, 0) == 2 | means(i, 0) == 0)
            cout << "mean :" << means(i, 0) << endl;
        }

        // cout << "start computing stds" << endl;
        for (int i = 0; i < num_snp; i++) {
          stds(i, 0) = 1 / sqrt((means(i, 0) * (1 - (0.5 * means(i, 0)))));
        }

        // std::cout << "means Matrix shape: (" << means.rows() << ", "
        //           << means.cols() << ")" << std::endl;
        // std::cout << "means.row(0): " << means.col(0) << std::endl;
        //
        // std::cout << "stds Matrix shape: (" << stds.rows() << ", "
        //           << stds.cols() << ")" << std::endl;
        // std::cout << "stds.row(0): " << stds.col(0) << std::endl;
        g = allgen_mail[bin_index];
        g.segment_size_hori = floor(log(n_samples) / log(3)) - 2;
        g.Nsegments_hori = ceil(annotation.jack_bin[jack_index][bin_index] *
                                1.0 / (g.segment_size_hori * 1.0));
        g.p.resize(g.Nsegments_hori, std::vector<int>(n_samples));
        g.not_O_i.resize(annotation.jack_bin[jack_index][bin_index]);
        g.not_O_j.resize(n_samples);

        int p = g.Nsnp;
        int n = g.Nindv;
        int k = num_evec;
        bool fast_mode = true;
        bool memory_efficient = false;
        bool var_normalize = false;
        c.resize(p, k);
        x.resize(k, n);
        v.resize(p, k);
        sum2.resize(p, 1);
        sum.resize(p, 1);

        if (!fast_mode && !memory_efficient) {
          geno_matrix.resize(p, n);
          g.generate_eigen_geno(geno_matrix, var_normalize);
        }

        // TODO: Initialization of c with gaussian distribution
        c = MatrixXdr::Random(p, k);

        // Initial intermediate data structures
        blocksize = k;
        hsegsize = g.segment_size_hori; // = log_3(n)
        int hsize = pow(3, hsegsize);
        int vsegsize = g.segment_size_ver; // = log_3(p)
        int vsize = pow(3, vsegsize);
        bool exclude_sel_snp;

        partialsums = new double[blocksize];
        sum_op = new double[blocksize];
        yint_e = new double[hsize * blocksize];
        yint_m = new double[hsize * blocksize];
        memset(yint_m, 0, hsize * blocksize * sizeof(double));
        memset(yint_e, 0, hsize * blocksize * sizeof(double));

        y_e = new double *[g.Nindv];
        for (int i = 0; i < g.Nindv; i++) {
          y_e[i] = new double[blocksize];
          memset(y_e[i], 0, blocksize * sizeof(double));
        }

        y_m = new double *[hsegsize];
        for (int i = 0; i < hsegsize; i++)
          y_m[i] = new double[blocksize];

        output = compute_XXz(num_snp, all_zb, means, stds, mask, num_evec,
                             n_samples, sel_snp_local_index, sum_op, g, yint_m,
                             y_m, p, yint_e, y_e, partialsums);

        std::cout << "output.row(0)" << output.row(0) << std::endl;

        if (remove_self_inter == true) {
          if (sel_snp_jack == jack_index &&
              annotation.annot_bool[snp_index - 1][bin_index] == 1) {
            exclude_sel_snp = true;
          }
        }

        MatrixXdr scaled_pheno;

        int env_index = 0;
        if (bin_index == gxgbin) {
          MatrixXdr env_all_zb =
              all_zb.array().colwise() * Enviro.col(env_index).array();
          output_env = compute_XXz(num_snp, env_all_zb, means, stds, mask,
                                   num_evec, n_samples, sel_snp_local_index,
                                   sum_op, g, yint_m, y_m, p, yint_e, y_e,
                                   partialsums); // z is the random vector

          output_env =
              output_env.array().colwise() * Enviro.col(env_index).array();
          for (int z_index = 0; z_index < num_evec; z_index++) {
            XXz.col(((annotation.n_bin) * num_evec) + z_index) +=
                output_env.col(z_index); // Boyang: change env_index to
                                         // bin_index; v3: change bin_index to 1
          }
          // cout << "finish XXz compuation" << endl;
        }

        ///
        if (bin_index == gxgbin) { // Boyang: v3: add condition check
          //    cout << "start scale_pheno compuation" << endl;
          scaled_pheno = pheno.array() * Enviro.col(env_index).array();
          // cout << "finish sacled pheno computation" << endl;
          MatrixXdr temp = compute_XXy(num_snp, scaled_pheno, means, stds, mask,
                                       sel_snp_local_index, n_samples, sum_op,
                                       g, yint_m, y_m, p, yint_e, y_e,
                                       partialsums); // Here is the memory error
          // cout << "finish computing XXy" << endl;
          temp = temp.array() * Enviro.col(env_index).array();
          wt.col(annotation.n_bin) +=
              temp; // Boyang: v3 change env_index to bin_index; v3:
                    // change bin_index to 0
          if (both_side_cov == false)
            yXXy(annotation.n_bin, 0) += compute_yXXy(
                num_snp, scaled_pheno, means, stds, sel_snp_local_index, sum_op,
                g, yint_m, y_m, p, partialsums);
        }
        exclude_sel_snp = false; // v5: move exclude_sel_snp outside of the
                                 // conditional block above

        ////// wt
        wt.col(bin_index) += compute_XXy(
            num_snp, pheno, means, stds, mask, sel_snp_local_index, n_samples,
            sum_op, g, yint_m, y_m, p, yint_e, y_e, partialsums);

        for (int z_index = 0; z_index < num_evec; z_index++) {
          XXz.col((bin_index * num_evec) + z_index) +=
              output.col(z_index); /// save whole sample

          if (both_side_cov == true) {
            vec1 = output.col(z_index);
            w1 = covariate.transpose() * vec1;
            w2 = Q * w1;
            w3 = covariate * w2;
            UXXz.col((bin_index * num_evec) + z_index) += w3;
          }
        }

        if (both_side_cov == true) {
          output =
              compute_XXUz(num_snp, num_evec, n_samples, means, stds, mask,
                           sum_op, g, yint_m, y_m, p, yint_e, y_e, partialsums);

          for (int z_index = 0; z_index < num_evec; z_index++) {
            XXUz.col((bin_index * num_evec) + z_index) +=
                output.col(z_index); /// save whole sample
          }
        }

        if (both_side_cov == false)
          yXXy(bin_index, 0) +=
              compute_yXXy(num_snp, pheno, means, stds, sel_snp_local_index,
                           sum_op, g, yint_m, y_m, p, partialsums);
        else
          yXXy(bin_index, 0) +=
              compute_yVXXVy(num_snp, pheno, means, stds, num_evec, sum_op, g,
                             yint_m, y_m, p, partialsums);
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

  cout << " Reading and computing  of all blocks are finished" << endl;

  //

  /// analytic se
  for (int i = 0; i < total_bin_num; i++) {
    wt.col(i) = wt.col(i) / annotation.len[i];
    cout << "bin " << i << " length: " << annotation.len[i] << endl;
    double dx = (wt.col(i).array() * pheno.array()).sum();
    cout << "bin " << i << " " << dx << endl; // cout<<yXXy(i,0)<<endl;
  }
  wt.col(total_bin_num) = pheno;

  ifstream ifs_2(name.c_str(), ios::in | ios::binary);
  read_header = true;
  global_snp_index = -1;

  for (int jack_index = 0; jack_index < jack_number; jack_index++) {

    int read_Nsnp = (jack_index < (jack_number - 1))
                        ? (annotation.step_size)
                        : (annotation.step_size + annotation.step_size_rem);

    for (int i = 0; i < annotation.n_bin; i++) {
      allgen_mail[i].segment_size_hori = floor(log(n_samples) / log(3)) - 2;
      allgen_mail[i].Nsegments_hori =
          ceil(annotation.jack_bin[jack_index][i] * 1.0 /
               (allgen_mail[i].segment_size_hori * 1.0));
      allgen_mail[i].p.resize(allgen_mail[i].Nsegments_hori,
                              std::vector<int>(n_samples));
      allgen_mail[i].not_O_i.resize(annotation.jack_bin[jack_index][i]);
      allgen_mail[i].not_O_j.resize(n_samples);
      allgen_mail[i].index = 0;
      allgen_mail[i].Nsnp = annotation.jack_bin[jack_index][i];
      allgen_mail[i].Nindv = n_samples;

      allgen_mail[i].columnsum.resize(annotation.jack_bin[jack_index][i], 1);
      for (int index_temp = 0; index_temp < annotation.jack_bin[jack_index][i];
           index_temp++)
        allgen_mail[i].columnsum[index_temp] = 0;
    }

    // bool use_1col_annot = false;
    // if (use_1col_annot == true)
    //   read_bed_1colannot(ifs_2, missing, read_Nsnp);
    // else
    read_bed2(ifs_2, missing, read_Nsnp, allgen_mail, n_samples, n_snps,
              global_snp_index, annotation, snp_index);
    read_header = false;
    MatrixXdr means; //(p,1)
    MatrixXdr stds;  //(p,1)
    genotype g;
    for (int bin_index = 0; bin_index < annotation.n_bin; bin_index++) {
      int num_snp = allgen_mail[bin_index].index;

      if (num_snp != 0) {
        stds.resize(num_snp, 1);
        means.resize(num_snp, 1);

        for (int i = 0; i < num_snp; i++)
          means(i, 0) = (double)allgen_mail[bin_index].columnsum[i] / n_samples;

        for (int i = 0; i < num_snp; i++)
          stds(i, 0) = 1 / sqrt((means(i, 0) * (1 - (0.5 * means(i, 0)))));
        g = allgen_mail[bin_index];
        g.segment_size_hori = floor(log(n_samples) / log(3)) - 2;
        g.Nsegments_hori = ceil(annotation.jack_bin[jack_index][bin_index] *
                                1.0 / (g.segment_size_hori * 1.0));
        g.p.resize(g.Nsegments_hori, std::vector<int>(n_samples));
        g.not_O_i.resize(annotation.jack_bin[jack_index][bin_index]);
        g.not_O_j.resize(n_samples);

        /// duplicate
        int p = g.Nsnp;
        int n = g.Nindv;
        MatrixXdr c; //(p,k)
        MatrixXdr x; //(k,n)
        MatrixXdr v; //(p,k)
        MatrixXdr sum2;
        MatrixXdr sum;
        MatrixXdr geno_matrix; //(p,n)
        int k = num_evec;
        bool fast_mode = true;
        bool memory_efficient = false;
        bool var_normalize = false;
        c.resize(p, k);
        x.resize(k, n);
        v.resize(p, k);
        sum2.resize(p, 1);
        sum.resize(p, 1);

        if (!fast_mode && !memory_efficient) {
          geno_matrix.resize(p, n);
          g.generate_eigen_geno(geno_matrix, var_normalize);
        }

        // TODO: Initialization of c with gaussian distribution
        c = MatrixXdr::Random(p, k);

        // Initial intermediate data structures
        int hsize = pow(3, hsegsize);
        int vsegsize = g.segment_size_ver; // = log_3(p)
        int vsize = pow(3, vsegsize);
        bool exclude_sel_snp;

        partialsums = new double[blocksize];
        sum_op = new double[blocksize];
        yint_e = new double[hsize * blocksize];
        yint_m = new double[hsize * blocksize];
        memset(yint_m, 0, hsize * blocksize * sizeof(double));
        memset(yint_e, 0, hsize * blocksize * sizeof(double));

        y_e = new double *[g.Nindv];
        for (int i = 0; i < g.Nindv; i++) {
          y_e[i] = new double[blocksize];
          memset(y_e[i], 0, blocksize * sizeof(double));
        }

        y_m = new double *[hsegsize];
        for (int i = 0; i < hsegsize; i++)
          y_m[i] = new double[blocksize];
        //// end duplicate

        // total_bin_num=annotation.n_bin+nongen_Nbin;
        MatrixXdr val_temp;
        for (int i = 0; i < (total_bin_num + 1); i++) {
          // cout<<"test "<<i<<endl;
          MatrixXdr val_temp = compute_XXy(
              num_snp, wt.col(i), means, stds, mask, sel_snp_local_index,
              n_samples, sum_op, g, yint_m, y_m, p, yint_e, y_e, partialsums);
          vt.col((bin_index * (total_bin_num + 1)) + i) +=
              val_temp / annotation.len[bin_index]; // Boyang: what is vt??
        }

        if (remove_self_inter == true) // self-interaction
          if (sel_snp_jack == jack_index &&
              annotation.annot_bool[snp_index - 1][bin_index] == 1)
            exclude_sel_snp = true;

        MatrixXdr scaled_vec;
        if (bin_index == gxgbin) {
          for (int i = 0; i < (total_bin_num + 1); i++) {
            scaled_vec = wt.col(i).array() *
                         Enviro.col(0).array(); // change env_index to 0
            MatrixXdr temp = compute_XXy(
                num_snp, scaled_vec, means, stds, mask, sel_snp_local_index,
                n_samples, sum_op, g, yint_m, y_m, p, yint_e, y_e, partialsums);
            temp =
                temp.array() * Enviro.col(0).array(); // change env_index to 0

            vt.col(((annotation.n_bin) * (total_bin_num + 1)) + i) +=
                temp /
                annotation.len[annotation.n_bin]; // Boyang: could be right; v3:
                                                  // change env_index to 0
          }
        }

        exclude_sel_snp = false;
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

  cout << "size of the bins :" << endl;
  for (int i = 0; i < annotation.n_bin + 1;
       i++) // Boyang: v3 change nongen_Nbin to 1
    cout << "bin " << i << " : " << annotation.len[i] << endl;

  annotation.n_bin = annotation.n_bin + 1; // Boyang: v3 change nongen_Nbin to 1
  //////////////////////////////////////
  ///////////////////////////////////// compute the elements of normal equation
  ///:

  /// normal equations LHS
  MatrixXdr A_trs(annotation.n_bin, annotation.n_bin);
  MatrixXdr b_trk(annotation.n_bin, 1);
  MatrixXdr c_yky(annotation.n_bin, 1);

  MatrixXdr X_l(annotation.n_bin + 1, annotation.n_bin + 1);
  MatrixXdr Y_r(annotation.n_bin + 1, 1);
  // int bin_index=0;
  int jack_index = jack_number;
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

      // cout<<"iiiiiiiiii"<<i<<endl;
      B1 = XXz.block(0, i * num_evec, n_samples, num_evec);
      B1 = all_zb.array() * B1.array();
      b_trk(i, 0) = B1.sum() / annotation.len[i] / num_evec;
    }

    c_yky(i, 0) = yXXy(i, 0) / annotation.len[i];

    if (both_side_cov == true) {
      B1 = XXz.block(0, i * num_evec, n_samples, num_evec);
      C1 = B1.array() * all_Uzb.array();
      C2 = C1.colwise().sum();
      tk_res = C2.sum();
      tk_res = tk_res / annotation.len[i] / num_evec;
      b_trk(i, 0) = Nindv_mask - tk_res;
    }
    for (int j = i; j < annotation.n_bin; j++) {
      B1 = XXz.block(0, i * num_evec, n_samples, num_evec);
      B2 = XXz.block(0, j * num_evec, n_samples, num_evec);
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

        B1 = XXUz.block(0, (i * num_evec), n_samples, num_evec);
        B2 = UXXz.block(0, (j * num_evec), n_samples, num_evec);
        C1 = B1.array() * B2.array();
        C2 = C1.colwise().sum();
        trkij_res3 = C2.sum();

        trkij += trkij_res3 - trkij_res1 - trkij_res1;
      }

      trkij = trkij / annotation.len[i] / annotation.len[j] / num_evec;
      A_trs(i, j) = trkij;
      A_trs(j, i) = trkij;
    }
  }

  ////solve normal equation as :

  X_l << A_trs, b_trk, b_trk.transpose(), NC;
  Y_r << c_yky, yy;

  MatrixXdr herit = X_l.colPivHouseholderQr().solve(Y_r);

  cout << "Xl" << X_l << endl;
  cout << "Yl" << Y_r << endl;

  for (int i = 0; i < (annotation.n_bin + 1); i++)
    point_est(i, 0) = herit(i, 0);

  cout << "sigms" << endl;
  cout << point_est << endl;

  double temp_sum = point_est.sum();
  double temp_sig = 0;
  for (int j = 0; j < annotation.n_bin; j++) {
    herit_est(j, 0) = point_est(j, 0) / temp_sum;
    temp_sig += herit_est(j, 0);
  }
  herit_est(annotation.n_bin, 0) = temp_sig;

  cout << "herit" << endl;
  cout << herit_est << endl;

  /////compute SE analytic version
  MatrixXdr e;
  e = MatrixXdr::Zero(n_samples, 1);
  for (int i = 0; i < (annotation.n_bin + 1); i++)
    e += herit(i, 0) * wt.col(i);

  double d = 0;
  for (int i = 0; i < annotation.n_bin; i++)
    d += herit(i, 0) * c_yky(i, 0);

  d += herit(annotation.n_bin, 0) * yy;
  cout << "d : " << d << endl;
  /// vt.col(t*annotation.n_bin + i) = X_tX_tw_i

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
      // double term2=(wt.col(k).array()*e.array()).sum();
      // double term3=(wt.col(l).array()*e.array()).sum();
      // double all=2*(sum1-term2-term3+d);
      double all = 2 * sum1;
      cov_q(k, l) = all;
    }
  }

  cout << cov_q << endl;

  MatrixXdr inver_X = X_l.inverse();

  // cout<<inver_X<<endl;
  MatrixXdr cov_sigma = inver_X * cov_q * inver_X;
  cout << cov_sigma << endl;

  cout << "se" << endl;
  for (int i = 0; i < (annotation.n_bin + 1); i++)
    cout << sqrt(cov_sigma(i, i)) << endl;
  std::ofstream outfile;
  string add_output = "OUTPUT_FILE_PATH";
  outfile.open(add_output.c_str(), std::ios_base::out);
  for (int j = 0; j <= annotation.n_bin; j++) {
    outfile << "sigma^2_" << j << ": " << point_est(j, 0)
            << " se: " << sqrt(cov_sigma(j, j)) << endl;
    cout << "sigma^2_" << j << ": " << point_est(j, 0)
         << " se: " << sqrt(cov_sigma(j, j)) << endl;
  }

  return (0);
}
