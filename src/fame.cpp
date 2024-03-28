/*  famer: An R Package implementation of the Fast Marginal Epistasis test
 *  Copyright (C) 2024  Julian Stamp
 *  This code is licensed under MIT license (see LICENSE.md for details)
 */

#include "fame.h"
#include "computation.h"
#include "count_data.h"
#include "fit_covariates.h"
#include "read_covariates.h"
#include "read_genotypes.h"
#include "read_phenotypes.h"
#include "set_block_parameters.h"

// [[Rcpp::export]]
Rcpp::List fame_cpp(std::string plink_file, std::string pheno_file,
                    std::string covariate_file, int n_randvecs,
                    int focal_snp_index, int n_blocks, int rand_seed) {

  auto start = std::chrono::high_resolution_clock::now();
  // Mailman algo variables. TODO: Give meaningful names
  MatrixXdr focal_snp_gtype;
  int hsegsize;
  double *partialsums;
  double *sum_op;
  double *yint_e;
  double *yint_m;
  double **y_e;
  double **y_m;

  MatrixXdr pheno_mask;
  MatrixXdr pheno;

  genotype genotype_block;

  MatrixXdr allelecount_means;
  MatrixXdr allelecount_stds;

  MatrixXdr random_vectors;
  MatrixXdr XXz;
  MatrixXdr Xy;
  MatrixXdr UXXz;
  MatrixXdr XXUz;
  MatrixXdr yXXy;

  MatrixXdr wt;
  MatrixXdr vt;

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

  int step_size = n_snps / n_blocks;
  int step_size_rem = n_snps % n_blocks;
  std::vector<int> block_sizes(n_blocks, step_size);
  block_sizes.back() += step_size_rem; // add remainder to last block

  vector<int> n_snps_variance_component = {n_snps, n_snps - 1};

  if (fam_lines != n_samples) {
    throw std::runtime_error(
        "Number of samples in fam file and pheno file do not match");
  }

  // read data into matrices pheno and pheno_mask
  read_phenotypes(n_samples, pheno_file, pheno, pheno_mask);

  focal_snp_gtype.resize(n_samples, 1);

  // keep track of how much was read for block wise reading of genotypes
  int global_snp_index = -1;
  read_focal_snp(bed_file, focal_snp_gtype, focal_snp_index, n_samples, n_snps,
                 global_snp_index);

  int focal_snp_block =
      std::min((focal_snp_index - 1) / step_size, n_blocks - 1);

  int focal_snp_local_index =
      focal_snp_index - (step_size * focal_snp_block) - 1;

  // Covariate handling - needs cleanup
  double y_sum;
  double y_mean;
  if (covariate_file != "") {
    bool snp_fix_ef = false;
    MatrixXdr covariate;
    int cov_num;
    cov_num = read_covariates(false, n_samples, covariate_file, covariate,
                              snp_fix_ef);
    if (snp_fix_ef == true)
      covariate.col(cov_num - 1) = focal_snp_gtype.col(0);

    pheno = fit_covariates(pheno_mask, pheno, n_samples, y_sum, y_mean,
                           covariate, cov_num);

  } else {
    y_sum = pheno.sum();
    y_mean = y_sum / pheno_mask.sum();
    for (int i = 0; i < n_samples; i++) {
      if (pheno(i, 0) != 0)
        pheno(i, 0) = pheno(i, 0) - y_mean; // center phenotype
    }
  }

  // random vector stuff
  random_vectors = MatrixXdr::Random(n_samples, n_randvecs);

  boost::mt19937 seedr;
  seedr.seed(rand_seed < 0 ? std::time(0) : rand_seed);
  boost::normal_distribution<> dist(0, 1);
  boost::variate_generator<boost::mt19937 &, boost::normal_distribution<>>
      z_vec(seedr, dist);

  for (int i = 0; i < n_randvecs; i++)
    for (int j = 0; j < n_samples; j++)
      random_vectors(j, i) = z_vec();

  for (int i = 0; i < n_randvecs; i++)
    for (int j = 0; j < n_samples; j++)
      random_vectors(j, i) = random_vectors(j, i) * pheno_mask(j, 0);

  // Computation?
  MatrixXdr output;
  MatrixXdr output_env;

  XXz = MatrixXdr::Zero(n_samples, 2 * n_randvecs);
  yXXy = MatrixXdr::Zero(2, 1);

  metaData metadata = set_metadata(n_samples, n_snps);
  ifstream bed_ifs(bed_file.c_str(), ios::in | ios::binary);
  global_snp_index = -1;

  MatrixXdr vec1;
  MatrixXdr w1;
  MatrixXdr w2;
  MatrixXdr w3;
  //////////// analytic se
  int total_bin_num = 2;
  wt = MatrixXdr::Zero(n_samples, total_bin_num + 1);
  vt = MatrixXdr::Zero(n_samples, (total_bin_num + 1) * (total_bin_num + 1));

  for (int block_index = 0; block_index < n_blocks; block_index++) {
    int read_Nsnp = (block_index < (n_blocks - 1))
                        ? (step_size)
                        : (step_size + step_size_rem);

    set_block_parameters(genotype_block, n_samples, block_sizes[block_index]);

    read_genotype_block(bed_ifs, read_Nsnp, genotype_block, n_samples, n_snps,
                        global_snp_index, metadata);

    int block_size = block_sizes[block_index];

    if (block_size != 0) {
      allelecount_stds.resize(block_size, 1);
      allelecount_means.resize(block_size, 1);
      for (int i = 0; i < block_size; i++) {
        allelecount_means(i, 0) =
            (double)genotype_block.columnsum[i] / n_samples;
        if (allelecount_means(i, 0) == 2 | allelecount_means(i, 0) == 0)
          cout << "mean :" << allelecount_means(i, 0) << endl;
      }

      for (int i = 0; i < block_size; i++) {
        allelecount_stds(i, 0) =
            1 / sqrt((allelecount_means(i, 0) *
                      (1 - (0.5 * allelecount_means(i, 0)))));
      }

      int p = genotype_block.n_snps;
      int n = genotype_block.n_samples;

      // Initial intermediate data structures
        allocate_memory(n_randvecs, genotype_block, hsegsize,
                        partialsums, sum_op, yint_e, yint_m, y_e, y_m);

        output = compute_XXz(block_size, random_vectors, allelecount_means,
                           allelecount_stds, pheno_mask, n_randvecs, n_samples,
                           focal_snp_local_index, sum_op, genotype_block,
                           yint_m, y_m, p, yint_e, y_e, partialsums, false);

      MatrixXdr scaled_pheno;

      bool in_gxg_block = (focal_snp_block == block_index);
      MatrixXdr env_all_zb = random_vectors.array().colwise() *
                             focal_snp_gtype.col(0).array();
      output_env =
          compute_XXz(block_size, env_all_zb, allelecount_means,
                      allelecount_stds, pheno_mask, n_randvecs, n_samples,
                      focal_snp_local_index, sum_op, genotype_block, yint_m,
                      y_m, p, yint_e, y_e, partialsums, in_gxg_block); // z

      output_env =
          output_env.array().colwise() * focal_snp_gtype.col(0).array();
      for (int z_index = 0; z_index < n_randvecs; z_index++) {
        XXz.col((n_randvecs) + z_index) += output_env.col(z_index);
      }

      scaled_pheno = pheno.array() * focal_snp_gtype.col(0).array();
      MatrixXdr temp = compute_XXy(
          block_size, scaled_pheno, allelecount_means, allelecount_stds,
          pheno_mask, focal_snp_local_index, n_samples, sum_op, genotype_block,
          yint_m, y_m, p, yint_e, y_e, partialsums,
          in_gxg_block);
      temp = temp.array() * focal_snp_gtype.col(0).array();
      wt.col(1) += temp;
      yXXy(1, 0) += compute_yXXy(block_size, scaled_pheno, allelecount_means,
                                 allelecount_stds, focal_snp_local_index,
                                 sum_op, genotype_block, yint_m, y_m, p,
                                 partialsums, in_gxg_block);

      wt.col(0) += compute_XXy(
          block_size, pheno, allelecount_means, allelecount_stds, pheno_mask,
          focal_snp_local_index, n_samples, sum_op, genotype_block, yint_m, y_m,
          p, yint_e, y_e, partialsums, false);

      for (int z_index = 0; z_index < n_randvecs; z_index++) {
        XXz.col(z_index) += output.col(z_index); /// save whole sample
      }

      yXXy(0, 0) +=
          compute_yXXy(block_size, pheno, allelecount_means, allelecount_stds,
                       focal_snp_local_index, sum_op, genotype_block, yint_m,
                       y_m, p, partialsums, false);

        MatrixXdr val_temp;
        for (int i = 0; i < (total_bin_num + 1); i++) {
            MatrixXdr val_temp = compute_XXy(
                    block_size, wt.col(i), allelecount_means, allelecount_stds, pheno_mask, focal_snp_local_index,
                    n_samples, sum_op, genotype_block, yint_m, y_m, p, yint_e, y_e,
                    partialsums, false);
            vt.col(i) += val_temp / n_snps_variance_component[0];
        }

        MatrixXdr scaled_vec;
        for (int i = 0; i < (total_bin_num + 1); i++) {
            scaled_vec = wt.col(i).array() * focal_snp_gtype.col(0).array();
            MatrixXdr temp = compute_XXy(block_size, scaled_vec, allelecount_means, allelecount_stds,
                                         pheno_mask, focal_snp_local_index,
                                         n_samples, sum_op, genotype_block, yint_m,
                                         y_m, p, yint_e, y_e, partialsums, false);
            temp = temp.array() * focal_snp_gtype.col(0).array();

            vt.col(((total_bin_num + 1)) + i) +=
                    temp / n_snps_variance_component[1];
        }

        deallocate_memory(hsegsize, partialsums, sum_op, yint_e, yint_m, y_e,
                          y_m, genotype_block);
    }
  }

  /// analytic se
  for (int i = 0; i < total_bin_num; i++) {
    wt.col(i) = wt.col(i) / n_snps_variance_component[i];
    double dx = (wt.col(i).array() * pheno.array()).sum();
  }
  wt.col(total_bin_num) = pheno;

  for (int i = 0; i < (total_bin_num + 1); i++) {
    vt.col((total_bin_num * (total_bin_num + 1)) + i) = wt.col(i);
  }

  int n_variance_components = 2;
  // compute the elements of normal equation:

  /// normal equations LHS
  MatrixXdr A_trs(n_variance_components, n_variance_components);
  MatrixXdr b_trk(n_variance_components, 1);
  MatrixXdr c_yky(n_variance_components, 1);

  MatrixXdr X_l(n_variance_components + 1, n_variance_components + 1);
  MatrixXdr Y_r(n_variance_components + 1, 1);
  int jack_index = n_blocks;
  MatrixXdr B1;
  MatrixXdr B2;
  MatrixXdr C1;
  MatrixXdr C2;
  double trkij;
  double yy = (pheno.array() * pheno.array()).sum();

  int Nindv_mask = pheno_mask.sum();
  int NC = Nindv_mask;

  MatrixXdr point_est;
  MatrixXdr herit_est;

  point_est.resize(n_variance_components + 1, 1);
  herit_est.resize(n_variance_components + 1, 1);

  double trkij_res1;
  double trkij_res2;
  double trkij_res3;
  double tk_res;

  for (int i = 0; i < n_variance_components; i++) {

    b_trk(i, 0) = Nindv_mask;

    if (i >= (n_variance_components - 1)) { // change nongen_Nbin to 1

      B1 = XXz.block(0, i * n_randvecs, n_samples, n_randvecs);
      B1 = random_vectors.array() * B1.array();
      b_trk(i, 0) = B1.sum() / n_snps_variance_component[i] / n_randvecs;
    }

    c_yky(i, 0) = yXXy(i, 0) / n_snps_variance_component[i];

    for (int j = i; j < n_variance_components; j++) {
      B1 = XXz.block(0, i * n_randvecs, n_samples, n_randvecs);
      B2 = XXz.block(0, j * n_randvecs, n_samples, n_randvecs);
      C1 = B1.array() * B2.array();
      C2 = C1.colwise().sum();
      trkij = C2.sum();

      trkij = trkij / n_snps_variance_component[i] /
              n_snps_variance_component[j] / n_randvecs;
      A_trs(i, j) = trkij;
      A_trs(j, i) = trkij;
    }
  }

  ////solve normal equation as :

  X_l << A_trs, b_trk, b_trk.transpose(), NC;
  Y_r << c_yky, yy;

  MatrixXdr herit = X_l.colPivHouseholderQr().solve(Y_r);

  for (int i = 0; i < (n_variance_components + 1); i++)
    point_est(i, 0) = herit(i, 0);

  double temp_sum = point_est.sum();
  double temp_sig = 0;
  for (int j = 0; j < n_variance_components; j++) {
    herit_est(j, 0) = point_est(j, 0) / temp_sum;
    temp_sig += herit_est(j, 0);
  }
  herit_est(n_variance_components, 0) = temp_sig;

  /////compute SE analytic version
  MatrixXdr e;
  e = MatrixXdr::Zero(n_samples, 1);
  for (int i = 0; i < (n_variance_components + 1); i++)
    e += herit(i, 0) * wt.col(i);

  double d = 0;
  for (int i = 0; i < n_variance_components; i++)
    d += herit(i, 0) * c_yky(i, 0);

  d += herit(n_variance_components, 0) * yy;

  MatrixXdr cov_q;
  cov_q.resize(n_variance_components + 1, n_variance_components + 1);
  for (int k = 0; k < (n_variance_components + 1); k++) {
    for (int l = 0; l < (n_variance_components + 1); l++) {
      double sum1 = 0;
      for (int t = 0; t < (n_variance_components + 1); t++) {

        MatrixXdr z1 = vt.col((t * (n_variance_components + 1)) + l);
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
  std::cout << "Execution time of main function: " << elapsed.count()
            << " seconds" << std::endl;
  return Rcpp::List::create(Rcpp::Named("Est") = point_est.col(0),
                            Rcpp::Named("SE") =
                                cov_sigma.diagonal().array().sqrt());
}

void
allocate_memory(int n_randvecs, const genotype &genotype_block, int &hsegsize,
                double *&partialsums, double *&sum_op, double *&yint_e,
                double *&yint_m, double **&y_e, double **&y_m) {
    hsegsize = genotype_block.segment_size_hori; // = log_3(n)
    int hsize = pow(3, hsegsize);
    int vsegsize = genotype_block.segment_size_ver; // = log_3(p)
    int vsize = pow(3, vsegsize);

    partialsums = new double[n_randvecs];
    sum_op = new double[n_randvecs];
    yint_e = new double[hsize * n_randvecs];
    yint_m = new double[hsize * n_randvecs];
    memset(yint_m, 0, hsize * n_randvecs * sizeof(double));
    memset(yint_e, 0, hsize * n_randvecs * sizeof(double));

    y_e = new double *[genotype_block.n_samples];
    for (int i = 0; i < genotype_block.n_samples; i++) {
      y_e[i] = new double[n_randvecs];
      memset(y_e[i], 0, n_randvecs * sizeof(double));
    }

    y_m = new double *[hsegsize];
    for (int i = 0; i < hsegsize; i++)
      y_m[i] = new double[n_randvecs];
}

void deallocate_memory(int hsegsize, double *partialsums, double *sum_op,
                       double *yint_e, double *yint_m, double **y_e,
                       double **y_m, genotype &genotype_block) {
    delete[] sum_op;
    delete[] partialsums;
    delete[] yint_e;
    delete[] yint_m;
    for (int i = 0; i < hsegsize; i++)
      delete[] y_m[i];
    delete[] y_m;

    for (int i = 0; i < genotype_block.n_samples; i++)
      delete[] y_e[i];
    delete[] y_e;

    std::vector<std::vector<int>>().swap(genotype_block.p);
    std::vector<std::vector<int>>().swap(genotype_block.not_O_j);
    std::vector<std::vector<int>>().swap(genotype_block.not_O_i);
    genotype_block.columnsum.clear();
    genotype_block.columnmeans.clear();
}
