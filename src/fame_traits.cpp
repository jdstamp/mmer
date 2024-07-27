#include "fame_traits.h"
#include "mme.h"
#include <highfive/H5Easy.hpp>

MatrixXdr draw_normal_effects(int n, float variance_component) {
  // Initialize a random number generator
  boost::mt19937 seedr;
  seedr.seed(std::time(0));
  boost::normal_distribution<> dist(0, std::sqrt(variance_component));
  boost::variate_generator<boost::mt19937 &, boost::normal_distribution<>>
      effect_size(seedr, dist);

  // Create a matrix to store the samples
  MatrixXdr samples(n, 1);

  // Draw the samples
  for (int i = 0; i < n; ++i) {
    samples(i, 0) = effect_size();
  }

  return samples;
}

// [[Rcpp::export]]
Rcpp::List fame_traits_cpp(std::string plink_file,
                           float additive_variance_component,
                           float gxg_variance_component, int focal_snp,
                           std::vector<int> gxg_indices,
                           std::vector<int> additive_indices) {

  std::string bim_file = plink_file + ".bim";
  std::string fam_file = plink_file + ".fam";
  std::string bed_file = plink_file + ".bed";
  int n_snps = count_snps_bim(bim_file);
  int n_samples = count_fam(fam_file);
  int n_additive = additive_indices.size();
  int n_gxg = gxg_indices.size();

  MatrixXdr additive_effects =
      draw_normal_effects(n_additive, additive_variance_component / n_additive);
  MatrixXdr epistatic_effects =
      draw_normal_effects(n_gxg, gxg_variance_component / n_gxg);

  // initialize the trait components for additive, epistatic and error
  // components
  MatrixXdr additive_component = MatrixXdr::Zero(n_samples, 1);
  MatrixXdr epistatic_component = MatrixXdr::Zero(n_samples, 1);
  MatrixXdr error_component = draw_normal_effects(
      n_samples, 1 - additive_variance_component - gxg_variance_component);

  // loop over the additive SNPs vector. For each SNP read the genotype
  // from the bed file and compute the dot product with the additive effect
  // corresponding to that SNP
  for (int i = 0; i < n_additive; ++i) {
    int snp_index = additive_indices[i];
    MatrixXdr snp_genotype;
    snp_genotype.resize(n_samples, 1);
    int global_snp_index = -1;
    read_focal_snp(bed_file, snp_genotype, snp_index, n_samples, n_snps,
                   global_snp_index);
    normalize_genotype(snp_genotype, n_samples);
    additive_component = additive_component.array() +
                         snp_genotype.array() * additive_effects(i, 0);
  }

  MatrixXdr focal_snp_genotype, epistatic_snp_genotype;
  focal_snp_genotype.resize(n_samples, 1);
  epistatic_snp_genotype.resize(n_samples, 1);
  int global_snp_index = -1;
  read_focal_snp(bed_file, focal_snp_genotype, focal_snp, n_samples, n_snps,
                 global_snp_index);
  normalize_genotype(focal_snp_genotype, n_samples);

  for (int i = 0; i < n_gxg; ++i) {
    int gxg_snp_index = gxg_indices[i];
    global_snp_index = -1;
    read_focal_snp(bed_file, epistatic_snp_genotype, gxg_snp_index, n_samples,
                   n_snps, global_snp_index);
    normalize_genotype(epistatic_snp_genotype, n_samples);
    epistatic_component =
        epistatic_component.array() + focal_snp_genotype.array() *
                                          epistatic_snp_genotype.array() *
                                          epistatic_effects(i, 0);
  }

  // compute the variances of the additive, the epistatic, and the error
  // components
  float additive_variance =
      (additive_component.array() - additive_component.array().mean())
          .square()
          .mean();
  float epistatic_variance =
      (epistatic_component.array() - epistatic_component.array().mean())
          .square()
          .mean();
  float error_variance =
      (error_component.array() - error_component.array().mean())
          .square()
          .mean();

  // sum the additive, epistatic, and error components to get the trait
  MatrixXdr trait = additive_component + epistatic_component + error_component;

  // Compute the updated variances
  additive_variance =
      (additive_component.array() - additive_component.array().mean())
          .square()
          .mean();
  epistatic_variance =
      (epistatic_component.array() - epistatic_component.array().mean())
          .square()
          .mean();
  error_variance = (error_component.array() - error_component.array().mean())
                       .square()
                       .mean();

  // print to Rccp console
  Rcpp::Rcout << "Additive variance: " << additive_variance << std::endl;
  Rcpp::Rcout << "Epistatic variance: " << epistatic_variance << std::endl;
  Rcpp::Rcout << "Error variance: " << error_variance << std::endl;

  // return a Rcpp list with the trait and the epistatic snps
  return Rcpp::List::create(
      Rcpp::Named("trait") = trait,
      Rcpp::Named("additive_variance") = additive_variance,
      Rcpp::Named("epistatic_variance") = epistatic_variance,
      Rcpp::Named("error_variance") = error_variance);
}
