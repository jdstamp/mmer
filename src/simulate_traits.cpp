#include "simulate_traits.h"
#include "mme.h"

std::vector<int> draw_random_ints(std::vector<int> numbers, int x) {

  // Initialize a random number generator
  std::random_device rd;
  std::mt19937 g(rd());

  // Shuffle the numbers
  std::shuffle(numbers.begin(), numbers.end(), g);

  // Resize the vector to keep only the first x elements
  numbers.resize(x);

  return numbers;
}

MatrixXdr draw_normal_effects(int n) {
  // Initialize a random number generator
  boost::mt19937 seedr;
  seedr.seed(std::time(0));
  boost::normal_distribution<> dist(0, 1);
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

MatrixXdr &scale_component(float target_variance, MatrixXdr &component) {
  float snp_variance =
      (component.array() - component.array().mean()).square().mean();
  component = component.array() * std::sqrt(target_variance / snp_variance);
  return component;
}

// [[Rcpp::export]]
Rcpp::List simulate_traits_cpp(std::string plink_file,
                               float additive_heritability,
                               float gxg_heritability,
                               std::vector<int> additive_snps,
                               std::vector<int> gxg_group_1,
                               std::vector<int> gxg_group_2) {

  std::string bim_file = plink_file + ".bim";
  std::string fam_file = plink_file + ".fam";
  std::string bed_file = plink_file + ".bed";
  int n_snps = count_snps_bim(bim_file);
  int n_samples = count_fam(fam_file);

  int n_additive_snps = additive_snps.size();
  MatrixXdr additive_effects = draw_normal_effects(n_additive_snps);

  int n_group_1 = gxg_group_1.size();
  int n_group_2 = gxg_group_2.size();
  int n_interactions = n_group_1 * n_group_2;
  MatrixXdr epistatic_effects = draw_normal_effects(n_interactions);

  MatrixXdr additive_component = MatrixXdr::Zero(n_samples, 1);
  MatrixXdr gxg_component = MatrixXdr::Zero(n_samples, 1);
  MatrixXdr error_component = draw_normal_effects(n_samples);

  for (int i = 0; i < n_additive_snps; ++i) {
    int snp_index = additive_snps[i];
    MatrixXdr snp_genotype;
    snp_genotype.resize(n_samples, 1);
    int global_snp_index = -1;
    read_focal_snp(bed_file, snp_genotype, snp_index, n_samples, n_snps,
                   global_snp_index);
    additive_component = additive_component.array() +
                         snp_genotype.array() * additive_effects(i, 0);
  }

  float g1_snp_gxg_heritability = gxg_heritability / n_group_1;
  for (int i = 0; i < n_group_1; ++i) {
    for (int j = 0; j < n_group_2; ++j) {
      int snp_index1 = gxg_group_1[i];
      int snp_index2 = gxg_group_2[j];
      MatrixXdr snp_gxg_component;
      MatrixXdr snp_genotype1, snp_genotype2;
      snp_genotype1.resize(n_samples, 1);
      snp_genotype2.resize(n_samples, 1);
      int global_snp_index = -1;
      read_focal_snp(bed_file, snp_genotype1, snp_index1, n_samples, n_snps,
                     global_snp_index);
      normalize_genotype(snp_genotype1, n_samples);
      global_snp_index = -1;
      read_focal_snp(bed_file, snp_genotype2, snp_index2, n_samples, n_snps,
                     global_snp_index);
      normalize_genotype(snp_genotype2, n_samples);
      snp_gxg_component = snp_genotype1.array() * snp_genotype2.array() *
                          epistatic_effects(i * n_group_2 + j, 0);
      snp_gxg_component =
          scale_component(g1_snp_gxg_heritability, snp_gxg_component);
      gxg_component = gxg_component.array() + snp_gxg_component.array();
    }
  }

  // scale the additive, epistatic, and error components to control variance of
  // the trait
  scale_component(1 - additive_heritability - gxg_heritability,
                  error_component);
  scale_component(additive_heritability, additive_component);
  scale_component(gxg_heritability, gxg_component);
  if (gxg_heritability <= 0) {
    gxg_component = MatrixXdr::Zero(n_samples, 1);
  }
  if (additive_heritability <= 0) {
    additive_component = MatrixXdr::Zero(n_samples, 1);
  }
  // sum the additive, epistatic, and error components to get the trait
  MatrixXdr trait = additive_component + gxg_component + error_component;

  // Compute the updated variances
  float additive_variance =
      (additive_component.array() - additive_component.array().mean())
          .square()
          .mean();
  float gxg_variance =
      (gxg_component.array() - gxg_component.array().mean()).square().mean();
  float error_variance =
      (error_component.array() - error_component.array().mean())
          .square()
          .mean();

  // return a Rcpp list with the trait and the epistatic snps
  return Rcpp::List::create(Rcpp::Named("trait") = trait,
                            Rcpp::Named("additive_variance") =
                                additive_variance,
                            Rcpp::Named("gxg_variance") = gxg_variance,
                            Rcpp::Named("error_variance") = error_variance);
}
