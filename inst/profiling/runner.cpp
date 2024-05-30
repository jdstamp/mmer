#include "fame.h"
#include <string>

int main() {
  std::string plink_file = "/Users/jds/Downloads/test100k/subset_20000";
  std::string pheno_file = "/Users/jds/Downloads/test100k/subset_20000.phen";
  std::string covariate_file = "";
  std::string mask_file = ""; // "/Users/jds/Downloads/test100k/hdf5_mask.h5"
  int n_randvecs = 100;
  int focal_snp_index = 1;
  int n_blocks = 100;
  int rand_seed = 123;
  std::vector<int> index_vector = {focal_snp_index};
  Rcpp::List results = fame_cpp(plink_file, pheno_file, covariate_file,
                                n_randvecs, n_blocks, rand_seed, index_vector,
                                mask_file);
  return 0;
}
