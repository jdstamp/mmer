#include "fame.h"
#include <string>

int main() {
  std::string plink_file = "/Users/jds/Downloads/test100k/test100k";
  std::string pheno_file = "/Users/jds/Downloads/test100k/h2_0.5.large.pheno";
  std::string annotation_file = "/Users/jds/Downloads/test100k/happnest.anno.txt";
  std::string covariate_file = "";
  int n_randvecs =100;
  int focal_snp_index = 1;
  int n_blocks = 100;
  Rcpp::List results = fame_cpp(plink_file,
                                pheno_file,
                                annotation_file,
                                covariate_file,
                                n_randvecs,
                                focal_snp_index,
                                n_blocks);
  return 0;
}
