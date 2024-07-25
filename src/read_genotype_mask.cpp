#include "read_genotype_mask.h"

// MatrixXdr indices_to_binary(const std::vector<int> &indices, int length) {
//   MatrixXdr binaryVec = MatrixXdr::Zero(length, 1);
//   for (int index : indices) {
//     if (index >= 0 && index < length) {
//       binaryVec(index, 0) = 1;
//     }
//   }
//   return binaryVec;
// }

void read_genotype_mask(const std::string &genotype_mask_file, int n_snps,
                        int gxg_i, const std::string &gxg_h5_dataset,
                        MatrixXdr &genotype_mask, int &n_gxg_snps) {
  std::vector<int> genotype_mask_indices;
  if (genotype_mask_file != "") {
    H5Easy::File file(genotype_mask_file, H5Easy::File::ReadOnly);
    genotype_mask_indices = H5Easy::load<std::vector<int>>(
        file, gxg_h5_dataset + "/" + std::to_string(gxg_i));
    //    genotype_mask = indices_to_binary(genotype_mask_indices, n_snps);
    for (int index : genotype_mask_indices) {
      if (index >= 0 && index < n_snps) {
        genotype_mask(index, 0) = 1;
      }
    }
    n_gxg_snps = genotype_mask.sum();
  } else {
    genotype_mask = MatrixXdr::Ones(n_snps, 1);
    n_gxg_snps = n_snps - 1;
  }
}