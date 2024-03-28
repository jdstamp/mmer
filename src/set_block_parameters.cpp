//
// Created by Julian Stamp on 3/22/24.
//
#include "set_block_parameters.h"

void set_block_parameters(std::vector <genotype> &allgen_mail, const int
&n_samples,
                          const int &block_size) {
  int i = 0;
    allgen_mail[i].segment_size_hori =
        floor(log(n_samples) / log(3)) - 2; // object of the mailman
    allgen_mail[i].Nsegments_hori =
        ceil(block_size * 1.0 /
             (allgen_mail[i].segment_size_hori * 1.0));
    allgen_mail[i].p.resize(allgen_mail[i].Nsegments_hori,
                            std::vector<int>(n_samples));
    allgen_mail[i].not_O_i.resize(block_size);
    allgen_mail[i].not_O_j.resize(n_samples);
    allgen_mail[i].block_size = 0;
    allgen_mail[i].Nsnp = block_size;
    allgen_mail[i].Nindv = n_samples;

    allgen_mail[i].columnsum.resize(
            block_size,
        1);
    for (int index_temp = 0; index_temp < block_size;
         index_temp++) {
      allgen_mail[i].columnsum[index_temp] = 0;
    }
}