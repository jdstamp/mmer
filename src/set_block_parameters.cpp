//
// Created by Julian Stamp on 3/22/24.
//

#include "read_annotation_file.h"
#include "genotype.h"
#include "fame.h"
#include "set_block_parameters.h"

void set_block_parameters(vector<genotype> &allgen_mail, int n_samples,
                          const annotationStruct &annotation, int block_index) {
  for (int i = 0; i < annotation.n_bin; i++) {
    allgen_mail[i].segment_size_hori =
        floor(log(n_samples) / log(3)) - 2; // object of the mailman
    allgen_mail[i].Nsegments_hori =
        ceil(annotation.jack_bin[block_index][i] * 1.0 /
             (allgen_mail[i].segment_size_hori * 1.0));
    allgen_mail[i].p.resize(allgen_mail[i].Nsegments_hori,
                            std::vector<int>(n_samples));
    allgen_mail[i].not_O_i.resize(annotation.jack_bin[block_index][i]);
    allgen_mail[i].not_O_j.resize(n_samples);
    allgen_mail[i].block_size = 0;
    allgen_mail[i].Nsnp = annotation.jack_bin[block_index][i];
    allgen_mail[i].Nindv = n_samples;

    allgen_mail[i].columnsum.resize(
        annotation.jack_bin[block_index][i],
        1);
    for (int index_temp = 0; index_temp < annotation.jack_bin[block_index][i];
         index_temp++) {
      allgen_mail[i].columnsum[index_temp] = 0;
    }
  }
}