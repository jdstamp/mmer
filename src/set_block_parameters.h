//
// Created by Julian Stamp on 3/22/24.
//

#ifndef FAMER_SET_BLOCK_PARAMETERS_H
#define FAMER_SET_BLOCK_PARAMETERS_H

#include "read_annotation_file.h"

void set_block_parameters(vector<genotype> &allgen_mail, const int &n_samples,
                          const annotationStruct &annotation, const int &block_index);

#endif //FAMER_SET_BLOCK_PARAMETERS_H
