#pragma once
#include <highfive/H5Easy.hpp>
#include <string>

#include "fame.h"

void read_genotype_mask(const std::string &genotype_mask_file, int n_snps,
                        int gxg_i, const std::string &gxg_h5_dataset,
                        MatrixXdr &genotype_mask, int &n_gxg_snps);