#include <RcppEigen.h>
#include <fstream>
#include <iostream>

#include "fame.h"
#include "set_metadata.h"

using namespace std;

void read_focal_snp(const string &filename, MatrixXdr &focal_genotype,
                    const int &focal_snp_index, const int &n_samples,
                    const int &n_snps, int &global_snp_index);

template <typename T>
static std::istream &binary_read(std::istream &stream, T &value);

float get_observed_allelefreq(const unsigned char *line,
                              const metaData &metadata);

int impute_genotype(const float &p_j);

void read_genotype_block(std::istream &ifs, const int &block_size,
                         genotype &genotype_block, const int &n_samples,
                         int &global_snp_index, const metaData &metadata);

void read_snp(std::istream &ifs, const int &n_samples, int &global_snp_index,
              const metaData &metadata, MatrixXdr &genotype_matrix);

void extract_plink_genotypes(int *y, const unsigned char &c,
                             const unsigned char &mask);

int encoding_to_allelecount(const int &value);

void normalize_genotype(MatrixXdr &focal_genotype, const int &n_samples);

int get_sample_block_size(const int &n_samples, const int &k, const int &ncol);

void encode_genotypes(genotype &genotype_block, int j, int val);

void read_masked_genotype_block(std::istream &ifs, const int &block_size,
                                genotype &genotype_block, const int &n_samples,
                                int &global_snp_index, const metaData &metadata,
                                const std::vector<int> &snp_indices);