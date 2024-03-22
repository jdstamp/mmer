#include <RcppEigen.h>
#include <fstream>
#include <iostream>

#include "fame.h"
#include "read_annotation_file.h"
#include "set_metadata.h"

using namespace std;

void read_focal_snp(string filename, MatrixXdr &focal_genotype, int focal_snp_index,
                    int n_samples, int n_snps, int &global_snp_index);

template <typename T>
static std::istream &binary_read(std::istream &stream, T &value);

float get_observed_allelefreq(const unsigned char *line, metaData metadata);

int impute_genotype(float p_j);

void read_genotype_block(std::istream &ifs, int num_snp, vector <genotype> &allgen_mail,
                         int n_samples, int n_snps, int &global_snp_index,
                         annotationStruct annotation);

void extract_plink_genotypes(int *y, unsigned char c, unsigned char mask);

int encoding_to_allelecount(int val);

int get_sample_block_size(int n_samples, int k, int ncol);