#include <RcppEigen.h>
#include <fstream>
#include <iostream>

#include "fame.h"
#include "read_annotation_file.h"
#include "set_metadata.h"

using namespace std;

void read_focal_snp(const string &filename, MatrixXdr &focal_genotype, const int &focal_snp_index,
                    const int &n_samples, const int &n_snps, int &global_snp_index);

template <typename T>
static std::istream &binary_read(std::istream &stream, T &value);

float get_observed_allelefreq(const unsigned char *line, const metaData &metadata);

int impute_genotype(const float &p_j);

void read_genotype_block(std::istream &ifs, const int &num_snp,
                         vector <genotype> &allgen_mail, const int &n_samples,
                         const int &n_snps, int &global_snp_index,
                         const annotationStruct &annotation,
                         const metaData &metadata);

void extract_plink_genotypes(int *y, const unsigned char &c, const unsigned char &mask);

int encoding_to_allelecount(const int &value);

int get_sample_block_size(const int &n_samples, const int &k, const int &ncol);