#include <RcppEigen.h>
#include <fstream>
#include <iostream>

#include "fame.h"
#include "read_annotation_file.h"
#include "set_metadata.h"

using namespace std;

void read_genotypes(string filename, MatrixXdr &genotypes, bool allow_missing,
                    int num_snp, int n_samples, int n_snps, int &global_snp_index);

template <typename T>
static std::istream &binary_read(std::istream &stream, T &value);

float get_observed_pj(const unsigned char *line, metaData metadata);

int simulate2_geno_from_random(float p_j);

void read_bed2(std::istream &ifs, bool allow_missing, int num_snp,
               vector<genotype> &allgen_mail, int n_samples, int n_snps,
               int &global_snp_index, annotationStruct annotation,
               int selected_snp_index);
