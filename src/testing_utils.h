#ifndef FAMER_TESTING_UTILS_H
#define FAMER_TESTING_UTILS_H

#endif // FAMER_TESTING_UTILS_H

#include "allocate_memory.h"
#include "computation.h"
#include "compute_block_stats.h"
#include "compute_covariance_q.h"
#include "compute_mom_components.h"
#include "fame.h"
#include "genotype.h"
#include "initialize_random_vectors.h"
#include "read_genotypes.h"
#include "read_phenotypes.h"
#include "set_block_parameters.h"
#include "set_metadata.h"

#include <fstream>
#include <iostream>
#include <unistd.h>

// this should come from one configuration file
extern std::string testdata_dir;
extern std::string test_bed;
extern std::string test_csv;
extern std::string test_pheno;

extern double tolerance;
extern int n_samples;
extern int block_size;
extern metaData metadata;

MatrixXdr readCSVToMatrixXdr(const std::string &filename);