#include "fame.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

int read_covariates(bool std, int Nind, std::string filename,
                    std::string covname, MatrixXdr &covariate, bool snp_fix_ef);
