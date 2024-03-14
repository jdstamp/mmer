#include "fame.h"
#include <RcppEigen.h>
#include <fstream>

using namespace std;
using namespace Eigen;

MatrixXdr read_phenotypes(int Nind, std::string filename, MatrixXdr &mask);
