#include "genotype.h"

using namespace std;

template <typename T>
static std::istream &binary_read(std::istream &stream, T &value) {
  return stream.read(reinterpret_cast<char *>(&value), sizeof(T));
}

double genotype::get_col_mean(int snpindex) const {
  double temp = columnmeans[snpindex];
  return temp;
}

double genotype::get_col_std(int snpindex) const {
  double p_i = get_col_mean(snpindex);
  if (p_i == 0 || p_i == 2)
    return 1.0;
  double temp = sqrt(p_i * (1 - (0.5 * p_i)));
  return temp;
}
