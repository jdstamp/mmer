#include "count_samples.h"

int count_samples(std::string filename) {
  ifstream ifs(filename.c_str(), ios::in);

  std::string line;
  int i = 0;
  while (std::getline(ifs, line)) {
    i++;
  }
  int n_samples = i - 1;
  return (n_samples);
}

int count_fam(std::string filename) {
  ifstream ifs(filename.c_str(), ios::in);

  std::string line;
  int i = 0;
  while (std::getline(ifs, line)) {
    i++;
  }
  return i;
}
