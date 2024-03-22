#include "count_data.h"

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
  std::ifstream ifs(filename.c_str(), ios::in);

  std::string line;
  int i = 0;
  while (std::getline(ifs, line)) {
    i++;
  }
  return i;
}

int count_snps_bim(std::string filename) {
  std::ifstream bim_file(filename.c_str(), ios::in);
  if (!bim_file.is_open()) {
    std::cerr << "Error reading file " << filename << std::endl;
    exit(1);
  }
  std::string line;
  int n_snps = 0;
  int linenum = 0;
  while (std::getline(bim_file, line)) {
    linenum++;
    char c = line[0];
    if (c == '#')
      continue;
    if (line.empty())
      continue;
    n_snps++;
  }
  bim_file.close();
  return (n_snps);
}
