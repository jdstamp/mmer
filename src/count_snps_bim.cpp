#include "count_snps_bim.h"

int count_snps_bim(std::string filename) {
  std::ifstream bim_file(filename.c_str());
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
