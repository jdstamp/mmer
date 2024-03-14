#include "read_phenotypes.h"

MatrixXdr read_phenotypes(int Nind, std::string filename, MatrixXdr &mask) {
  MatrixXdr pheno;

  ifstream ifs(filename.c_str(), ios::in);

  std::string line;
  std::istringstream in;
  int phenocount = 0;
  // read header
  std::getline(ifs, line);
  in.str(line);
  string b;
  while (in >> b) {
    if (b != "FID" && b != "IID")
      phenocount++;
  }
  pheno.resize(Nind, phenocount);
  mask.resize(Nind, phenocount);
  int i = 0;
  while (std::getline(ifs, line)) {
    in.clear();
    in.str(line);
    string temp;
    in >> temp;
    in >> temp;
    for (int j = 0; j < phenocount; j++) {
      in >> temp;
      double cur = atof(temp.c_str());
      if (temp == "NA" || cur == -9) {
        pheno(i, j) = 0;
        mask(i, j) = 0;
      } else {
        pheno(i, j) = atof(temp.c_str());
        mask(i, j) = 1;
      }
    }
    i++;
  }
  return (pheno);
}
