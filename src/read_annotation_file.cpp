#include "read_annotation_file.h"

annotationStruct read_annotation_file(string filename, int Nsnp, int Njack) {
  vector<bool> snp_annot;
  std::vector<std::vector<int>> jack_bin;
  vector<int> len;
  std::vector<std::vector<bool>> annot_bool;
  vector<int> selected_snps_vec;
  int step_size;
  int Nbin = 160;
  int step_size_rem;
  ifstream inp(filename.c_str());
  if (!inp.is_open()) {
    cerr << "Error reading file " << filename << endl;
    exit(1);
  }
  string line;
  int j = 0;
  int linenum = 0;
  int num_parti;
  stringstream check1(line);
  string intermediate;
  vector<string> tokens;
  while (std::getline(inp, line)) {
    char c = line[0];
    if (c == '#')
      continue;
    istringstream ss(line);
    if (line.empty())
      continue;
    j++;
    // cout<<line<<endl;

    stringstream check1(line);
    string intermediate;
    vector<string> tokens;
    // Tokenizing w.r.t. space ' '
    while (getline(check1, intermediate, ' ')) {
      tokens.push_back(intermediate);
    }
    if (linenum == 0) {
      num_parti = tokens.size();
      Nbin = num_parti;
      snp_annot.resize(Nbin, 0);
      jack_bin.resize(Njack, vector<int>(Nbin, 0));
      len.resize(num_parti, 0);
    }
    int index_annot = 0;
    for (int i = 0; i < tokens.size(); i++) {
      snp_annot[i] = 0;
      if (tokens[i] == "1") {
        len[i]++;
        snp_annot[i] = 1;
      }
    }
    annot_bool.push_back(snp_annot);
    linenum++;
  }

  if (Nsnp != linenum) {
    cout << "Number of the rows in bim file and annotation file does not match"
         << endl;
  }

  Nsnp = linenum;
  selected_snps_vec.resize(num_parti, 0);
  for (int i = 0; i < num_parti; i++) {
    cout << len[i] << " SNPs in " << i << "-th bin" << endl;
    selected_snps_vec[i] =
        len[i]; // Boyang: change selected_snps to selected_snps_vec
    cout << " Number of selected SNPs w.r.t  annot file : "
         << selected_snps_vec[i] << endl;
  }

  step_size = Nsnp / Njack;
  step_size_rem = Nsnp % Njack;
  cout << "Number of SNPs per block : " << step_size << endl;
  jack_bin.resize(Njack, vector<int>(Nbin, 0));
  int temp;
  for (int i = 0; i < Nsnp; i++)
    for (int j = 0; j < Nbin; j++)
      if (annot_bool[i][j] == 1) {
        temp = i / step_size;
        if (temp >= Njack)
          temp = Njack - 1;
        jack_bin[temp][j]++;
      }
  annotationStruct annotation =
      annotationStruct{step_size,         Nbin,          annot_bool, len,
                       selected_snps_vec, step_size_rem, jack_bin};
  return (annotation);
}
