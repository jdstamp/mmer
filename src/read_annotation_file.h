#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

struct annotationStruct {
  int step_size;
  int n_bin;
  std::vector<std::vector<bool>> annot_bool;
  vector<int> len;
  vector<int> selected_snps_vec;
  int step_size_rem;
  std::vector<std::vector<int>> jack_bin;
};

annotationStruct read_annotation_file(string filename, int n_snps, int n_blocks);
