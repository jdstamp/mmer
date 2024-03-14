#include <cmath>

struct metaData {
  int wordsize;
  int unitsize;
  unsigned int unitsperword;
  unsigned char mask;
  int ncol;
  int nrow;
};

metaData set_metadata(int n_samples);
