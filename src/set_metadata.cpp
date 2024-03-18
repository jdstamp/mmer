#include "set_metadata.h"

metaData set_metadata(int n_samples, int n_snps) {
  int wordsize = sizeof(char) * 8;
  int unitsize = 2;
  unsigned int unitsperword = wordsize / unitsize;
  unsigned char mask = 0;
  for (int i = 0; i < unitsize; i++)
    mask = mask | (0x1 << i);
  int ncol = ceil(1.0 * n_samples / unitsperword);
  int nrow = n_snps;
  metaData data =
      metaData{wordsize, unitsize, unitsperword, mask, ncol, nrow, n_samples, n_snps};
  return (data);
}
