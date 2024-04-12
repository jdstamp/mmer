#include "read_genotypes.h"

void read_focal_snp(const string &filename, MatrixXdr &focal_genotype,
                    const int &focal_snp_index, const int &n_samples,
                    const int &n_snps, int &global_snp_index) {
  ifstream ifs(filename.c_str(), ios::in | ios::binary);
  char magic[3];
  metaData metadata = set_metadata(n_samples, n_snps);
  unsigned char *gtype;
  gtype = new unsigned char[metadata.ncol];

  if (global_snp_index < 0) {
    binary_read(ifs, magic);
  }

  int y[4];

  for (int i = 0; i <= focal_snp_index; i++) {
    global_snp_index++;
    ifs.read(reinterpret_cast<char *>(gtype),
             metadata.ncol * sizeof(unsigned char));
    // skip to the block_size of interest
    if (i == (focal_snp_index)) {
      float p_j = get_observed_allelefreq(gtype, metadata);
      for (int k = 0; k < metadata.ncol; k++) {
        unsigned char c = gtype[k];
        unsigned char mask = metadata.mask;
        extract_plink_genotypes(y, c, mask);
        int j0 = k * metadata.unitsperword;
        int ncol = metadata.ncol;
        int lmax = get_sample_block_size(n_samples, k, ncol);
        for (int l = 0; l < lmax; l++) {
          int j = j0 + l;
          int val = encoding_to_allelecount(y[l]);
          // impute missing genotype
          val = (val == -1) ? impute_genotype(p_j) : val;
          focal_genotype(j, 0) = val;
        }
      }
    }
  }

    delete[] gtype;
}

void normalize_genotype(MatrixXdr &focal_genotype, const int &n_samples) {
    double mean_sel_snp = focal_genotype.array().sum() / n_samples;
    double sd_sel_snp = sqrt((mean_sel_snp * (1 - (0.5 * mean_sel_snp))));
    focal_genotype.array() = focal_genotype.array() - mean_sel_snp;
    focal_genotype.array() = focal_genotype.array() / sd_sel_snp;
}

int get_sample_block_size(const int &n_samples, const int &k, const int &ncol) {
  // Handle number of individuals not being a multiple of 4
  int lmax = 4;
  if (k == ncol - 1) {
    lmax = n_samples % 4;
    lmax = (lmax == 0) ? 4 : lmax;
  }
  return lmax;
}

template <typename T>
static std::istream &binary_read(std::istream &stream, T &value) {
  return stream.read(reinterpret_cast<char *>(&value), sizeof(T));
}

float get_observed_allelefreq(const unsigned char *line,
                              const metaData &metadata) {
  int y[4];
  int observed_sum = 0;
  int observed_ct = 0;
  for (int k = 0; k < metadata.ncol; k++) {
    unsigned char c = line[k];
    unsigned char mask = metadata.mask;
    extract_plink_genotypes(y, c, mask);
    int lmax = get_sample_block_size(metadata.n_samples, k, metadata.ncol);
    for (int l = 0; l < lmax; l++) {
      int val = encoding_to_allelecount(y[l]);
      if (val > -1) {
        observed_sum += val;
        observed_ct++;
      }
    }
  }
  return observed_sum * 0.5 / observed_ct;
}

void extract_plink_genotypes(int *y, const unsigned char &c,
                             const unsigned char &mask) {
  // Extract PLINK genotypes
  y[0] = (c)&mask;
  y[1] = (c >> 2) & mask;
  y[2] = (c >> 4) & mask;
  y[3] = (c >> 6) & mask;
}

int impute_genotype(const float &p_j) {
  float rval = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
  float dist_pj[3] = {(1 - p_j) * (1 - p_j), 2 * p_j * (1 - p_j), p_j * p_j};
  if (rval < dist_pj[0])
    return 0;
  else if (rval >= dist_pj[0] && rval < (dist_pj[0] + dist_pj[1]))
    return 1;
  else
    return 2;
}

void read_genotype_block(std::istream &ifs, const int &block_size,
                         genotype &genotype_block, const int &n_samples,
                         int &global_snp_index, const metaData &metadata) {
  char magic[3];

  unsigned char *gtype;
  gtype = new unsigned char[metadata.ncol];

  if (global_snp_index < 0) {
    binary_read(ifs, magic);
  }
  int y[4];

  for (int i = 0; i < block_size; i++) {
    global_snp_index++;
    ifs.read(reinterpret_cast<char *>(gtype),
             metadata.ncol * sizeof(unsigned char));
    // TODO: compute p_j in preprocessing for all SNPs and only look up when
    //  needed?
    //    float p_j = get_observed_allelefreq(gtype, metadata);

    for (int k = 0; k < metadata.ncol; k++) {
      unsigned char c = gtype[k];
      unsigned char mask = metadata.mask;
      extract_plink_genotypes(y, c, mask);
      int j0 = k * metadata.unitsperword;
      int ncol = metadata.ncol;
      int lmax = get_sample_block_size(n_samples, k, ncol);
      for (int l = 0; l < lmax; l++) {
        int j = j0 + l;
        int val = encoding_to_allelecount(y[l]);
        // impute missing genotype
        val = (val == -1)
                  ? impute_genotype(get_observed_allelefreq(gtype, metadata))
                  : val;
        int snp_index;
        snp_index = genotype_block.block_size;
        int horiz_seg_no = snp_index / genotype_block.segment_size_hori;
        genotype_block.p[horiz_seg_no][j] =
            3 * genotype_block.p[horiz_seg_no][j] + val;
        // computing sum for every snp to compute mean
        genotype_block.columnsum[snp_index] += val;
      }
    }

    genotype_block.block_size++;
  }
  delete[] gtype;
}

int encoding_to_allelecount(const int &value) {
  // Extract  PLINK coded genotype and convert into 0/1/2
  // PLINK coding:
  // 00->0
  // 01->missing
  // 10->1
  // 11->2
  switch (value) {
  case 0:
    return 0;
  case 1:
    return -1;
  case 2:
    return 1;
  case 3:
    return 2;
  default:
    // Handle invalid input
    return -1; // or any other default value you prefer
  }
}
