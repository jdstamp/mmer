#include "read_genotypes.h"

void read_genotypes(string filename, MatrixXdr &genotypes, bool allow_missing,
                    int num_snp, int n_samples, int &global_snp_index) {
  ifstream ifs(filename.c_str(), ios::in | ios::binary);
  char magic[3];
  metaData metadata = set_metadata(n_samples);
  unsigned char *gtype;
  gtype = new unsigned char[metadata.ncol];

  // if(read_header)
  binary_read(ifs, magic);

  int sum = 0;

  // Note that the coding of 0 and 2 can get flipped relative to plink because
  // plink uses allele frequency (minor) allele to code a SNP as 0 or 1. This
  // flipping does not matter for results.
  int y[4];

  for (int i = 0; i < num_snp; i++) {
    global_snp_index++;
    ifs.read(reinterpret_cast<char *>(gtype),
             metadata.ncol * sizeof(unsigned char));
    float p_j = get_observed_pj(gtype, metadata);
    for (int k = 0; k < metadata.ncol; k++) {
      unsigned char c = gtype[k];
      // Extract PLINK genotypes
      y[0] = (c)&metadata.mask;
      y[1] = (c >> 2) & metadata.mask;
      y[2] = (c >> 4) & metadata.mask;
      y[3] = (c >> 6) & metadata.mask;
      int j0 = k * metadata.unitsperword;
      // Handle number of individuals not being a multiple of 4
      int lmax = 4;
      if (k == metadata.ncol - 1) {
        lmax = n_samples % 4;
        lmax = (lmax == 0) ? 4 : lmax;
      }
      for (int l = 0; l < lmax; l++) {
        int j = j0 + l;
        // Extract  PLINK coded genotype and convert into 0/1/2
        // PLINK coding:
        // 00->0
        // 01->missing
        // 10->1
        // 11->2
        int val = y[l];
        if (val == 1 && !allow_missing) {
          val = simulate2_geno_from_random(p_j);
          val++;
          val = (val == 1) ? 0 : val;
          //                                 val=0;
        }
        val--;
        val = (val < 0) ? 0 : val;
        sum += val;

        if (i == (num_snp - 1))
          genotypes(j, 0) = val;
      }
    }
  }

  sum = 0;
  delete[] gtype;
}

template <typename T>
static std::istream &binary_read(std::istream &stream, T &value) {
  return stream.read(reinterpret_cast<char *>(&value), sizeof(T));
}

float get_observed_pj(const unsigned char *line, metaData metadata) {
  int y[4];
  int observed_sum = 0;
  int observed_ct = 0;

  for (int k = 0; k < metadata.ncol; k++) {
    unsigned char c = line[k];
    y[0] = (c)&metadata.mask;
    y[1] = (c >> 2) & metadata.mask;
    y[2] = (c >> 4) & metadata.mask;
    y[3] = (c >> 6) & metadata.mask;
    int j0 = k * metadata.unitsperword;
    int lmax = 4;
    if (k == metadata.ncol - 1) {
      lmax = metadata.nrow % 4;
      lmax = (lmax == 0) ? 4 : lmax;
    }
    for (int l = 0; l < lmax; l++) {
      int j = j0 + l;
      // Extract  PLINK coded genotype and convert into 0/1/2
      // // PLINK coding:
      // // 00->0
      // // 01->missing
      // // 10->1
      // // 11->2
      int val = y[l];
      val--;
      if (val != 0) {
        val = (val < 0) ? 0 : val;
        observed_sum += val;
        observed_ct++;
      }
    }
  }
  return observed_sum * 0.5 / observed_ct;
}

int simulate2_geno_from_random(float p_j) {
  float rval = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
  float dist_pj[3] = {(1 - p_j) * (1 - p_j), 2 * p_j * (1 - p_j), p_j * p_j};
  if (rval < dist_pj[0])
    return 0;
  else if (rval >= dist_pj[0] && rval < (dist_pj[0] + dist_pj[1]))
    return 1;
  else
    return 2;
}

void read_bed2(std::istream &ifs, bool allow_missing, int num_snp,
               vector<genotype> &allgen_mail, int n_samples,
               int &global_snp_index, annotationStruct annotation,
               int selected_snp_index) {
  // ifstream ifs (filename.c_str(), ios::in|ios::binary);

  bool read_sel_snp = false;
  char magic[3];
  metaData metadata = set_metadata(n_samples);

  unsigned char *gtype;
  gtype = new unsigned char[metadata.ncol];

  // if (read_header)
  binary_read(ifs, magic);

  int sum = 0;

  // Note that the coding of 0 and 2 can get flipped relative to plink because
  // plink uses allele frequency (minor) allele to code a SNP as 0 or 1. This
  // flipping does not matter for results.
  int y[4];

  int bin_pointer;

  vector<int> pointer_bins;
  cout << "In read_bed2: about to read: " << num_snp << "number of snps"
       << endl;
  for (int i = 0; i < num_snp;
       i++) {           // Boyang: num_snp: snps in the current bin
    global_snp_index++; // Boyang: global_snp_index starts from 0

    ifs.read(reinterpret_cast<char *>(gtype),
             metadata.ncol * sizeof(unsigned char));
    float p_j = get_observed_pj(gtype, metadata);

    pointer_bins.clear();
    for (int bin_index = 0; bin_index < annotation.n_bin; bin_index++)
      if (annotation.annot_bool[global_snp_index][bin_index] ==
          1) // Boyang: checking annotation !!!
        pointer_bins.push_back(bin_index);

    for (int k = 0; k < metadata.ncol; k++) {
      unsigned char c = gtype[k];
      // Extract PLINK genotypes
      y[0] = (c)&metadata.mask;
      y[1] = (c >> 2) & metadata.mask;
      y[2] = (c >> 4) & metadata.mask;
      y[3] = (c >> 6) & metadata.mask;
      int j0 = k * metadata.unitsperword;
      // Handle number of individuals not being a multiple of 4
      int lmax = 4;
      if (k == metadata.ncol - 1) {
        lmax = n_samples % 4;
        lmax = (lmax == 0) ? 4 : lmax;
      }
      for (int l = 0; l < lmax; l++) {
        int j = j0 + l;
        // Extract  PLINK coded genotype and convert into 0/1/2
        // PLINK coding:
        // 00->0
        // 01->missing
        // 10->1
        // 11->2
        int val = y[l];
        if (val == 1 && !allow_missing) {
          val = simulate2_geno_from_random(p_j);
          val++;
          val = (val == 1) ? 0 : val;
          // val=0;
        }
        val--;
        val = (val < 0) ? 0 : val;
        sum += val;
        for (int bin_index = 0; bin_index < pointer_bins.size();
             bin_index++) { // !!!
          bin_pointer = pointer_bins[bin_index];
          int snp_index;
          //       cout << "in read_bed2: bin_index: " << bin_index <<
          //       "snp_index: "<< allgen_mail[bin_pointer].index << endl; //
          //       step by step update

          snp_index = allgen_mail[bin_pointer].index;
          int horiz_seg_no =
              snp_index / allgen_mail[bin_pointer].segment_size_hori;
          allgen_mail[bin_pointer].p[horiz_seg_no][j] =
              3 * allgen_mail[bin_pointer].p[horiz_seg_no][j] + val;
          // computing sum for every snp to compute mean
          allgen_mail[bin_pointer].columnsum[snp_index] += val;

          // Ali's index fix
          int sel_snp_local_index;
          if (global_snp_index == (selected_snp_index - 1) &
              read_sel_snp == false) {
            sel_snp_local_index = snp_index;
            cout << "################" << endl;
            cout << "current global_snp_index is: " << global_snp_index << endl;
            for (int itmp = 0; itmp < pointer_bins.size(); itmp++)
              cout << "pointer_bins  is " << pointer_bins[itmp]
                   << "size is: " << allgen_mail[pointer_bins[itmp]].index
                   << endl;
            cout << "Adjust bin_index is: " << bin_index << endl;
            cout << "Adjusted local index of target SNP: "
                 << sel_snp_local_index << endl;
            // sel_snp_bin=bin_index;
            read_sel_snp = true;
          }
        }
      }
    }

    for (int bin_index = 0; bin_index < pointer_bins.size(); bin_index++) {
      bin_pointer = pointer_bins[bin_index];
      allgen_mail[bin_pointer].index++;
    }
  }

  sum = 0;
  delete[] gtype;
}
