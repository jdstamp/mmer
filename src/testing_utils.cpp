#include "testing_utils.h"

// this should come from one configuration file
std::string testdata_dir = "../../inst/testdata/";
std::string checkdata_dir = "../../famer/testdata/";

std::string test_bed = testdata_dir + "test.bed";
std::string test_csv = testdata_dir + "test.csv";
std::string test_pheno = testdata_dir + "test_h2_0.5.pheno";

bool is_test = true;

double tolerance = 1e-6;
int n_samples = 200;
int block_size = 10;
metaData metadata = set_metadata(n_samples, block_size);

MatrixXdr readCSVToMatrixXdr(const std::string &filename) {
  std::ifstream data(filename);
  std::string line;
  std::vector<double> values;
  int rows = 0;
  int cols = 0;

  // Skip the header line
  std::getline(data, line);

  // Read data, line by line
  while (std::getline(data, line)) {
    std::stringstream lineStream(line);
    std::string cell;
    int temp_cols = 0;

    // Read each cell
    while (std::getline(lineStream, cell, ',')) {
      values.push_back(std::stod(cell));
      temp_cols++;
    }

    // Update the number of columns
    if (cols == 0) {
      cols = temp_cols;
    }

    rows++;
  }

  // Convert the vector of values into a MatrixXdr
  MatrixXdr mat(rows, cols);
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      mat(i, j) = values[i * cols + j];
  return mat;
}

bool fileExists(const std::string &path) {
  std::ifstream file(path);
  return file.good();
}

void correctTestFiles(std::string &test_csv, std::string &test_bed,
                      std::string &test_pheno) {
  is_test = fileExists(test_csv);
  if (!is_test) {
    test_bed = checkdata_dir + "test.bed";
    test_csv = checkdata_dir + "test.csv";
    test_pheno = checkdata_dir + "test_h2_0.5.pheno";
  }
}
