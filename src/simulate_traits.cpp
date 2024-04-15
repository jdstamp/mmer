#include "simulate_traits.h"
#include "fame.h"
#include <highfive/H5Easy.hpp>

std::vector<int> draw_random_ints(std::vector<int> numbers, int x) {

    // Initialize a random number generator
    std::random_device rd;
    std::mt19937 g(rd());

    // Shuffle the numbers
    std::shuffle(numbers.begin(), numbers.end(), g);

    // Resize the vector to keep only the first x elements
    numbers.resize(x);

    return numbers;
}

MatrixXdr draw_normal_effects(int n) {
    // Initialize a random number generator
    boost::mt19937 seedr;
    seedr.seed(std::time(0));
    boost::normal_distribution<> dist(0,1);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > effect_size(seedr, dist);

    // Create a matrix to store the samples
    MatrixXdr samples(n, 1);

    // Draw the samples
    for (int i = 0; i < n; ++i) {
        samples(i, 0) = effect_size();
    }

    return samples;
}

// [[Rcpp::export]]
std::vector<int> readH5File(const std::string& filename, const std::string& datasetName) {
    // Open the HDF5 file in read-only mode
    H5Easy::File file(filename, H5Easy::File::ReadOnly);
    // Load the dataset
    std::vector<int> dataset = H5Easy::load<std::vector<int>>(file, datasetName);
    return dataset;
}

// [[Rcpp::export]]
void replaceH5Dataset(const std::string& filename, const std::string& datasetName, const std::vector<int>& newData) {
    // Open the HDF5 file in read/write mode
    H5Easy::File file(filename, H5Easy::File::ReadWrite);
    // Overwrite the dataset with the new data
    H5Easy::dump(file, datasetName, newData, H5Easy::DumpMode::Overwrite);
}

// [[Rcpp::export]]
Rcpp::List
simulate_traits_cpp(std::string plink_file, float heritability, float rho,
                    int n_additive_snps, std::vector<int> gxg_group_1,
                    std::vector<int> gxg_group_2) {

    std::string bim_file = plink_file + ".bim";
    std::string fam_file = plink_file + ".fam";
    std::string bed_file = plink_file + ".bed";
    int n_snps = count_snps_bim(bim_file);
    int n_samples = count_fam(fam_file);
    int n_group_1 = gxg_group_1.size();
    int n_group_2 = gxg_group_2.size();
    int n_remaining_causal = n_additive_snps - n_group_1 - n_group_2;

    std::vector<int> epistatic_snps = gxg_group_1;
    epistatic_snps.insert(epistatic_snps.end(), gxg_group_2.begin(), gxg_group_2.end());
    std::sort(epistatic_snps.begin(), epistatic_snps.end());
    std::vector<int> difference;
    std::vector<int> numbers(n_snps);
    std::iota(numbers.begin(), numbers.end(), 0); // Fill with 0, 1, ..., n-1

    std::set_difference(numbers.begin(), numbers.end(), epistatic_snps.begin(), epistatic_snps.end(), std::back_inserter(difference));

    std::vector<int> additive_snps = draw_random_ints(difference, n_remaining_causal);
    additive_snps.insert(additive_snps.end(), epistatic_snps.begin(), epistatic_snps.end());
    std::sort(additive_snps.begin(), additive_snps.end());

    MatrixXdr additive_effects = draw_normal_effects(n_additive_snps);
    int n_interactions = n_group_1 * n_group_2;
    MatrixXdr epistatic_effects = draw_normal_effects(n_interactions);

    // initialize the trait components for additive, epistatic and error
    // components
    MatrixXdr additive_component = MatrixXdr::Zero(n_samples, 1);
    MatrixXdr epistatic_component = MatrixXdr::Zero(n_samples, 1);
    MatrixXdr error_component = draw_normal_effects(n_samples);

    // loop over the additive SNPs vector. For each SNP read the genotype
    // from the bed file and compute the dot product with the additive effect
    // corresponding to that SNP
    for (int i = 0; i < n_additive_snps; ++i) {
        int snp_index = additive_snps[i];
        MatrixXdr snp_genotype;
        snp_genotype.resize(n_samples, 1);
        int global_snp_index = -1;
        read_focal_snp(bed_file, snp_genotype, snp_index, n_samples, n_snps,
                       global_snp_index);
        additive_component = additive_component.array() + snp_genotype.array() *
                additive_effects(i, 0);
    }

    // loop over the first half of the epistatic SNPs vector and in an inner
    // loop over the second half of the epistatic SNP vector in an upper
    // triangular fashion
    for (int i = 0; i < n_group_1; ++i) {
        for (int j = 0; j < n_group_2; ++j) {
            int snp_index1 = gxg_group_1[i];
            int snp_index2 = gxg_group_2[j];
            MatrixXdr snp_genotype1, snp_genotype2;
            snp_genotype1.resize(n_samples, 1);
            snp_genotype2.resize(n_samples, 1);
            int global_snp_index = -1;
            read_focal_snp(bed_file, snp_genotype1, snp_index1, n_samples, n_snps,
                           global_snp_index);
            normalize_genotype(snp_genotype1, n_samples);
            global_snp_index = -1;
            read_focal_snp(bed_file, snp_genotype2, snp_index2, n_samples, n_snps,
                           global_snp_index);
            normalize_genotype(snp_genotype2, n_samples);
            epistatic_component = epistatic_component.array()
                    + snp_genotype1.array() * snp_genotype2.array() * epistatic_effects(i*n_group_2 + j, 0);
        }
    }

    // compute the variances of the additive, the epistatic, and the error
    // components
    float additive_variance = (additive_component.array() -
            additive_component.array().mean()).square().mean();
    float epistatic_variance = (epistatic_component.array() -
            epistatic_component.array().mean()).square().mean();
    float error_variance = (error_component.array() -
            error_component.array().mean()).square().mean();

    // scale the additive, epistatic, and error components to control variance of the trait
    error_component = error_component.array() * std::sqrt((1 - heritability) / error_variance);
    additive_component = additive_component.array() * std::sqrt(heritability * rho / additive_variance);
    epistatic_component = epistatic_component.array() * std::sqrt(heritability * (1 - rho) / epistatic_variance);

    // sum the additive, epistatic, and error components to get the trait
    MatrixXdr trait = additive_component + epistatic_component + error_component;

    // Compute the updated variances
    additive_variance = (additive_component.array() -
                       additive_component.array().mean()).square().mean();
    epistatic_variance = (epistatic_component.array() -
                       epistatic_component.array().mean()).square().mean();
    error_variance = (error_component.array() -
                       error_component.array().mean()).square().mean();

    // return a Rcpp list with the trait and the epistatic snps
    return Rcpp::List::create(Rcpp::Named("trait") = trait,
                              Rcpp::Named("epistatic_snps") = epistatic_snps,
                              Rcpp::Named("additive_variance") = additive_variance,
                              Rcpp::Named("epistatic_variance") = epistatic_variance,
                              Rcpp::Named("error_variance") = error_variance);
}
