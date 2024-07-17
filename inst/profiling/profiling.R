library(famer)

plink_file <- "/Users/jds/data/ukbb/c12_100k-samples_010k-snps_imputed"
pheno_file <- "/Users/jds/data/ukbb/c12_100k-samples_010k-snps_imputed.phen"
# plink_file <- "/users/jstamp1/data/jstamp1/ukbb/c12_100k-samples_010k-snps_imputed"
# pheno_file <- "/users/jstamp1/data/jstamp1/ukbb/c12_100k-samples_010k-snps_imputed.phen"
covariate_file <- ""
mask_file <- "" # "/Users/jds/Downloads/test100k/hdf5_mask.h5"
n_randvecs <- 10
n_blocks <- 100
rand_seed <- 123

chunksize <- 30
n_threads <- 10
log_level <- "DEBUG"

snp_indices <- 1:1000

Rprof(memory.profiling = TRUE)

result <- fame(plink_file,
               pheno_file,
               covariate_file,
               mask_file,
               snp_indices,
               chunksize,
               n_randvecs,
               n_blocks,
               n_threads,
               rand_seed,
               log_level)

Rprof(NULL)
# Summarize the profiling results
summaryRprof(memory = "both")

