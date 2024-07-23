library(famer)

plink_file <- "/Users/jds/data/ukbb/c12_100k-samples_010k-snps_imputed"
pheno_file <- "/Users/jds/data/ukbb/c12_100k-samples_010k-snps_imputed.phen"
# plink_file <- "/users/jstamp1/data/jstamp1/ukbb/c12_100k-samples_010k-snps_imputed"
# pheno_file <- "/users/jstamp1/data/jstamp1/ukbb/c12_100k-samples_010k-snps_imputed.phen"
covariate_file <- ""
mask_file <- ""
# mask_file <- "/Users/jds/data/ukbb/c12_100k-samples_010k-snps_masksize500.h5"
log_level <- "DEBUG"
n_blocks <- 100
rand_seed <- 123

snp_indices <- 1:90

chunksize <- c(30)
n_threads <- c(10)
n_randvecs <- c(100)

parameter_grid <- expand.grid(chunksize = chunksize,
                              n_threads = n_threads,
                              n_randvecs = n_randvecs)

duration <- parameter_grid
duration$average_duration <- rep(0, nrow(parameter_grid))

# Rprof(memory.profiling = TRUE)

for (i in 1:nrow(parameter_grid)) {
  result <- fame(plink_file,
                 pheno_file,
                 covariate_file,
                 mask_file,
                 snp_indices,
                 parameter_grid$chunksize[i],
                 parameter_grid$n_randvecs[i],
                 n_blocks,
                 parameter_grid$n_threads[i],
                 rand_seed,
                 log_level)
  duration$average_duration[i] <- result$average_duration
}

print(duration)

# Rprof(NULL)
# # Summarize the profiling results
# summaryRprof(memory = "both")

