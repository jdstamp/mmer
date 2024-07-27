#' Estimate the memory requirements
#'
#' Get an approximate estimation of the memory requirements
#' of the mme routine based on current input parameters
#'
#' @param n_samples Number of samples
#' @param n_snps Number of SNPs
#' @param n_blocks Number of genotype blocks
#' @param n_randvecs Number of random vectors
#' @param chunksize Number of focal SNPS per chunk
#'
#' @return Approximate memory requirement in giga bytes
#' @noRd
approximate_memory_requirements <- function(n_samples,
                                            n_snps,
                                            n_blocks,
                                            n_randvecs,
                                            chunksize) {
  # memory of persisted XXz like intermediary results
  trace_estimates <- n_samples * n_randvecs * (chunksize + 1)
  quadratic_forms <- n_samples * 10 * (9 + 3)

  # dynamic memory of streaming blocks
  n_encoded <- ceiling(n_snps / n_blocks)
  genotype_blocks <- n_samples * n_encoded
  rand_vectors <- n_samples * n_randvecs
  block_stats <- 2 * n_encoded
  phenotype <- n_samples
  vc_estimate <- 2 * n_snps * 3 # point estimate and se for each component
  total <- genotype_blocks + rand_vectors + trace_estimates + quadratic_forms +
    block_stats + phenotype + vc_estimate
  return(total * 8 / 1024 / 1024 / 1024)
}
