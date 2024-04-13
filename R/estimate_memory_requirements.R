#' Estimate the memory requirements
#'
#' Get an approximate estimation of the memory requirements
#' of the fame routine based on current input parameters
#'
#' @param n_samples Number of samples
#' @param n_snps Number of SNPs
#' @param n_blocks Number of genotype blocks
#' @param n_randvecs Number of random vectors
#'
#' @return Approximate memory requirement in giga bytes
#' @noRd
approximate_memory_requirements <- function(n_samples,
                                            n_snps,
                                            n_blocks,
                                            n_randvecs) {
    block_size <- ceiling(n_snps / n_blocks)
  genotype_blocks <- n_samples * block_size
  rand_vectors <- n_samples * n_randvecs
  trace_estimates <- 2 * n_samples * n_randvecs # additive and epistatic
  quadratic_forms <- n_samples * (3 + 9)
  block_stats <- 2 * block_size
  phenotype <- n_samples
  vc_estimate <- 2 * n_snps * 3 # point estimate and se for each component
  total <- genotype_blocks + rand_vectors + trace_estimates + quadratic_forms +
    block_stats + phenotype + vc_estimate
  return(total * 8 / 1024 / 1024 / 1024)
}
