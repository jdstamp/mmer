#' Fast Marginal Epistasis Test implementation
#'
#' @param plink_file File path for plink genotype dataset without *.bed extension.
#' @param pheno_file File path for plink *.pheno data.
#' @param covariate_file File path for covariates data.
#' @param n_randvecs Integer. Number of random vectors.
#' @param variant_indices List of indices for the focal variants.
#' @param n_blocks Integer representing the number of blocks the SNPs will be read in.
#' @param rand_seed Integer to seed generation of random vectors. Only positive values are considered.
#'
#' @return A list of P values and PVEs
#' @name fame
#' @useDynLib famer
#' @import Rcpp
#' @import RcppEigen
#' @export
fame <- function(plink_file, pheno_file, covariate_file, n_randvecs,
variant_indices = NULL, n_blocks, rand_seed = -1) {
  result <- fame_cpp(plink_file, pheno_file, covariate_file, n_randvecs,
  n_blocks, rand_seed, variant_indices)
  return(result)
}
