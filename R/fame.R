#' Fast Marginal Epistasis Test implementation
#'
#' @param plink_file File path for plink genotype dataset without *.bed extension.
#' @param pheno_file File path for plink *.pheno data.
#' @param covariate_file File path for covariates data.
#' @param n_randvecs Integer. Number of random vectors.
#' @param focal_snp_index Integer  of the focal SNP.
#' @param n_blocks Integer representing the number of blocks the SNPs will be read in.
#'
#' @return A list of P values and PVEs
#' @name fame
#' @useDynLib famer
#' @import Rcpp
#' @import RcppEigen
#' @export
fame <- function(plink_file, pheno_file, covariate_file, n_randvecs, focal_snp_index, n_blocks) {
  result <- fame_cpp(plink_file, pheno_file, covariate_file, n_randvecs, focal_snp_index, n_blocks)
  return(result)
}
