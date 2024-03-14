#' Fast Marginal Epistasis Test implementation
#'
#' @param plink_file File path for plink genotype dataset without *.bed extension.
#' @param pheno_file File path for plink *.pheno data.
#' @param annotation_file File path for plink *.pheno data.
#' @param covariate_file File path for covariates data.
#' @param gxg_bin Integer.
#' @param num_evec Integer. Number of random vectors.
#' @param snp_index Integer.
#' @param jack_number Integer.
#'
#' @return A list of P values and PVEs
#' @name fame
#' @useDynLib famer
#' @import Rcpp
#' @import RcppEigen
#' @export
fame <- function(plink_file, pheno_file, annotation_file, covariate_file, gxg_bin, num_evec, snp_index, jack_number) {
  fame_cpp(plink_file, pheno_file, annotation_file, covariate_file, gxg_bin, num_evec, snp_index, jack_number)
  return(1)
}
