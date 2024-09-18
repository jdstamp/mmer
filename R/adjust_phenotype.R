#' Adjust Phenotype by Covariates
#' 
#' This function adjusts the phenotype by the covariates.
#' @param pheno_file File path to phenotype.
#' @param covariate_file File path to covariates.
#' 
#' @export
adjust_phenotype <- function(pheno_file, covariate_file) {
  residuals <- get_residuals(pheno_file, covariate_file)
  # Read the phenotype file
  header <- readLines(pheno_file, n = 1)
  header <- strsplit(header, "\\s+")[[1]]
  pheno_data <- read.table(pheno_file, header = TRUE)
  pheno_data <- cbind(pheno_data[, 1:2], residuals)
  colnames(pheno_data) <- header
  return(pheno_data)
}
