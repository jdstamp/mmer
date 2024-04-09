#' Fast Marginal Epistasis Test implementation
#'
#' @param plink_file File path for plink genotype dataset without *.bed extension.
#' @param pheno_file File path for plink *.pheno data.
#' @param covariate_file File path for covariates data.
#' @param mask_file File path for gxg mask data.
#' @param n_randvecs Integer. Number of random vectors.
#' @param gxg_indices List of indices for the focal variants.
#' @param n_blocks Integer representing the number of blocks the SNPs will be read in.
#' @param rand_seed Integer to seed generation of random vectors. Only positive values are considered.
#'
#' @return A list of P values and PVEs
#' @name fame
#' @useDynLib famer
#' @import Rcpp
#' @import RcppEigen
#' @import dplyr
#' @importFrom stats pnorm
#' @export
fame <- function(plink_file, pheno_file, covariate_file, mask_file = NULL,
n_randvecs = 10, gxg_indices = NULL, n_blocks = 100, rand_seed = -1) {
  logging::logReset()
  logging::basicConfig(level = "DEBUG")
  log <- logging::getLogger("fame")

  bim_file <- paste0(plink_file, ".bim")
  fam_file <- paste0(plink_file, ".fam")
  n_snps <- count_snps_bim(bim_file)
  n_samples <- count_samples(pheno_file)
  n_fam_lines <- count_fam(fam_file)

  if(n_samples != n_fam_lines) {
    stop("Number of samples in fam file and pheno file do not match.")
  }

  mem_req <- n_samples * ceiling(n_snps / n_blocks) * 8 / 1024 / 1024 / 1024
  log$debug("Estimated memory requirement: %.2f GB per block.", mem_req)

  if(is.null(gxg_indices)) {
    gxg_indices <- c(1:n_snps)
  }

  result <- fame_cpp(plink_file, pheno_file, covariate_file, n_randvecs,
  n_blocks, rand_seed, gxg_indices - 1, mask_file) # R is 1-indexed, C++ is 0-indexed
  z_score <- result$Est / result$SE
  p_values <- (1 - pnorm(z_score))
  pve <- result$Est / apply(result$Est, 1, sum)
  id <- sprintf("gxg_%d", gxg_indices)
  summary <- data.frame(id = id, p = p_values[,2], pve = pve[,2])
  pve <- cbind(id, pve)
  colnames(pve) <- c("id", "grm", "gxg", "error")
  result$p <- p_values
  result$pve <- as_tibble(pve)
  result$summary <- as_tibble(summary)
  return(result)
}
