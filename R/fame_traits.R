#' Simulate traits from plink genotypes according to FAME paper
#'
#' @param plink_file Path to the plink data set without file extension
#' @param trait_output_file Path to the file that stores the simulated trait
#' @param additive_variance_component Float between 0 and 1 that quantifies the additive heritability of the trait
#' @param gxg_variance_component Variance component of the focal SNP
#' @param focal_snp Index of the focal SNP
#' @param gxg_indices Vector of SNP indices that are chosen to be epistatic with the focal SNP
#' @param additive_indices Vector of SNP indices that are chosen to have additive effects
#'
#' @return None
#' @useDynLib famer
#' @export
#' @import genio
#' @import dplyr
#' @importFrom utils write.table
fame_traits <- function(plink_file,
                        trait_output_file,
                        additive_variance_component,
                        gxg_variance_component,
                        focal_snp,
                        gxg_indices,
                        additive_indices) {
  sim <- fame_traits_cpp(
    plink_file,
    additive_variance_component,
    gxg_variance_component,
    focal_snp,
    gxg_indices - 1,
    additive_indices - 1
  )
  fam_data <- read_fam(paste0(plink_file, ".fam"), verbose = FALSE)
  pheno_data <- data.frame(FID = fam_data$fam,
                           IID = fam_data$id,
                           TRAIT = sim$trait)
  write.table(
    pheno_data,
    file = trait_output_file,
    sep = " ",
    quote = FALSE,
    row.names = FALSE
  )
}
