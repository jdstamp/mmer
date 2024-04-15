#' Simulate traits from plink genotypes
#'
#' @param plink_file Path to the plink data set without file extension
#' @param output_file Path to the file that stores the simulated trait
#' @param heritability Float between 0 and 1 that quantifies the heritability of the trait
#' @param rho Fraction of heritability due to additive only effects
#' @param n_causal Integers number of causal SNPs
#' @param gxg_indices Vector of SNP indices that are chosen to be epistatic
#'
#' @return None
#' @useDynLib famer
#' @export
#' @import genio
#' @import dplyr
#' @importFrom utils write.table
simulate_traits <- function(plink_file,
                            output_file,
                            heritability,
                            rho,
                            n_causal,
                            gxg_indices) {
  gxg <- get_groups(gxg_indices)
  sim <- simulate_traits_cpp(plink_file,
                             heritability,
                             rho,
                             n_causal,
                             gxg$group_1 - 1,
                             gxg$group_2 - 1
                             )
  fam_data <- read_fam(paste0(plink_file, ".fam"), verbose = FALSE)
  pheno_data <- data.frame(FID = fam_data$fam,
                           IID = fam_data$id,
                           TRAIT = sim$trait)
  write.table(
    pheno_data,
    file = output_file,
    sep = " ",
    quote = FALSE,
    row.names = FALSE
  )
}

