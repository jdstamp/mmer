#' Simulate a trait from plink genotypes
#'
#' @param plink_file Path to the plink data set without file extension
#' @param output_file Path to the file that stores the simulated trait
#' @param additive_heritability Float between 0 and 1 that quantifies the heritability of the trait
#' @param gxg_heritability Fraction of heritability due to additive only effects
#' @param additive_indices Vector of SNP indices that are chosen to contribute to the additive heritability
#' @param gxg_indices_1 Vector of SNP indices that are chosen to be epistatic in group 1
#' @param gxg_indices_2 Vector of SNP indices that are chosen to be epistatic in group 2
#' @param log_level Log level.
#'
#' @return None
#' @useDynLib mmer
#' @export
#' @import genio
#' @import dplyr
#' @importFrom utils write.table
simulate_a_trait <- function(plink_file,
                            output_file,
                            additive_heritability,
                            gxg_heritability,
                            additive_indices,
                            gxg_indices_1,
                            gxg_indices_2,
                            log_level = "WARNING") {

   if (additive_heritability + gxg_heritability > 1) {
     stop("Additive heritability and gxg heritability should sum to less than 1")
   } else if (additive_heritability < 0 || gxg_heritability < 0) {
     stop("Heritabilities should be positive")
   }
  logging::basicConfig(level = log_level)
  log <- logging::getLogger("mmer::simulate_a_trait")
  sim <- simulate_a_trait_cpp(plink_file,
                             additive_heritability,
                             gxg_heritability,
                             additive_indices - 1,
                             gxg_indices_1 - 1,
                             gxg_indices_2 - 1
                             )
  log$info("Simulated traits with additive variance %.2f, gxg variance %.2f, and error variance %.2f", sim$additive_variance, sim$gxg_variance, sim$error_variance)
  fam_data <- read_fam(paste0(plink_file, ".fam"), verbose = FALSE)
  pheno_data <- data.frame(FID = fam_data$fam,
                           IID = fam_data$id,
                           TRAIT = sim$trait - mean(sim$trait)
                           )
  write.table(
    pheno_data,
    file = output_file,
    sep = " ",
    quote = FALSE,
    row.names = FALSE
  )
}

