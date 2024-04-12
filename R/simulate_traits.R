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
                             gxg$group_1,
                             gxg$group_2
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

#' Split vector into groups
#'
#' @param gxg_indices vector on SNP indices
#' @return named list with indeces for the two groups
#' @export
get_groups <- function(gxg_indices) {
  mid <- length(gxg_indices) %/% 2
  gxg_group_1 <- gxg_indices[1:mid]
  gxg_group_2 <- gxg_indices[(mid+1):length(gxg_indices)]
  g <- list(group_1 = gxg_group_1,
            group_2 = gxg_group_2)
  return(g)
}

#' Update mask
#'
#' @param mask_file Genotype mask
#' @param gxg_group_1 snps in group 1
#' @param gxg_group_2 snps in group 2
#' @export
update_mask <- function(mask_file, gxg_group_1, gxg_group_2) {
  # update the mask file
  for (gt in gxg_group_1) {
    h5_ds <- sprintf("gxg/%d", gt)
    m <- readH5File(mask_file, h5_ds)
    m <- sample(m, size = length(m))
    new_m <- unique(c(gxg_group_2, m))
    new_m <- sort(new_m[1:length(m)])
    replaceH5Dataset(mask_file, h5_ds, new_m)
  }
}

