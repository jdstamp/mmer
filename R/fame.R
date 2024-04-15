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
#' @importFrom tidyr pivot_longer
#' @export
fame <-
  function(plink_file,
           pheno_file,
           covariate_file,
           mask_file = NULL,
           n_randvecs = 10,
           gxg_indices = NULL,
           n_blocks = 100,
           rand_seed = -1) {
    logging::logReset()
    logging::basicConfig(level = "DEBUG")
    log <- logging::getLogger("fame")

    bim_file <- paste0(plink_file, ".bim")
    fam_file <- paste0(plink_file, ".fam")
    n_snps <- count_snps_bim(bim_file)
    n_samples <- count_samples(pheno_file)
    n_fam_lines <- count_fam(fam_file)

    if (n_samples != n_fam_lines) {
      stop("Number of samples in fam file and pheno file do not match.")
    }

    mem_req <- approximate_memory_requirements(n_samples,
                                               n_snps,
                                               n_blocks,
                                               n_randvecs)
    log$debug("Estimated memory requirement: %.2f GB per block.", mem_req)

    if (is.null(gxg_indices)) {
      gxg_indices <- c(1:n_snps)
    }

    result <-
      fame_cpp(
        plink_file,
        pheno_file,
        covariate_file,
        n_randvecs,
        n_blocks,
        rand_seed,
        gxg_indices - 1,
        # R is 1-indexed, C++ is 0-indexed
        mask_file
      )

    z_score <- abs(result$vc_estimate / result$vc_se)
    p_values <- 2 * (1 - pnorm(z_score))

    pve <- result$vc_estimate / apply(result$vc_estimate, 1, sum)
    id <- sprintf("gxg_%d", gxg_indices)
    vc_names <- c("id", "grm", "gxg", "error")
    component_col <- "component"
    summary <-
      data.frame(
        id = id,
        p = p_values[, 2],
        pve = pve[, 2],
        vc = result$vc_estimate[, 2],
        se = result$vc_se[, 2]
      )
    pve <- cbind(id, as.data.frame(pve))
    p_values <- cbind(id, as.data.frame(p_values))
    vc <- cbind(id, as.data.frame(result$vc_estimate))
    se <- cbind(id, as.data.frame(result$vc_se))
    colnames(pve) <- vc_names
    colnames(p_values) <- vc_names
    colnames(vc) <- vc_names
    colnames(se) <- vc_names
    result$p <- as_tibble(p_values)
    result$pve <- pivot_output(pve,
                               component_col,
                               "pve",
                               vc_names[2:4])
    result$vc_estimate <- pivot_output(vc,
                                       component_col,
                                       "vc_estimate",
                                       vc_names[2:4])
    result$vc_se <- pivot_output(se,
                                 component_col,
                                 "vc_se",
                                 vc_names[2:4])
    result$summary <- as_tibble(summary)
    return(result)
  }


pivot_output <- function(df, names_to, values_to, cols) {
  as_tibble(df) %>% pivot_longer(cols = cols,
                                 names_to = names_to,
                                 values_to = values_to)
}
