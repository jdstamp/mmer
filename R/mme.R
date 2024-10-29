#' Fast Marginal Epistasis Test implementation
#'
#' @param plink_file File path for plink genotype dataset without *.bed extension.
#' @param pheno_file File path for plink *.pheno data.
#' @param mask_file File path for gxg mask data.
#' @param gxg_indices List of indices for the focal variants.
#' @param chunksize Integer representing the number of SNPs analyzed with the same set of random vectors.
#' @param n_randvecs Integer. Number of random vectors.
#' @param n_blocks Integer representing the number of blocks the SNPs will be read in.
#' @param n_threads Integer representing the number of threads that are setup for OMP.
#' @param gxg_h5_group String of the hdf5 group with the gxg mask. These SNPs will be included.
#' @param ld_h5_group String of the hdf5 group with the ld mask. These SNPs will be excluded.
#' @param rand_seed Integer to seed generation of random vectors. Only positive values are considered.
#' @param log_level Log level.
#'
#' @return A list of P values and PVEs
#' @name mme
#' @useDynLib mmer
#' @import Rcpp
#' @import RcppEigen
#' @import dplyr
#' @importFrom stats pnorm
#' @importFrom tidyr pivot_longer
#' @importFrom progress progress_bar
#' @export
mme <-
  function(plink_file,
           pheno_file,
           mask_file = NULL,
           gxg_indices = NULL,
           chunksize = NULL,
           n_randvecs = 10,
           n_blocks = 100,
           n_threads = 1,
           gxg_h5_group = "gxg",
           ld_h5_group = "ld",
           rand_seed = -1,
           log_level = "WARNING") {
    logging::logReset()
    logging::basicConfig(level = log_level)
    log <- logging::getLogger("mme")

    n_gxg_indices <- length(gxg_indices)
    log$debug("Number of gxg indices: %d", n_gxg_indices)

    bim_file <- paste0(plink_file, ".bim")
    fam_file <- paste0(plink_file, ".fam")
    n_snps <- count_snps_bim(bim_file)
    n_samples <- count_samples(pheno_file)
    n_fam_lines <- count_fam(fam_file)

    log$debug("Dataset: %s", plink_file)
    log$debug("Trait file: %s", pheno_file)
    log$debug("Mask file: %s", mask_file)
    log$debug("Number of samples: %d", n_samples)
    log$debug("Number of SNPs: %d", n_snps)
    log$debug("Number of random vectors: %d", n_randvecs)
    log$debug("Number of blocks: %d", n_blocks)

    if (check_openmp()) {
      log$info("openMP is enabled")
      log$info("Number of requested threads: %d", n_threads)
    }

    if (n_samples != n_fam_lines) {
      stop("Number of samples in fam file and pheno file do not match.")
    }

    mem_req <- approximate_memory_requirements(n_samples, n_snps, n_blocks, n_randvecs, chunksize)
    log$debug("Estimated memory requirement: %.2f GB", mem_req)

    if (is.null(gxg_indices)) {
      gxg_indices <- c(1:n_snps)
    }

    if (is.null(chunksize)) {
      n_chunks <- ceiling(n_gxg_indices / n_threads)
      log$debug("No chunksize specified. Using %d chunks.", n_chunks)
    } else {
      n_chunks <- ceiling(n_gxg_indices / chunksize)
      log$debug("Chunksize set to %d. Using %d chunks.", chunksize, n_chunks)
    }

    if (n_chunks > 1) {
      chunks <- split(gxg_indices,
                      cut(seq_along(gxg_indices), n_chunks, labels = FALSE))
    } else {
      chunks <- list(`1` = gxg_indices)
    }

    VC <- NULL
    SE <- NULL
    TIME <- NULL

    pb <- progress_bar$new(
      format = "processing chunks [:bar] :percent elapsed: :elapsed eta: :eta",
      total = n_chunks,
      clear = FALSE,
      width = 60
    )
    pb$tick(0)
    for (i in seq_along(chunks)) {
      chunk <- chunks[[i]]
      result <-
        mme_cpp(
          plink_file,
          pheno_file,
          mask_file,
          n_randvecs,
          n_blocks,
          rand_seed,
          chunk - 1,
          # R is 1-indexed, C++ is 0-indexed
          n_threads,
          gxg_h5_group,
          ld_h5_group
        )
      VC <- rbind(VC, result$vc_estimate)
      SE <- rbind(SE, result$vc_se)
      TIME <- c(TIME, result$duration)
      pb$tick()
    }
    total_duration <- sum(TIME)
    average_duration <- total_duration / n_gxg_indices
    log$debug("Total computation time: %f seconds",
                  total_duration)
    log$debug("Average computation time per SNP: %f seconds",
              average_duration)

    z_score <- (VC / SE)
    p_values <- (1 - pnorm(z_score)) # one sided test

    pve <- VC / apply(VC, 1, sum)

    bim <- read.delim(bim_file, header = FALSE)
    colnames(bim) <- c("chromosome", "variant_id", "genetic_distance", "position", "allele1", "allele2")

    id <- bim$variant_id[gxg_indices]
    chromosome <- bim$chromosome[gxg_indices]
    position <- bim$position[gxg_indices]
    vc_names <- c("id", "grm", "gxg", "error")
    component_col <- "component"
    summary <-
      data.frame(
        id = id,
        index = gxg_indices,
        chromosome = chromosome,
        position = position,
        p = p_values[, 2],
        pve = pve[, 2],
        vc = VC[, 2],
        se = SE[, 2]
      )
    pve <- cbind(id, as.data.frame(pve))
    p_values <- cbind(id, as.data.frame(p_values))
    vc <- cbind(id, as.data.frame(VC))
    se <- cbind(id, as.data.frame(SE))
    colnames(pve) <- vc_names
    colnames(p_values) <- vc_names
    colnames(vc) <- vc_names
    colnames(se) <- vc_names
    result$p <- as_tibble(p_values)
    result$pve <- pivot_output(pve, component_col, "pve", vc_names[2:4])
    result$vc_estimate <- pivot_output(vc, component_col, "vc_estimate", vc_names[2:4])
    result$vc_se <- pivot_output(se, component_col, "vc_se", vc_names[2:4])
    result$summary <- as_tibble(summary)
    result$average_duration <- average_duration
    return(result)
  }


pivot_output <- function(df, names_to, values_to, cols) {
  as_tibble(df) %>% pivot_longer(cols = all_of(cols),
                                 names_to = names_to,
                                 values_to = values_to)
}
