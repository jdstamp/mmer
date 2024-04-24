#' Update mask
#'
#' Update the mask for two groups of snps.
#'
#' @param mask_file Genotype mask
#' @param gxg_group_1 snps in group 1 in 1-based indexing
#' @param gxg_group_2 snps in group 2 in 1-based indexing
#' @export
update_mask <- function(mask_file, gxg_group_1, gxg_group_2) {
  # convert to C++ 0-based indexing
  gxg_group_1 <- gxg_group_1 - 1
  gxg_group_2 <- gxg_group_2 - 1
  update_single_group(mask_file, gxg_group_1, gxg_group_2)
  update_single_group(mask_file, gxg_group_2, gxg_group_1)
}


#' Single group update
#'
#' Update the mask for a single group of snps.
#'
#' @param mask_file Genotype mask
#' @param gxg_group_1 snps in group 1
#' @param gxg_group_2 snps in group 2
#' @noRd
update_single_group <- function(mask_file, gxg_group_1, gxg_group_2) {
  for (gt in gxg_group_1) {
    h5_ds <- sprintf("gxg/%d", gt)
    m <- readH5File(mask_file, h5_ds)
    m <- sample(m, size = length(m)) # shuffle the mask for randomized update
    new_m <- unique(c(gxg_group_2, m))
    new_m <- sort(new_m[1:length(m)])
    replaceH5Dataset(mask_file, h5_ds, new_m)
  }
}
