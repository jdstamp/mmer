#' Update mask
#'
#' @param mask_file Genotype mask
#' @param gxg_group_1 snps in group 1
#' @param gxg_group_2 snps in group 2
#' @export
update_mask <- function(mask_file, gxg_group_1, gxg_group_2) {
  # update the mask file
  for (gt in gxg_group_1) {
    h5_ds <- sprintf("gxg/%d", gt - 1)
    m <- readH5File(mask_file, h5_ds)
    m <- sample(m, size = length(m))
    new_m <- unique(c(gxg_group_2, m))
    new_m <- sort(new_m[1:length(m)])
    replaceH5Dataset(mask_file, h5_ds, new_m)
  }
}
