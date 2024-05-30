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
