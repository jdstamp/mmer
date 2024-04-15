#' Read the HDF5 file dataset
#'
#' Read the HDF5 file dataset
#'
#' @param filename File path to hdf5 file
#' @param datasetName name of the dataset to read
#'
#' @return Vector of indices
#' @export
read_h5_file <- function(filename, datasetName) {
    return(readH5File(filename, datasetName))
}
