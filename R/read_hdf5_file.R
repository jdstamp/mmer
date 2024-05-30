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

#' Write the HDF5 file dataset
#'
#' Write the HDF5 file dataset
#'
#' @param filename File path to hdf5 file
#' @param dataset_name name of the dataset to read
#' @param new_data Vector to write to data set
#'
#' @return Writes to file
#' @export
write_h5_dataset <- function(filename, dataset_name, new_data) {
  replaceH5Dataset(filename, dataset_name, new_data)
}
