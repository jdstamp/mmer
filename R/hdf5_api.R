#' Create an HDF5 File
#'
#' This function creates a new, empty HDF5 file at the specified location.
#'
#' @param hdf5_file A character string specifying the path and name of the HDF5 file to be created.
#'
#' @return No return value; the function creates the HDF5 file at the specified location.
#'
#' @examples
#' \dontrun{
#' # Create an empty HDF5 file
#' create_hdf5_file("example.h5")
#' }
#' 
#' @export
create_hdf5_file <- function(hdf5_file) {
  createH5File(hdf5_file)
}

#' Write Data to an HDF5 Dataset
#'
#' This function writes new data to an existing HDF5 file. If the dataset already exists,
#' it will be replaced with the new data.
#'
#' @param file_name A character string specifying the path to the HDF5 file.
#' @param dataset_name A character string specifying the name of the dataset to be written in the HDF5 file.
#' @param new_data The new data to write to the dataset. The data must be compatible with the dataset's structure.
#'
#' @return No return value; the function modifies the specified dataset in the HDF5 file.
#'
#' @examples
#' \dontrun{
#' # Write new data to a dataset named 'my_dataset' in the HDF5 file
#' write_hdf5_dataset("example.h5", "my_dataset", new_data)
#' }
#'
#' @export
write_hdf5_dataset <- function(file_name, dataset_name, new_data) {
  replaceH5Dataset(file_name, dataset_name, new_data)
}

#' Read Dataset from an HDF5 File
#'
#' This function reads a dataset from an existing HDF5 file.
#'
#' @param file_name A character string specifying the path to the HDF5 file.
#' @param dataset_name A character string specifying the name of the dataset within the HDF5 file to read.
#'
#' @return The content of the dataset from the HDF5 file, typically in the form of an R object.
#'
#' @examples
#' \dontrun{
#' # Read a dataset named 'my_dataset' from an HDF5 file
#' data <- read_hdf5_dataset("example.h5", "my_dataset")
#' }
#' @export
read_hdf5_dataset <- function(file_name, dataset_name) {
  return(readH5File(file_name, dataset_name))
}
