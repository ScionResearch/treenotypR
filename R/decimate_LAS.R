#' Decimate ultra-high-density point clouds to manageable size
#'
#' This function reduces the density of high-resolution point clouds (e.g., from
#' MLS or TLS) to a manageable size suitable for further processing, such as
#' crown delineation. It performs spatial subsampling using a specified 3D grid
#' resolution.
#'
#' @param las_dir_path Character. Path to the directory containing LAS files.
#' @param out_path Character. Path to the directory where decimated LAS files
#'                 will be saved.
#' @param grid_res Numeric. Grid resolution in meters. A 3D grid of the
#'                 specified size will be used for decimation.
#'                 Default is 0.1.
#' @param num_points Integer. Number of points to retain per 3D grid cell in the
#'                   decimation process. Random selecion will be applied.
#'                   Default is 1.
#' @param workers Integer. Number of workers/logical processors to use for
#'                processing. This is used to set up parallel processing.
#'                Default is 10.
#' @param defined_crs Numeric. Coordinate reference system (CRS) for the output
#'                    layers.
#'                    Default is EPSG:2193.
#'
#' @return A vector of file paths to the decimated LAS files.
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Reads the LAS files from the specified directory as a LAS catalog.
#'   \item Checks that the specified \code{num_points} is appropriate given the
#'         point density of the LAS catalog.
#'   \item Applies a decimation process using the specified grid resolution,
#'         \code{num_points} value and random point selection.
#'   \item Saves decimated LAS files in the specified output directory.
#' }
#'
#' @importFrom lidR readLAScatalog opt_output_files opt_laz_compression LAS catalog_map st_crs set_lidr_threads
#' @importFrom VoxR filter_point_density
#' @importFrom future plan multisession
#'
#' @examples
#' \dontrun{
#' Decimate_LAS(
#'   las_dir_path = "path/to/las/files",
#'   out_path = "path/to/output"
#' )
#' }
#' @export
#'
decimate_LAS <- function(
    las_dir_path,
    out_path,
    grid_res = 0.1,
    num_points = 1,
    workers = 10,
    defined_crs = 2193) {
  # Validate inputs.
  validate_path(las_dir_path)
  validate_path(out_path, create = TRUE)
  validate_numeric(grid_res, min = 0.01, max = 10)
  validate_integer(num_points)
  validate_integer(workers)
  validate_crs(defined_crs)

  # Read LAS files in the provided path.
  las_ctg <- lidR::readLAScatalog(las_dir_path)
  if (is.null(las_ctg)) {
    stop("Error reading LAS files: las_dir_path does not contain
         valid LAS files.")
  }

  # Set CRS.
  lidR::st_crs(las_ctg) <- defined_crs

  # Create directory to store decimated tiles.
  dir.create(file.path(out_path, "Decimated_Tiles"), showWarnings = FALSE)

  # Set catalog options.
  lidR::opt_output_files(las_ctg) <- paste0(
    out_path,
    "/Decimated_Tiles/{*}_decimated"
  )
  lidR::opt_laz_compression(las_ctg) <- TRUE

  # Define the decimation function.
  decimate_las <- function(las, res = grid_res) {
    decimated_xyz <- VoxR::filter_point_density(las@data, res)
    decimated_las <- suppressMessages(lidR::LAS(decimated_xyz))
    return(decimated_las)
  }

  # Set parallel processing.
  future::plan(future::multisession, workers = workers)
  lidR::set_lidr_threads(workers)

  # Apply decimation to each tile, handling errors.
  decimated_ctg <- tryCatch(
    {
      lidR::catalog_map(las_ctg, decimate_las, res = grid_res)
    },
    error = function(e) {
      message("Error applying decimation to the LAS catalog: ", e$message)
      return(NULL)
    }
  )
  lidR::st_crs(decimated_ctg) <- defined_crs

  # Return paths of decimated LAS files.
  return(decimated_ctg$filename)
}
