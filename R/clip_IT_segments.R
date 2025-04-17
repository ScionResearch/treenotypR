#' Clip individual tree segments from a point cloud using crown polygons
#'
#' This function reads tree crown polygons and then clips individual tree
#' segments based on the provided polygons. Each tree segment is saved as an
#' individual LAS file.
#'
#' @param las_dir_path Character. Path to the directory containing LAS files.
#' @param crowns_path Character. Path to the shapefile containing tree crown
#'                    polygons.Crowns must have a unique ID, see below.
#' @param out_path Character. Path to the directory where output files will be
#'                 saved.A new directory will be created here.
#' @param treeID Character. Name of the column in crowns file that uniquely
#'               identifies each tree.
#' @param workers Integer. Number of workers/logical processors to use for
#'                processing. This is used to set up parallel processing.
#'                Default is 10.
#' @param defined_crs Numeric. Coordinate reference system (CRS) for output
#'                    layers.
#'                    Default is EPSG:2193.
#'
#' @return Clipped individual tree segments will be saved in a newly created
#'         directory.A list containing paths to these segments will also be
#'         returned.
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Reads the tree crown polygons from the provided path.
#'   \item Reads the LAS files as a catalog.
#'   \item Creates an output directory to store the clipped tree segments.
#'   \item Clips individual tree segments from the LAS files, saving each tree
#'         segment as a separate file.
#'   \item Returns the file paths to the clipped LAS files, each representing an
#'         individual tree segment.
#' }
#'
#' @importFrom lidR readLAScatalog opt_output_files opt_laz_compression clip_roi st_crs
#' @importFrom sf st_read st_transform
#' @importFrom future plan multisession
#'
#' @examples
#' \dontrun{
#' clip_IT_segments(
#'   las_dir_path = "path/to/las/files",
#'   crowns_path = "path/to/crowns.shp",
#'   out_path = "path/to/output"
#' )
#' }
#'
#' @export
#'
clip_IT_segments <- function(
    las_dir_path,
    crowns_path,
    out_path,
    treeID = "treeID",
    defined_crs = 2193,
    workers = 10) {
  # Validate input.
  validate_path(las_dir_path)
  validate_path(crowns_path)
  validate_path(out_path, create = TRUE)
  validate_character(treeID)
  validate_crs(defined_crs)
  validate_integer(workers)

  # Create a directory to store individual tree segments.
  output_dir <- file.path(out_path, "Tree_Segments")
  counter <- 1
  while (dir.exists(output_dir)) {
    output_dir <- file.path(out_path, paste0("Tree_Segments_", counter))
    counter <- counter + 1
  }
  dir.create(output_dir)

  las_ctg <- tryCatch(
    {
      las_ctg <- lidR::readLAScatalog(las_dir_path)
      if (is.null(las_ctg)) stop("Error reading LAS catalog: las_path does
      not contain valid LAS files.")
      las_ctg
    },
    error = function(e) {
      message("Error reading LAS catalog: ", e$message)
      return(NULL)
    }
  )

  # Set output options.
  lidR::st_crs(las_ctg) <- defined_crs
  lidR::opt_output_files(las_ctg) <- paste0(output_dir, "/{", treeID, "}")
  lidR::opt_laz_compression(las_ctg) <- TRUE

  # Set up parallel processing and clip segments.
  future::plan(future::multisession, workers = workers)
  lidR::set_lidr_threads(workers)

  # Read and prep crowns filr.
  crowns <- sf::read_sf(crowns_path)
  crowns <- sf::st_transform(crowns, defined_crs)

  # Attempt to clip the LAS segments using tryCatch to handle potential errors.
  clipped_las_segments <- tryCatch(
    {
      lidR::clip_roi(las_ctg, crowns)
    },
    error = function(e) {
      message("Error clipping LAS catalog: ", e$message)
      return(NULL)
    }
  )

  # Return paths of clipped LAS files.
  return(clipped_las_segments$filename)
}
