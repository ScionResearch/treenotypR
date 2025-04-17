#' Create Density Raster and Calculate Point Density Statistics
#'
#' This function creates a density raster for a LAS catalog and calculates
#' point density statistics, providing insights on point distribution.
#'
#' @param las_path Character. Path to the directory containing LAS files.
#' @param raster_res Numeric. Resolution (in meters) of the density grid.
#'                     Default is 1 meter.
#' @param defined_crs Integer. EPSG code for the Coordinate Reference System.
#'                    Default is EPSG:2193.
#'
#' @return A list containing:
#' \item{raster}{The created density raster.}
#' \item{mean_density}{Mean density value across the raster.}
#' \item{sd_density}{Standard deviation of the density values.}
#'
#' @details
#' \itemize{
#'   \item Generates a raster of point densities (points per square meter).
#'   \item Calculates mean and standard deviation of point density values.
#' }
#'
#' @importFrom lidR readLAScatalog st_crs grid_density
#' @importFrom raster crs
#' @importFrom stats sd
#'
#' @examples
#' \dontrun{
#' density_stats <- create_density_raster(
#'   las_path = "path/to/las_catalog",
#'   out_path = "output/directory",
#'   raster_res = 1,
#'   defined_crs = 2193
#' )
#' print(density_stats)
#' }
#'
#' @export
#'
create_density_raster <- function(
    las_path,
    raster_res = 1,
    defined_crs = 2193) {
  # Validate inputs.
  validate_path(las_path)
  validate_numeric(raster_res, min = 0.01, max = 10)
  validate_crs(defined_crs)

  # Load LAS catalog with error handling.
  las_ctg <- tryCatch(
    {
      las_ctg <- lidR::readLAScatalog(las_path)
      if (is.null(las_ctg)) stop("Error reading LAS catalog: las_path does
      not contain valid LAS files.")
      las_ctg
    },
    error = function(e) {
      message("Error reading LAS catalog: ", e$message)
      return(NULL)
    }
  )

  # If LAS catalog loaded successfully, set CRS and proceed.
  if (!is.null(las_ctg)) {
    sf::st_crs(las_ctg) <- defined_crs

    # Create density raster with error handling.
    density_raster <- tryCatch(
      {
        dens_raster <- lidR::grid_density(las_ctg, res = raster_res)
        raster::crs(dens_raster) <- defined_crs
        dens_raster
      },
      error = function(e) {
        message("Error creating density raster: ", e$message)
        return(NULL)
      }
    )

    # Calculate density statistics if raster was successfully created.
    if (!is.null(density_raster)) {
      raster_values <- density_raster[]
      mean_density <- round(mean(raster_values, na.rm = TRUE))
      sd_density <- round(stats::sd(raster_values, na.rm = TRUE))

      # Return the raster and density statistics.
      return(list(
        raster = density_raster,
        mean_density = mean_density,
        sd_density = sd_density
      ))
    }
  }

  # Return NULL if there was an error and no results could be generated.
  return(NULL)
}
