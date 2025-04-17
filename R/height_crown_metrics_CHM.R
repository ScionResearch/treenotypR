#' Calculate height and crown metrics from CHM
#'
#' This function calculates height metrics for each crown polygon, excluding
#' pixels below a specified minimum height threshold. Optionally, it calculates
#' the crown area for each polygon.
#'
#' @param chm_path Character. The file path to the canopy height model (CHM)
#'                 raster.
#' @param crowns_path Character. The file path to the crown polygons.
#' @param height Logical. If \code{TRUE}, calculates height metrics.
#'               Default is \code{TRUE}.
#' @param crown Logical. If \code{TRUE}, calculates crown area.
#'              Default is \code{TRUE}.
#' @param min_ht Numeric. The minimum height threshold (in meters) for
#'               filtering pixels in CHM. Pixels with heights below this will be
#'               excluded from metric calculations.
#'               Default is 2.
#' @param treeID Character. Name of the column in the crowns data that
#'               identifies individual trees.
#'               Default is "treeID".
#' @param defined_crs Integer. EPSG code for the Coordinate Reference System to
#'                    be assigned to both the CHM and crowns.
#'                    Default is EPSG:2193.
#'
#' @return A dataframe with calculated height metrics (if specified) and crown
#'         area (if specified) for each crown polygon.
#'
#' @details
#' \itemize{
#'   \item{Height metrics}{Calculates statistics including max, mean, standard
#'                         deviation, coefficient of variation, and specified
#'                         quantiles (5th, 50th, and 95th percentiles) for each
#'                         crown polygon, excluding pixels below the
#'                         \code{min_ht} threshold.}
#'   \item{Crown metrics}{Calculates the area of each crown polygon if
#'                        \code{crown = TRUE}.}
#' }
#'
#' @importFrom terra rast crs ifel extract
#' @importFrom sf st_read st_transform st_area st_centroid st_coordinates
#' @importFrom dplyr select mutate rename group_by summarize ungroup left_join rename row_number everything across
#' @importFrom tibble tibble rownames_to_column
#' @importFrom stats quantile sd
#' @importFrom utils globalVariables
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' crown_metrics <- get_tree_metrics_CHM(
#'   chm_path = "path/to/chm.tif",
#'   crowns_path = "path/to/crowns.shp",
#'   height = TRUE,
#'   crown = TRUE,
#'   min_ht = 2,
#'   treeID = "treeID",
#'   defined_crs = 2193
#' )
#' print(crown_metrics)
#' }
#'
#' @export
#'
#' @name get_tree_metrics_CHM
#'
utils::globalVariables(c("ID", "raster_values"))

get_tree_metrics_CHM <- function(chm_path,
                                 crowns_path,
                                 height = TRUE,
                                 crown = TRUE,
                                 min_ht = 2,
                                 treeID = "treeID",
                                 defined_crs = 2193) {
  # Validate inputs
  validate_path(chm_path)
  validate_path(crowns_path)
  validate_numeric(min_ht, min = 0.01, max = 100)
  validate_logical(height)
  validate_logical(crown)
  validate_crs(defined_crs)

  # Read and preprocess CHM raster
  chm <- terra::rast(chm_path)
  terra::crs(chm) <- paste0("EPSG:", defined_crs)
  masked_raster <- terra::ifel(chm < min_ht, NA, chm)

  # Read crown polygons and set CRS
  crowns <- sf::st_read(crowns_path, quiet = TRUE)
  crowns <- sf::st_transform(crowns, defined_crs)
  crowns_df <- crowns %>%
    tibble::tibble() %>%
    dplyr::select({{ treeID }}) %>%
    dplyr::mutate(ID = dplyr::row_number())

  # Initialize output tables
  metrics_ht <- crowns_df[, 1]
  metrics_crown <- crowns_df[, 1]

  # Calculate height metrics if enabled
  if (height) {
    extracted_values <- terra::extract(masked_raster, crowns) %>%
      tibble::tibble() %>%
      dplyr::rename(raster_values = !!tools::file_path_sans_ext(
        basename(chm_path)
      ))

    height_metrics <- extracted_values %>%
      dplyr::group_by(ID) %>%
      dplyr::summarize(
        zmax = max(raster_values, na.rm = TRUE),
        zmean = mean(raster_values, na.rm = TRUE),
        zsd = stats::sd(raster_values, na.rm = TRUE),
        zcv = stats::sd(raster_values, na.rm = TRUE) / mean(raster_values,
          na.rm = TRUE
        ),
        zq5 = stats::quantile(raster_values, 0.05, na.rm = TRUE),
        zq50 = stats::quantile(raster_values, 0.5, na.rm = TRUE),
        zq95 = stats::quantile(raster_values, 0.95, na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(crowns_df, by = "ID") %>%
      dplyr::select({{ treeID }}, dplyr::everything(), -ID)
  }

  # Calculate crown area and centroids if enabled
  if (crown) {
    crown_centroids <- sf::st_centroid(crowns)
    crown_coords <- sf::st_coordinates(crown_centroids)

    metrics_crown <- tibble::tibble(
      treeID = crowns[[treeID]],
      crown_area = as.numeric(sf::st_area(crowns)),
      crownX = crown_coords[, 1],
      crownY = crown_coords[, 2]
    )
  }

  # Combine results
  metrics <- crowns_df %>%
    dplyr::select(-ID) %>%
    dplyr::left_join(height_metrics, by = treeID) %>%
    dplyr::left_join(metrics_crown, by = treeID) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), round, 2))

  # Return combined metrics
  return(metrics)
}
