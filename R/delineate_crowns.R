#' Delineate tree crowns on the canopy height model (CHM)
#'
#' This function delineates individual tree crowns from a canopy height model
#' (CHM) raster using a watershed segmentation approach. The resulting polygons
#' represent crown boundaries and are matched to the provided stem map.
#'
#' @param chm_path Character. Path to the CHM raster file.
#' @param stem_map_path Character. Path to the stem map (point shapefile).
#' @param min_ht Numeric. Minimum height threshold in meters for crown
#'               delineation. Any pixels having height below this value will be
#'               excluded from the crown delineation process.
#' @param smooth Logical. If \code{TRUE}, smooths the edges of delineated
#'               crowns.
#'               Default is \code{TRUE}.
#' @param defined_crs Numeric. Coordinate reference system (CRS) for output
#'                    layers.
#'                    Default is EPSG:2193.
#'
#' @return A polygon object containing the delineated crown polygons with their
#'         matched stem IDs.
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Reads and processes the CHM raster, applying focal filters to fill
#'         gaps and smooth raster).
#'   \item Uses watershed segmentation to delineate tree crowns using stem
#'         map locations of trees as seed points.
#'   \item Label the delineated crowns with matching stem ID.
#'   \item Optionally (if specified) smooths the edges of the crown polygons.
#' }
#'
#' @importFrom raster raster crs
#' @importFrom terra focal
#' @importFrom sf st_read st_transform as_Spatial st_as_sf
#' @importFrom dplyr mutate select left_join rename
#' @importFrom tidyr drop_na
#' @importFrom smoothr fill_holes smooth
#' @importFrom ForestTools mcws
#' @importFrom tibble tibble
#' @importFrom utils globalVariables
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' segs_sf <- delineate_crowns(
#'   chm_path = "path/to/chm.tif",
#'   stem_map_path = "path/to/stem_map.shp"
#' )
#' ggplot2::ggplot(data = crowns) +
#'   ggplot2::geom_sf() +
#'   ggplot2::geom_sf_text(ggplot2::aes(label = treeID),
#'     size = 3, color = "blue"
#'   ) +
#'   ggplot2::theme_bw() +
#'   ggplot2::labs(title = "Crowns with tree IDs")
#' }
#'
#' @export
#'
#' @name delineate_crowns
#'
utils::globalVariables(c("%>%", ":="))

delineate_crowns <- function(chm_path,
                             stem_map_path,
                             min_ht = 2,
                             treeID = "treeID",
                             smooth = TRUE,
                             defined_crs = 2193) {
  # Validate inputs.
  validate_path(chm_path)
  validate_path(stem_map_path)
  validate_character(treeID)
  validate_numeric(min_ht, min = 0, max = 100)
  validate_crs(defined_crs)
  validate_logical(smooth)

  # Set moving window for focal operations.
  ws <- matrix(1, 3, 3)

  # Function for filling NoData.
  fill.na <- function(x, i = 5) {
    if (is.na(x)[i]) {
      return(mean(x, na.rm = TRUE))
    } else {
      return(x[i])
    }
  }

  # Read and preprocess CHM raster.
  chm <- raster::raster(chm_path)
  chm <- terra::focal(chm, ws, fun = fill.na)
  chm <- terra::focal(chm, w = ws, fun = mean, na.rm = TRUE)

  # Set CRS for CHM raster.
  raster::crs(chm) <- defined_crs

  # Read stem map and match CRS.
  stem_map <- sf::st_read(stem_map_path, quiet = TRUE)
  stem_map <- sf::st_transform(stem_map, defined_crs)

  # Check if the specified column exists in stem_map and rename if necessary.
  if (treeID %in% colnames(stem_map)) {
    stem_map <- stem_map %>%
      dplyr::rename_with(~ paste0(treeID, ".original"), .cols = treeID)
  }

  # Add a new 'treeID' column with row numbers as unique identifiers.
  stem_map <- stem_map %>%
    dplyr::mutate(!!treeID := dplyr::row_number())
  stem_map_df <- tibble::tibble(stem_map) %>%
    dplyr::select(-c("geometry"))

  # Convert stem map to spatial object for mcws function.
  stem_map_sp <- sf::as_Spatial(stem_map)

  # Delineate crowns using watershed segmentation.
  segs <- tryCatch(
    {
      ForestTools::mcws(stem_map_sp, chm,
        minHeight = min_ht,
        format = "POLYGON"
      )
    },
    error = function(e) {
      stop("Error in watershed segmentation: ", e$message)
    }
  )

  segs_sf <- sf::st_as_sf(segs)

  # Join crown polygons with stem map IDs and filter out NA values.
  segs_sf <- segs_sf %>%
    dplyr::left_join(stem_map_df, by = "treeID") %>%
    dplyr::select(-treeID) %>%
    dplyr::rename(!!treeID := !!dplyr::sym(paste0(treeID, ".original"))) %>%
    tidyr::drop_na()

  # Fill holes in the polygons.
  segs_sf <- smoothr::fill_holes(segs_sf, threshold = 5)

  # Smooth edges if enabled
  if (smooth) {
    segs_sf <- smoothr::smooth(segs_sf, method = "ksmooth", smoothness = 2L)
  }

  # Ensure the final output has the correct CRS.
  segs_sf <- sf::st_transform(segs_sf, defined_crs)

  return(segs_sf)
}
