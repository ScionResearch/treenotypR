#' Calculate height and crown metrics from an individual tree point cloud
#'
#' This function calculates height and crown metrics from a point cloud.
#' This is specifically for extracting metrics from airborne laser scanning
#' (ALS) data, including unmanned aerial vehicle laser scanning (ULS).
#' Users can toggle height and crown metric calculations separately and apply
#' a minimum height filter to exclude points below a specified height.
#'
#' @param las_path Character. The file path to the LAS file.
#' @param height Logical. If \code{TRUE}, calculates height metrics.
#'               Default is \code{TRUE}.
#' @param crown Logical. If \code{TRUE}, calculates crown metrics.
#'              Default is \code{TRUE}.
#' @param min_ht Numeric. The minimum height (in meters) to filter points.
#'               Points with heights below this will be excluded.
#'               Default is 0.3.
#' @param defined_crs Integer. Coordinate Reference System (CRS) to assign to
#'                    the LAS point cloud.
#'                    Defaults is EPSG:2193.
#'
#' @return A dataframe with calculated height metrics (if specified) and crown
#'         metrics (if specified).
#'
#' @importFrom tibble tibble
#' @importFrom dplyr select mutate bind_cols across relocate
#' @importFrom lidR readLAS filter_poi st_concave_hull st_crs cloud_metrics
#' @importFrom sf st_convex_hull st_centroid st_coordinates st_area
#' @importFrom geometry convhulln
#' @importFrom stats quantile
#' @importFrom tools file_path_sans_ext
#' @importFrom utils globalVariables
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' las <- readLAS("path/to/lasfile.las")
#' metrics_ALS <- get_tree_metrics(las_path,
#'     min_ht = 2
#'     height = TRUE,
#'     crown = TRUE)
#' print(metrics_ALS)
#' }
#'
#' @export
#'
#' @name get_tree_metrics_ALS
#'
utils::globalVariables(c("zq99", "zq95", "file_id"))

get_tree_metrics_ALS <- function(las_path,
                                 min_ht = 0.3,
                                 height = TRUE,
                                 crown = TRUE,
                                 defined_crs = 2193) {
  # Validate inputs.
  validate_path(las_path)
  validate_numeric(min_ht, min = 0.01, max = 100)
  validate_logical(height)
  validate_logical(crown)
  validate_crs(defined_crs)

  # Initialize empty tibble for results.
  results <- tibble::tibble()

  # Try reading the LAS file with error handling.
  las <- tryCatch(
    {
      las <- lidR::readLAS(las_path)
      if (is.null(las)) stop("LAS file could not be read or is empty.")
      las
    },
    error = function(e) {
      message("Error reading LAS file: ", e$message)
      return(NULL)
    }
  )

  # Proceed to filtering and crs setting if LAS file was successfully read.
  if (!is.null(las)) {
    las <- tryCatch(
      {
        las <- lidR::filter_poi(las, Z >= min_ht)
        lidR::projection(las) <- defined_crs
        las
      },
      error = function(e) {
        message("Error filtering LAS points or setting CRS: ", e$message)
        return(NULL)
      }
    )
  }

  # Calculate height metrics if enabled.
  if (height) {
    results <- tryCatch(
      {
        # Calculate cloud metrics and additional quantile.
        point_metrics <- lidR::cloud_metrics(las, ~ stdmetrics(
          X, Y, Z, Intensity, ReturnNumber, Classification,
          dz = 1, th = min_ht, zmin = 0
        )) %>%
          data.frame() %>%
          dplyr::mutate(
            zq99 = stats::quantile(las$Z, probs = 0.9999, na.rm = TRUE),
            dplyr::across(dplyr::where(is.numeric), round, 2)
          ) %>%
          dplyr::relocate(zq99, .after = zq95)
      },
      error = function(e) {
        message("Error calculating height metrics: ", e$message)
        return(NULL)
      }
    )
  }

  # Calculate crown metrics if enabled.
  if (crown) {
    # Extract unique 3D and 2D coordinates from the LiDAR point cloud data.
    xyz <- unique(as.matrix(las@data[, c("X", "Y", "Z")]))
    xy <- unique(as.matrix(las@data[, c("X", "Y")]))

    # Try creating the 2D concave hull.
    concavehull_2D <- tryCatch(
      {
        lidR::st_concave_hull(las)
      },
      error = function(e) {
        message("Error creating 2D concave hull: ", e$message)
        return(NULL)
      }
    )

    # Try creating the 2D convex hull
    convexhull_2D <- tryCatch(
      {
        geometry::convhulln(xy, "FA")
      },
      error = function(e) {
        message("Error creating 2D convex hull: ", e$message)
        return(NULL)
      }
    )

    # Try creating the 3D convex hull if 2D was successful.
    convexhull_3D <- tryCatch(
      {
        geometry::convhulln(xyz, "FA")
      },
      error = function(e) {
        message("Error creating 3D convex hull: ", e$message)
        return(NULL)
      }
    )

    # Calculate crown metrics
    crown_metrics <- tibble::tibble(
      crownX = as.numeric(sf::st_coordinates(
        sf::st_centroid(concavehull_2D)
      )[1]),
      crownY = as.numeric(sf::st_coordinates(
        sf::st_centroid(concavehull_2D)
      )[2]),
      concavehull2D_area = as.numeric(sf::st_area(concavehull_2D)),
      convexhull2D_area = convexhull_2D$vol,
      convexhull3D_volume = convexhull_3D$vol,
      convexhull3D_area = convexhull_3D$area
    ) %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), round, 2))

    # Bind crown metrics to the existing results dataframe.
    results <- dplyr::bind_cols(results, crown_metrics)
  }

  # Add file basename to results.
  results <- results %>%
    dplyr::mutate(file_id = tools::file_path_sans_ext(
      base::basename(las_path)
    )) %>%
    dplyr::select(file_id, dplyr::everything())


  return(results)
}
