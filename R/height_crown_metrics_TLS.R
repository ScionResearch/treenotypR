#' Calculate height and crown metrics from an individual tree point cloud
#'
#' This function calculates height and crown metrics from a point cloud.
#' This is specifically for extracting metrics from mobile laser scanning (MLS)
#' or terrestrial laser scanning (TLS) data.
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
#' @details
#' This function provides the following metrics:
#' \itemize{
#'   \item Height Metrics: If \code{height = TRUE}, the function calculates:
#'     \itemize{
#'       \item \code{maxHt}: Maximum height of the point cloud.
#'       \item \code{meanHt}: Mean height.
#'       \item \code{modeHt}: Median height.
#'       \item \code{stHt}: Standard deviation of heights.
#'       \item \code{CVHt}: Coefficient of variation of heights.
#'       \item Height percentiles: 5\%, 10\%, 20\%, 25\%, 30\%, 40\%, 50\%,
#'             60\%, 70\%, 75\%, 80\%, 90\%, 95\%, and 99\%.
#'     }
#'
#'   \item Crown Metrics: If \code{crown = TRUE}, the function calculates:
#'     \itemize{
#'       \item \code{crownX}, \code{crownY}: Centroid location of the crown.
#'       \item \code{concavehull2D_area}: Area of the 2D concave hull.
#'       \item \code{convexhull2D_area}: Area of the 2D convex hull.
#'       \item \code{convexhull3D_volume}: Volume of the 3D convex hull.
#'       \item \code{convexhull3D_area}: Surface area of the 3D convex hull.
#'     }
#'     These metrics offer insights into canopy spread and crown volume.
#' }
#'
#' @importFrom tibble tibble
#' @importFrom dplyr select mutate bind_cols across
#' @importFrom tidyr drop_na
#' @importFrom lidR readLAS filter_poi st_concave_hull projection
#' @importFrom sf st_convex_hull st_centroid st_coordinates st_area
#' @importFrom geometry convhulln
#' @importFrom stats quantile
#' @importFrom tools file_path_sans_ext
#' @importFrom moments skewness kurtosis
#' @importFrom utils globalVariables
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' path <- readLAS("path/to/lasfile.las")
#' metrics <- get_tree_metrics_TLS(
#'   las_path = path,
#'   height = TRUE,
#'   crown = TRUE,
#'   min_ht = 2
#' )
#' print(metrics)
#' }
#'
#' @export
#'
#' @name get_tree_metrics_TLS
#'
utils::globalVariables(c("percentHt", "file_id"))

get_tree_metrics_TLS <- function(las_path,
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
      las <- lidR::readLAS(las_path, select = "xyz")
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
        lidR::st_crs(las) <- defined_crs
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
    data <- tibble::tibble(las@data$Z) %>%
      tidyr::drop_na() %>%
      unlist(use.names = FALSE)

    height_metrics <- tibble::tibble(
      zmax = max(data),
      zmean = mean(data),
      zsd = stats::sd(data),
      zcv = stats::sd(data) / mean(data),
      zskew = moments::skewness(data),
      zkurt = moments::kurtosis(data),
      percentHt = list(stats::quantile(data,
        probs = c(
          0.05, 0.1, 0.2, 0.25, 0.3, 0.4,
          0.5, 0.6, 0.7, 0.75, 0.8, 0.9,
          0.95, 0.99
        )
      ))
    )

    # Process quantiles into individual columns.
    quantile_names <- paste0(
      "zq",
      c(
        5, 10, 20, 25, 30, 40, 50,
        60, 70, 75, 80, 90, 95, 99
      )
    )
    quantiles <- dplyr::as_tibble(stats::setNames(
      as.list(height_metrics$percentHt[[1]]),
      quantile_names
    ))
    results <- height_metrics %>%
      dplyr::select(-percentHt) %>%
      dplyr::bind_cols(quantiles) %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), round, 2))
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

    # Only proceed if the concave hull was successfully created.
    if (!is.null(concavehull_2D)) {
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
    }

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
      basename(las_path)
    )) %>%
    dplyr::select(file_id, dplyr::everything())

  return(results)
}
