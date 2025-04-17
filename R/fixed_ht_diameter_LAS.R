#' Calculate diameter at a specified fixed height from a LAS object
#'
#' This function calculates the diameter of a tree stem at a specified fixed
#' height from an individual tree point cloud.
#'
#' @param stem_las A LAS object representing a tree stem point cloud.
#' @param fixed_ht Numeric. The fixed height (in meters) at which to calculate
#'                 the diameter.
#'                 Default is 1.4.
#' @param seg_ht Numeric. The vertical interval (in meters) above and below
#'               \code{fixed_ht} used to define the segment height range.
#'               Default is 0.2.
#' @param n_points Integer. Minimum number of points required in the segment
#'                 for diameter estimation.
#'                 Default is 25.
#' @param quantiles Numeric vector of length 2. Quantiles for filtering points
#'                  within the segment based on their distance to the shape
#'                  center. The values must be between 0 and 1.
#'                  Default is \code{c(0.1, 0.9)}.
#'
#' @return A tibble with columns: x, y, z, n_points, and diameter (in meters)
#'         representing the calculated values at the specified height.
#'
#' @details
#' This function calculates the diameter of a tree stem at a specified fixed
#' height by performing the following steps:
#' \itemize{
#'   \item Extracts a vertical segment of points from the point cloud within a
#'         range of \code{fixed_ht} ± \code{seg_ht}/2 meters.
#'   \item Filters points within this segment based on their distances from the
#'         estimated shape center, using quantiles specified in
#'         \code{quantiles} to remove outliers.
#'   \item Fits a shape to the filtered points to estimate the diameter at the
#'         specified height.
#'   \item Returns a dataframe containing the fitted shape's \code{x} and
#'         \code{y} center coordinates, the specified \code{z} height, the count
#'         of points \code{n_points} used, and the calculated \code{diameter}.
#' }
#'
#' @importFrom lidR filter_poi
#' @importFrom dplyr mutate select
#' @importFrom utils globalVariables
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' stem_points <- detect_stem(
#'   las_path = "path/to/your/file.laz",
#'   stem_only = TRUE
#' )
#' diameter_at_fixed_ht <- fixedht_diameter_LAS(stem_points,
#'   fixed_ht = 1.4,
#'   seg_ht = 0.25
#' )
#' }
#'
#' @export
#'
#' @name fixedht_diameter_LAS
#'
utils::globalVariables(c("%>%", "x", "y", "z", "diameter", "Z"))

fixedht_diameter_LAS <- function(stem_las,
                                 fixed_ht = 1.4,
                                 seg_ht = 0.2,
                                 n_points = 25,
                                 quantiles = c(0.1, 0.9)) {
  # Validate inputs.
  validate_las(stem_las, required_cols = c("X", "Y", "Z", "Stem"))
  validate_numeric(fixed_ht, min = 0.1, max = 100)
  validate_numeric(seg_ht, min = 0.01, max = 10)

  # Define the height range for the segment.
  seg_min <- fixed_ht - seg_ht
  seg_max <- fixed_ht + seg_ht

  # Extract points within the specified height range, handling the errors.
  segment <- tryCatch(
    {
      lidR::filter_poi(stem_las, Z > seg_min & Z < seg_max)
    },
    error = function(e) {
      message("Error in filtering segment point cloud: ", e$message)
      return(NULL)
    }
  )
  # Calculate diameter for the segment.
  profile <- segment_diameter(segment,
    n_points = n_points,
    quantiles = quantiles
  ) %>%
    dplyr::mutate(z = fixed_ht) %>%
    dplyr::select(x, y, z, diameter)

  return(profile)
}
