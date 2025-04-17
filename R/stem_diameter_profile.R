#' Calculate stem diameter profile
#'
#' This function calculates the diameter profile of a tree stem from an
#' individual tree point cloud.
#'
#' @param stem_las A LAS object that has stem points detected and classified.
#' @param seg_ht Numeric. Height of each segment (profile interval) in meters.
#' @param min_ht Numeric. Minimum stem height for segmentation in meters.
#'               Points below this height will be excluded. Default is 0.2.
#' @param max_ht Numeric or character. Maximum height for segmentation,
#'               either as "zmax" (to use the maximum height in the LAS file) or
#'               a specific numeric value. Default is \code{"zmax"}.
#' @param min_diameter Numeric. Minimum allowable stem diameter in meters.
#'                     Defaults to 0.05.
#' @param max_diameter Numeric. Maximum allowable stem diameter in meters.
#'                     Defaults to 0.6.
#' @param quantiles Numeric vector of length 2. Quantiles for point filtering
#'                  within each segment, used to exclude outliers in diameter
#'                  calculation.
#'                  Default is \code{c(0.1, 0.9)}.
#' @param n_points Integer. Minimum number of points per segment required for
#'                 diameter estimation. Default is 25.
#'
#' @return A dataframe with columns: `ID`, `x`, `y`, `z`, `points`, and
#'         `diameter`(in meters), representing the diameter profile of the stem.
#'
#' @details
#' This function calculates a diameter profile for a tree stem by performing the
#' following steps:
#' \itemize{
#'   \item Divides the point cloud into vertical segments, each with a height
#'         specified by \code{seg_ht}. Only points within the height
#'         range specified by \code{min_ht} and \code{max_ht} are used.
#'   \item Fits a shape to the filtered points in each segment and calculates
#'         the segment diameter. Only segments with at least \code{n_points}
#'         are considered for diameter estimation.
#'   \item Returns a data frame where each row represents a segment, including
#'         the estimated `x` and `y` center coordinates,`z` height (in meters),
#'         point count, and calculated `diameter` (in meters).
#' }
#'
#' @importFrom dplyr mutate filter select bind_rows arrange across
#' @importFrom lidR filter_poi
#' @importFrom tibble rowid_to_column tibble
#' @importFrom utils globalVariables
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' stem_points <- detect_stem(
#'   las_path = "path/to/your/file.laz",
#'   stem_only = TRUE
#' )
#' d_profile <- stem_diameter_profile(stem_points,
#'   seg_ht = 0.25,
#'   min_ht = 0.3,
#'   max_ht = "zmax",
#'   quantiles = c(0.2, 0.85),
#'   n_points = 25,
#'   max_diameter = 0.5,
#'   min_diameter = 0.1
#' )
#' }
#'
#' @export
#'
#' @name stem_diameter_profile
#'
utils::globalVariables("diameter_diff")

stem_diameter_profile <- function(stem_las,
                                  seg_ht = 0.2,
                                  min_ht = 0.3,
                                  max_ht = "zmax",
                                  min_diameter = 0.05,
                                  max_diameter = 0.6,
                                  quantiles = c(0.1, 0.9),
                                  n_points = 25) {
  # Validate inputs.
  validate_las(stem_las, required_cols = c("X", "Y", "Z", "Stem"))
  validate_numeric(seg_ht, min = 0.01, max = 5)
  validate_numeric(min_ht, min = 0.1, max = 100)
  validate_numeric(min_diameter, min = 0.02, max = 5)
  validate_numeric(max_diameter, min = 0.02, max = 5)

  # Validate diameter range.
  if (min_diameter >= max_diameter) {
    stop("Invalid diameter range: min_diameter must be less than max_diameter.")
  }

  # Set up parameters for z-axis segmentation.
  seg_n <- round((1.4 - min(stem_las$Z)) / seg_ht, 0)
  z_min <- round(1.4 - (seg_n * seg_ht), 4) - (seg_ht / 2)
  zmax <- ifelse(max_ht != "zmax", as.numeric(max_ht),
    round(max(stem_las$Z), 4)
  )

  validate_numeric(zmax, min = 0.1, max = 100)

  # Validate height range.
  if (min_ht >= zmax) {
    stop("Invalid Input: height range: min_ht must be less than max_ht.")
  }

  # Initialize an empty tibble for the diameter profile results.
  diameter_profile_df <- tibble::tibble()

  # Loop over segments to calculate diameter profiles.
  while (z_min < zmax) {
    z_max <- z_min + seg_ht
    segment <- lidR::filter_poi(stem_las, Z > z_min & Z < z_max)

    # Calculate diameter profile for the segment.
    profile <- segment_diameter(segment,
      n_points = n_points,
      quantiles = quantiles
    )

    # Append profile with average z-coordinate.
    diameter_profile_df <- profile %>%
      dplyr::mutate(z = (z_min + z_max) / 2, .after = y) %>%
      dplyr::bind_rows(diameter_profile_df)

    # Increment to next segment.
    z_min <- z_max
  }

  # Post-process: filter diameters by range and smooth out anomalies.
  diameter_profile_df <- diameter_profile_df %>%
    dplyr::mutate(
      diameter = ifelse(diameter >= (max_diameter * 100) | diameter <= (
        min_diameter * 100), NA, diameter),
      diameter_diff = abs(diameter - dplyr::lag(diameter)) / dplyr::lag(
        diameter
      ),
      diameter = ifelse(diameter_diff >= 0.5, NA, diameter),
      dplyr::across(dplyr::where(is.numeric), round, 2)
    ) %>%
    dplyr::select(x, y, z, n_points, diameter) %>%
    dplyr::arrange(z) %>%
    tibble::rowid_to_column("seg_id")

  return(diameter_profile_df)
}
