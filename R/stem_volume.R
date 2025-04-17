#' Calculate cumulative volume of a tree stem
#'
#' This function calculates the cumulative volume of a tree stem based on its
#' diameter profile, using the formula for the volume of a cylinder for each
#' segment. Missing sectional diameter values can be optionally filled or
#' dropped.
#'
#' @param diameter_profile A dataframe containing the diameter profile of a tree
#'                         with columns `z` (height) and `diameter`.
#' @param drop.NoData Logical. If \code{TRUE}, rows with missing diameter values
#'                    are omitted.
#'                    Default is \code{FALSE}.
#' @param fill.NoData A vector with a logical value and a character string. The
#'                    logical value indicates whether to fill missing diameters,
#'                    and the string specifies the fill direction ("up", "down",
#'                    or "updown"). `up` uses the preceding value, `down` uses
#'                    the next value, and `updown` fills as needed from both
#'                    directions.
#'                    Default is \code{c(TRUE, "updown")}.
#'
#' @return A dataframe with columns `vol_min_ht`, where the minimum height of
#'         the tree for which volume was calculated, `vol_max_ht`, where the
#'         maximum height of the tree for which volume was calculated, and
#'         `volume`, the cumulative volume (in cubic meters) of the tree stem
#'         within the height range.
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Checks that `diameter_profile` contains necessary columns and
#'         validates input parameters.
#'   \item Identifies the last valid height with data continuity and removes
#'         rows above this height.
#'   \item Determines the height of each segment based on consecutive `z`
#'         values.
#'   \item Fills or drops missing diameter values according to the user input.
#'   \item Calculates the volume for each segment and then computes the
#'         cumulative volume (in cubic meters) up to each height (in meters).
#' }
#'
#' @importFrom dplyr filter arrange mutate select slice first pull
#' @importFrom tidyr fill drop_na replace_na
#' @importFrom utils globalVariables
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' stem_points <- detect_stem(
#'   las_path = "path/to/your/file.laz",
#'   stem_only = TRUE
#' )
#' d_profile <- stem_diameter_profile(stem_points)
#' tree_volume <- stem_volume(d_profile)
#' }
#'
#' @export
#'
#' @name stem_volume
#'
utils::globalVariables(c("volume", "cum_vol", "vol_min_ht"))

stem_volume <- function(diameter_profile,
                        drop.NoData = FALSE,
                        fill.NoData = c(TRUE, "updown")) {
  # Validate inputs.
  validate_dataframe(diameter_profile, required_cols = c("z", "diameter"))
  validate_logical(drop.NoData)
  validate_vector(fill.NoData, 2, "fill.NoData")

  # Remove heights above the last continuous data point.
  nan <- dplyr::filter(diameter_profile, !is.na(diameter)) %>%
    dplyr::arrange(z) %>%
    dplyr::mutate(
      index_diff = seg_id - dplyr::lag(seg_id),
      z_diff = z - dplyr::lag(z)
    )

  nan_z <- nan %>%
    dplyr::slice(which.max(z)) %>%
    dplyr::filter(index_diff == 1) %>%
    dplyr::pull(z)

  diameter_profile <- dplyr::filter(diameter_profile, z <= nan_z)

  # Calculate segment height.
  seg_ht <- nan %>%
    tidyr::drop_na() %>%
    dplyr::filter(index_diff == 1) %>%
    dplyr::slice(1) %>%
    dplyr::pull(z_diff)

  # Fill gaps if specified.
  if (as.logical(fill.NoData[1])) {
    diameter_profile <- diameter_profile %>%
      tidyr::fill(diameter, .direction = fill.NoData[2])
  }

  # Drop NA values if specified.
  if (drop.NoData) {
    diameter_profile <- tidyr::drop_na(diameter_profile)
  }

  # Calculate cumulative volume.
  tree_volume <- diameter_profile %>%
    dplyr::mutate(
      volume = round(pi * (diameter / 200)^2 * as.numeric(seg_ht), 2),
      cum_vol = cumsum(tidyr::replace_na(volume, 0)),
      vol_min_ht = min(z, na.rm = TRUE)
    ) %>%
    dplyr::slice(which.max(cum_vol)) %>%
    dplyr::select(vol_min_ht,
      vol_max_ht = z,
      volume = cum_vol
    )

  return(tree_volume)
}
