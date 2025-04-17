#' Calculate diameter at a specified fixed height from a diameter profile
#'
#' This function retrieves the estimated diameter of a tree stem at a specified
#' fixed height from a pre-computed diameter profile.
#'
#' @param diameter_profile A dataframe containing the diameter profile of a tree
#'                         with columns `z` (height) and `diameter`.
#' @param fixed_ht Numeric. The fixed height (in meters) at which to retrieve
#'                 the diameter.
#'                 Default is 1.4 meters.
#'
#' @return A dataframe with columns: `x`, `y`, `z`, `n_points`, and `diameter`
#'         (in meters), representing the closest diameter at the specified
#'         height.
#'
#' @details
#' This function identifies the closest diameter to a specified height by:
#' \itemize{
#'   \item Calculates distance between each \code{z} value in
#'         \code{diameter_profile} and the specified \code{fixed_ht} to
#'         determine proximity.
#'   \item Selects the smallest absolute difference, which represents the
#'         diameter measurement closest to \code{fixed_ht}.
#'   \item Returns a dataframe with the relevant details, including: `x`, `y`,
#'         `z`, `n_points`, and `diameter` (in meters).
#' }
#'
#' @importFrom dplyr mutate slice select
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
#' diameter_at_fixed_ht <- fixedht_diameter_DF(d_profile,
#'   fixed_ht = 1.4
#' )
#' }
#'
#' @export
#'
#' @name fixedht_diameter_DF
#'
utils::globalVariables(c("%>%", "z", "z_abs", "seg_id", "n_points"))

fixedht_diameter_DF <- function(diameter_profile,
                                fixed_ht = 1.4) {
  # Validate inputs.
  validate_dataframe(diameter_profile, required_cols = c("diameter", "z"))
  validate_numeric(fixed_ht, min = 0.1, max = 100)

  # Find the row closest to the specified fixed height.
  profile <- diameter_profile %>%
    dplyr::mutate(z_abs = abs(z - fixed_ht)) %>%
    dplyr::slice(which.min(z_abs)) %>%
    dplyr::select(-c(z_abs, seg_id, n_points))

  return(profile)
}
