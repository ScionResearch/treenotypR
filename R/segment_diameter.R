#' Get diameter from a stem segment point cloud
#'
#' This function takes a segment of LAS points as input and fits a cylinder
#' to the points, filtering them based on specified distance quantiles and
#' calculating the shape parameters (center coordinates and diameter).
#' It returns `x`, `y` coordinates at the center and `diameter` if the fitting
#' is successful.
#'
#' @param segment A LAS object representing a segment of stem points.
#' @param quantiles Numeric vector of length 2. Quantiles for point filtering
#'                  within segments based on distance from the shape center.
#'                  Values must be between 0 and 1.
#'                  Default is \code{c(0.1, 0.9)}.
#' @param n_points Integer. Minimum number of points required per segment for
#'                 diameter estimation.
#'                 Default is 25.
#'
#' @return A tibble with columns: `x`, `y`, and `diameter` (in meters) if
#'         fitting is successful, or a tibble with `NA` values if fitting fails.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Verifies if the segment contains the minimum required number of
#'         points as specified by \code{n_points}.
#'   \item Attempts to fit a shape to the points in the segment using the
#'         RANSAC algorithm via \code{TreeLS::shapeFit()}.
#'   \item Filters points based on their distances to the initial shape center
#'         within the quantile range defined by \code{quantiles}.
#'   \item Refits a shape to the filtered points.
#'   \item If the refit is successful, returns the shape parameters (center
#'         \code{x}, \code{y}, and \code{diameter}). Otherwise, returns
#'         \code{NA} values.
#' }
#'
#' @importFrom TreeLS shapeFit
#' @importFrom circular lsfit.circle
#' @importFrom tibble as_tibble_row tibble
#' @importFrom dplyr mutate filter select
#' @importFrom stats quantile
#' @importFrom utils globalVariables
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' segment <- readLAS("path/to/las_segment.laz")
#' diameter <- segment_diameter(segment,
#'   n_points = 30,
#'   quantiles = c(0.05, 0.95)
#' )
#' }
#'
#' @export
#'
#' @name segment_diameter
#'
utils::globalVariables(c("dist", "a", "b", "r"))

segment_diameter <- function(segment,
                             n_points = 25,
                             quantiles = c(0.1, 0.9)) {
  # Validate inputs.
  validate_las(segment, required_cols = c("X", "Y", "Z"))
  validate_integer(n_points)
  validate_vector(quantiles, 2, "quantiles")

  # Check if segment has enough points for analysis.
  seg_points <- nrow(segment@data)
  if (seg_points < n_points) {
    warning("Not enought points for shape fitting.")
    return(tibble::tibble(x = NA, y = NA, n_points = seg_points, diameter = NA))
  }

  # Try fitting a circle to the segment.
  pars <- tryCatch(
    TreeLS::shapeFit(segment, shape = "circle", algorithm = "ransac"),
    error = function(e) NULL
  )

  if (is.null(pars)) {
    # Return NA if shape fitting fails.
    warning("Shape fitting was unsuccessful.")
    return(tibble::tibble(x = NA, y = NA, n_points = seg_points, diameter = NA))
  }

  # Filter points based on distance from the circle center.
  segment_df <- tibble::tibble(segment@data) %>%
    dplyr::mutate(dist = sqrt((X - pars$X)^2 + (Y - pars$Y)^2)) %>%
    dplyr::filter(dplyr::between(
      dist, stats::quantile(dist, probs = quantiles[1]),
      stats::quantile(dist, probs = quantiles[2])
    )) %>%
    dplyr::select(X, Y)

  # Fit a circle to the quantile-filtered points.
  pars_qs <- tryCatch(
    circular::lsfit.circle(as.matrix(segment_df), units = "degree"),
    error = function(e) NULL
  )

  # Return circle parameters or NA values if fitting fails.
  if (!is.null(pars_qs)) {
    tibble::as_tibble_row(pars_qs$coefficients) %>%
      dplyr::select(x = a, y = b, diameter = r) %>%
      dplyr::mutate(n_points = nrow(segment_df), diameter = diameter * 200)
  } else {
    warning("Circle fitting was unsuccessful.")
    return(tibble::tibble(x = NA, y = NA, n_points = seg_points, diameter = NA))
  }
}
