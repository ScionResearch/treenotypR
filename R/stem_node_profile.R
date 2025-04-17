#' Detect nodes and internodes from a diameter profile
#'
#' This function analyzes the diameter profile of a tree stem to identify nodes
#' (areas of significant diameter increase) and internodes (areas of stability
#' or reduction in diameter). It provides detailed metrics for each identified
#' section, including height range and average diameter.
#'
#' @param diameter_profile A dataframe containing the diameter profile of a tree
#'                         with columns `z` (height) and `diameter` in meters.
#' @param tol Numeric. Tolerance threshold for identifying diameter changes in.
#'            centimeters. Sections with diameter changes larger than this
#'            threshold are classified as nodes, while sections with stable
#'            diameter are classified as internodes.
#'            Default is 1.
#' @param base Numeric. A base value used to inflate the initial segment
#'             diameter to assist in node detection. The base will be multiplied
#'             by tol and added to the diameter at the min_ht.
#'             Default is 3.
#' @param split A vector with a logical value and a character string. The
#'              character values must be either "Node" or "Internode". If the
#'              logical value is \code{TRUE}, only sections of type specified in
#'              the string (either "Node" or "Internode") are returned.
#'              Default is \code{c(FALSE, "Node")}.
#'
#' @return A dataframe containing columns: `section_id`, `z_min`, `z_max`,
#'         `z_avg`, `section_length`, `diameter_lower`, `diameter_upper`,
#'         `diameter_avg` and `section_type`. The z values and diameters will be
#'         in meters.
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Filters for continuous data, fills or drops missing values, and
#'         adjusts segment heights for consistent detection.
#'   \item Calculates the diameter difference for each segment and assigns a
#'         section type based on the tolerance threshold.
#'   \item Summarizes each section by height range, diameter range, average
#'         diameter, and section type.
#'   \item If specified in \code{split}, returns only sections of the specified
#'         type (either "Node" or "Internode").
#' }
#'
#' @importFrom dplyr filter arrange mutate select slice case_when first last pull bind_rows between group_by summarize
#' @importFrom tidyr fill drop_na
#' @importFrom purrr assign_in
#' @importFrom tibble rowid_to_column
#' @importFrom data.table rleid
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
#' n_profile <- node_profile(d_profile,
#'   tol = 1, base = 4,
#'   split = c(FALSE, "Node")
#' )
#' }
#'
#' @export
#'
#' @name stem_node_profile
#'
utils::globalVariables(c(
  "index_diff", "z_diff", "section_id", "diameter_lower",
  "diameter_upper", "diameter_avg", "node_id"
))

stem_node_profile <- function(diameter_profile,
                              tol = 1,
                              base = 4,
                              split = c(FALSE, "Node")) {
  # Validate inputs.
  validate_dataframe(diameter_profile,
    required_cols = c("z", "diameter", "seg_id")
  )
  validate_numeric(tol, min = 0.1, max = 50)
  validate_numeric(base, min = 1, max = 50)
  validate_vector(split, 2, "split")

  # Filter for valid diameter values and arrange by height.
  nan <- dplyr::filter(diameter_profile, !is.na(diameter)) %>%
    dplyr::arrange(z) %>%
    dplyr::mutate(
      index_diff = seg_id - dplyr::lag(seg_id),
      z_diff = z - dplyr::lag(z)
    )

  # Identify last height with continuous data.
  nan_z <- nan %>%
    dplyr::slice(which.min(index_diff)) %>%
    dplyr::pull(z)
  f_index <- which(diameter_profile$z == nan_z)

  # Segment height calculation.
  seg_ht <- nan %>%
    tidyr::drop_na() %>%
    dplyr::filter(index_diff == 1) %>%
    dplyr::slice(1) %>%
    dplyr::pull(z_diff)

  # Adjust base segment for node detection.
  base_d <- as.numeric(diameter_profile[f_index, "diameter"] + base * tol)
  base_z <- as.numeric(diameter_profile[f_index, "z"]) - seg_ht

  # Node detection.
  if (base_z <= 3) {
    diameter_profile <- diameter_profile %>%
      dplyr::slice(f_index) %>%
      purrr::assign_in(list(6, 1), base_d) %>%
      purrr::assign_in(list(4, 1), base_z) %>%
      purrr::assign_in(list(1, 1), 0) %>%
      dplyr::bind_rows(diameter_profile[f_index:nrow(diameter_profile), ])
  }

  # Calculate diameter changes and section type.
  node_metrics <- diameter_profile %>%
    dplyr::arrange(z) %>%
    dplyr::filter(z <= max(nan$z, na.rm = TRUE)) %>%
    dplyr::mutate(
      diameter_diff = c(0, diff(diameter)),
      section_type = dplyr::case_when(
        diameter_diff >= tol ~ "Node",
        dplyr::between(diameter_diff, -tol, tol) ~ "No change",
        diameter_diff <= -tol ~ "Internode"
      )
    ) %>%
    dplyr::filter(seg_id != 0) %>%
    tibble::rowid_to_column("section_id") %>%
    dplyr::group_by(section_id) %>%
    dplyr::summarize(
      z_min = min(z),
      z_max = max(z),
      diameter_lower = dplyr::first(diameter),
      diameter_upper = dplyr::last(diameter),
      diameter_avg = mean(diameter, na.rm = TRUE),
      section_type = dplyr::first(section_type),
      .groups = "drop"
    )

  # Fill "No change" sections by forward fill.
  node_metrics <- node_metrics %>%
    dplyr::mutate(section_type = ifelse(section_type == "No change",
      NA, section_type
    )) %>%
    tidyr::fill(section_type, .direction = "down") %>%
    dplyr::arrange(z_min) %>%
    dplyr::mutate(section_id = data.table::rleid(section_type)) %>%
    dplyr::group_by(section_id) %>%
    dplyr::summarize(
      z_min = min(z_min) - (seg_ht / 2),
      z_max = max(z_max) + (seg_ht / 2),
      z_avg = (z_min + z_max) / 2,
      section_length = z_max - z_min,
      diameter_lower = dplyr::first(diameter_lower),
      diameter_upper = dplyr::last(diameter_upper),
      diameter_avg = mean(diameter_avg, na.rm = TRUE),
      section_type = dplyr::first(section_type),
      .groups = "drop"
    ) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), round, 2))

  # Add unique IDs for nodes and internodes.
  node_ids <- node_metrics %>%
    dplyr::filter(section_type == "Node") %>%
    tibble::rowid_to_column("node_id") %>%
    dplyr::mutate(node_id = paste0("N", node_id))

  node_metrics <- node_metrics %>%
    dplyr::filter(section_type == "Internode") %>%
    tibble::rowid_to_column("node_id") %>%
    dplyr::mutate(node_id = paste0("IN", node_id)) %>%
    dplyr::bind_rows(node_ids) %>%
    dplyr::arrange(section_id)

  # Filter by section type if specified.
  if (as.logical(split[1])) {
    node_metrics <- node_metrics %>%
      dplyr::filter(section_type == split[2]) %>%
      tibble::rowid_to_column("section_id")
  }

  return(node_metrics)
}
