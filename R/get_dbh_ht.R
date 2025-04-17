#' Detect the nearest internode height where DBH must be measured
#' (PlotSafe Guideline Replication)
#'
#' This function identifies the nearest internode in a stem profile where DBH
#' must be measured according to the PlotSafe Guidelines. Typically breast
#' height is 1.4 meters for direct measurements. If no direct measurement is
#' possible, this function estimates the nearest height to 1.4 meters where no
#' nodes/swellings are present.
#'
#' @param node_profile A dataframe containing the node profile with columns
#'                     `z_min`, `z_max`, and `section_type`, where
#'                     `section_type` denotes "Node" or "Internode".
#' @param bh Numeric. The target height in meters to find the closest internode.
#'           Default is 1.4 meters.
#' @param min_ht Numeric. The minimum height in meters to consider for internode
#'               detection.
#'               Default is 0.2 meters.
#' @param max_ht Numeric. The maximum height in meters to consider for internode
#'               detection.
#'               Default is 2.9 meters.
#' @param tol1 Numeric. Tolerance range (in meters) around the target height for
#'             a direct internode measurement.
#'             Default is 0.1.
#' @param tol2 Numeric. Additional tolerance range (in meters) for
#'             finding overlapping internode segments if no direct measurement
#'             is possible and split is required.
#'             Default is 0.3.
#'
#' @return A tibble with columns:`dbh_class`, which is the classification of the
#'        DBH measurement type ("DirectM", "DirectSS", "DirectSS", or "Split"),
#         `nearest_int1` which is the nearest internode height or range as a
#'        list, and `nearest_int2`, an additional internode height or range if a
#'        split measurement is needed.
#' \itemize{
#'   \item DirectM: Diameter is measured at breast height (1.4 meters) directly
#'         over an internode.
#'   \item Direct SS: Diameter is measured directly at a distance of `tol1`
#'         above or below breast height, over an internode (1.3 to 1.5 m).
#'   \item Split: To avoid a long welling, diameter measurements are taken at
#'         two locations (internodes) over the internodes and then averaged.
#'   \item DirectLS: When diameter cannot be measured using the above methods
#'         due to the presence of a long node or swollen area, a direct
#'         measurement is taken between 0.2 and 2.9 meters over an internode.
#' }
#'
#' @details
#' This function detects the nearest internode height to a specified target
#' height, typically used to find the closest internode around breast height
#' (1.4 meters). The function operates as follows:
#' \itemize{
#'   \item  Attempts to identify an internode that overlaps with the specified
#'          height (`bh`) within a small tolerance (`tol1`).
#'   \item If no direct measurement is available, the function searches for a
#'         pair of internode sections, one above and one below the target
#'         height, to calculate an estimated height.
#'   \item If neither a direct nor split measurement is possible, the function
#'         finds the closest available internode below the target height.
#' }
#'
#' @importFrom dplyr select mutate if_else across rowwise ungroup slice arrange filter
#' @importFrom tibble tibble rowid_to_column
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
#' n_profile <- stem_node_profile(d_profile,
#'   tol = 1, base = 4,
#'   split = c(FALSE, "Node")
#' )
#' dbh_ht <- variable_DBH_ht(n_profile,
#'   bh = 1.4, min_ht = 0.2,
#'   max_ht = 2.9, tol1 = 0.1, tol2 = 0.3
#' )
#' }
#'
#' @export
#'
#' @name get_dbh_ht
#'
utils::globalVariables(c(
  "%>%", "z_min", "z_max", "section_type", "z_min_abs",
  "z_max_abs", "overlap_min_max", "overlap_14",
  "z_lower_abs", "z_higher_abs", "overlap_13_15",
  "dbh_class", "nearest_ht", "sect", "id", "overlap"
))

get_dbh_ht <- function(node_profile,
                       bh = 1.4,
                       min_ht = 0.2,
                       max_ht = 2.9,
                       tol1 = 0.1,
                       tol2 = 0.3) {
  # Validate inputs.
  validate_dataframe(node_profile,
    required_cols = c("z_min", "z_max", "section_type")
  )
  validate_numeric(bh, min = 0.2, max = 2.9)
  validate_numeric(min_ht, min = 0.1, max = 10)
  validate_numeric(max_ht, min = 0.2, max = 10)
  validate_numeric(tol1, min = 0.01, max = 0.5)
  validate_numeric(tol2, min = 0.01, max = 1)

  # Prep data.
  node_profile <- node_profile %>%
    dplyr::select(z_min, z_max, section_type) %>%
    dplyr::mutate(
      abs_min_ht = min_ht >= z_min & min_ht <= z_max,
      abs_max_ht = max_ht >= z_max & max_ht <= z_max,
      overlap_14 = bh >= z_min & bh <= z_max,
      z_min_abs = abs(bh - z_min),
      z_max_abs = abs(bh - z_max),
      z_lower_abs = dplyr::if_else(
        z_min_abs < z_max_abs, z_min_abs, z_max_abs
      ),
      z_higher_abs = dplyr::if_else(
        z_min_abs < z_max_abs, z_max_abs, z_min_abs
      ),
      nearest_ht = dplyr::if_else(z_min_abs > z_max_abs, z_max, z_min),
      sect = dplyr::if_else(z_min >= bh, "upper", "lower"),
      dplyr::across(dplyr::where(is.numeric), round, 4)
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      overlap_min_max = max(z_min, min_ht) <= min(z_max, max_ht),
      overlap_13_15 = max(z_min, (bh - tol1)) <= min(z_max, (bh + tol1))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(overlap_min_max == TRUE, section_type == "Internode") %>%
    dplyr::mutate(
      z_min = dplyr::if_else(z_min < min_ht, min_ht, z_min),
      z_max = dplyr::if_else(z_max > max_ht, max_ht, z_max),
      dbh_class = dplyr::if_else(any(TRUE %in% overlap_14), "DirectM", "Other")
    ) %>%
    dplyr::select(
      z_min, z_max, z_min_abs, z_max_abs, z_lower_abs, z_higher_abs,
      overlap_14, overlap_13_15, dbh_class, nearest_ht, sect
    )

  # Initialize result tibble.
  nrst_int_ht <- tibble::tibble()

  while (nrow(nrst_int_ht) != 1) {
    # Direct measurement detection around 1.4 m.
    direct_m <- node_profile %>%
      dplyr::filter(dplyr::if_any(
        c(overlap_14, overlap_13_15),
        ~ . == TRUE
      )) %>%
      dplyr::slice(which.min(z_lower_abs)) %>%
      dplyr::mutate(dbh_class = dplyr::if_else(dbh_class == "Other",
        "DirectSS", dbh_class
      ))

    if (nrow(direct_m) >= 1) {
      nearest_int1 <- direct_m %>%
        dplyr::select(z_min, z_max) %>%
        unlist(use.names = FALSE)
      nrst_int_ht <- tibble::tibble(
        dbh_class = direct_m$dbh_class,
        nearest_int1 = list(nearest_int1),
        nearest_int2 = NA
      )
      break
    }

    # Estimate split heights if no direct measurement is available.
    node_profile_asc <- node_profile %>%
      dplyr::arrange(z_lower_abs) %>%
      tibble::rowid_to_column("id")
    for (i in seq_len(nrow(node_profile_asc))) {
      split_ht1 <- dplyr::filter(node_profile_asc, id == i)
      range_lower <- as.numeric(split_ht1$z_lower_abs)
      range_higher <- as.numeric(split_ht1$z_higher_abs)
      split_sect1 <- as.character(split_ht1$sect)

      # Find overlapping range for nearest split height.
      node_profile_asc <- dplyr::filter(
        node_profile_asc, id > i,
        sect != split_sect1
      )
      split_ht2 <- node_profile_asc %>%
        dplyr::rowwise() %>%
        dplyr::mutate(overlap = min(z_higher_abs, range_higher) - max(
          z_lower_abs, range_lower
        )) %>%
        dplyr::ungroup() %>%
        dplyr::filter(overlap >= (0 - tol2)) %>%
        dplyr::slice(which.max(overlap))
      if (nrow(split_ht2) == 0) next

      split_ht1 <- split_ht1 %>%
        dplyr::select(z_min, z_max) %>%
        unlist(use.names = FALSE)
      split_ht2 <- split_ht2 %>%
        dplyr::select(z_min, z_max) %>%
        unlist(use.names = FALSE)
      nrst_int_ht <- tibble::tibble(
        dbh_class = "Split",
        nearest_int1 = list(split_ht1),
        nearest_int2 = list(split_ht2)
      )
    }
    if (nrow(nrst_int_ht) != 0) break

    # Fallback for a direct measurement below target height.
    direct_ls <- node_profile %>%
      dplyr::filter(nearest_ht <= 1.4) %>%
      dplyr::slice(which.min(z_lower_abs)) %>%
      dplyr::select(z_min, z_max) %>%
      unlist(use.names = FALSE)
    if (length(direct_ls) != 0) {
      nrst_int_ht <- tibble::tibble(
        dbh_class = "DirectLS",
        nearest_int1 = list(direct_ls),
        nearest_int2 = NA
      )
    } else {
      nrst_int_ht <- tibble::tibble(
        dbh_class = "DirectLS",
        nearest_int1 = list(c(1.3, 1.3)),
        nearest_int2 = NA
      )
    }
  }

  # Add tolerances.
  nrst_int_ht <- nrst_int_ht %>%
    dplyr::mutate(
      breast_ht = bh,
      directM_tol = ifelse(dbh_class == "DirectSS", tol1, NA),
      splitM_tol = ifelse(dbh_class == "Split", tol2, NA)
    )

  return(nrst_int_ht)
}
