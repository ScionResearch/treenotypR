#' Estimate diameter at variable height
#'
#' This function estimates the diameter of a tree stem at variable heights from
#' a diameter profile, using on pre-computed internode information.
#'
#' @param diameter_profile A dataframe containing the diameter profile of a tree
#'                         with columns `z` (height) and `diameter`.
#' @param var_ht A dataframe generated from `variable_DBH_ht` with columns
#'               specifying breast height (`breast_ht`), tolerance
#'               (`splitM_tol`), internode ranges, and DBH classification.
#'
#' @return A dataframe with columns: `dbh_class`, `dbh_height`, and `diameter`,
#'         representing the estimated diameter (in meters) based on height
#'         criteria. If the required height measurements or internode data are
#'         unavailable, the function returns \code{NA} values in the output
#'         dataframe.
#'
#' @details
#' This function retrieves the estimated diameter at a height defined from a
#' specified diameter profile using internode ranges and t olerances:
#' \itemize{
#'   \item Calculates internode segment heights and nearest internode ranges.
#'   \item Determines the nearest height segments that overlap with the breast
#'         height (bh).
#'   \item Based on the `dbh_class`, selects the appropriate diameter measure:
#'         \code{"DirectM"} retrieves the diameter directly at breast height
#'         (1.4m), \code{"DirectSS"} finds the diameter closest to the nearest
#'         internode, and \code{"Split"} averages diameters within a split
#'         tolerance range.
#' }
#'
#' @importFrom dplyr filter arrange mutate select slice summarise bind_rows pull if_else across tibble
#' @importFrom tidyr drop_na
#' @importFrom purrr discard
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
#' n_profile <- stem_node_profile(d_profile)
#' dbh_ht <- get_dbh_ht(n_profile)
#' dbh <- variable_ht_dbh(d_profile, dbh_ht)
#' print(dbh)
#' }
#'
#' @export
#'
#' @name variable_ht_dbh
#'
globalVariables(c("z_abs_int1", "z_abs_split_int1", "z_abs_split_int2"))

variable_ht_dbh <- function(
    diameter_profile,
    var_ht) {
  # Validate inputs.
  validate_dataframe(diameter_profile, required_cols = c("z", "diameter"))
  validate_dataframe(var_ht, required_cols = c(
    "breast_ht", "splitM_tol",
    "dbh_class", "nearest_int1", "nearest_int2"
  ))

  # Retrieve necessary values from var_ht.
  bh <- round(var_ht$breast_ht, 2)
  tol2 <- round(var_ht$splitM_tol, 2)
  dbh_class <- as.character(var_ht$dbh_class)

  # Skip to last if dbh_class is DirectM.
  if (dbh_class == "DirectM") {
    nearest_int1 <- 0
  } else {
    # Segment height calculation.
    seg_ht <- diameter_profile %>%
      dplyr::filter(!is.na(diameter)) %>%
      dplyr::arrange(z) %>%
      dplyr::mutate(
        index_diff = dplyr::row_number() - dplyr::lag(dplyr::row_number()),
        z_diff = z - dplyr::lag(z)
      ) %>%
      tidyr::drop_na() %>%
      dplyr::filter(index_diff == 1) %>%
      dplyr::slice(1) %>%
      dplyr::pull(z_diff)

    # Calculate internode ranges.
    int1 <- sort(unlist(var_ht$nearest_int1) + c(seg_ht / 2, -seg_ht / 2))
    int2 <- sort(unlist(var_ht$nearest_int2) + c(seg_ht / 2, -seg_ht / 2))

    # Calculate nearest internode height.
    d_int1 <- abs(int1 - bh)
    d_int2 <- abs(int2 - bh)
    nearest_int1 <- ifelse(d_int1[1] <= d_int1[2], int1[1], int1[2])

    # Determine first and second split heights based on distances.
    split_int1 <- ifelse(d_int1[1] >= d_int2[1], int1[1], int2[1])
    split_int1_sect <- ifelse(split_int1 %in% int1, "int2", "int1")

    if (get(split_int1_sect)[1] == get(split_int1_sect)[2]) {
      split_int2 <- get(split_int1_sect)[1]
    } else {
      # Calculate ideal second split height as distance from split_int1
      p_ht <- 2.8 - split_int1

      # Define potential heights within tolerance range
      splits <- c(p_ht, p_ht - tol2, p_ht + tol2)

      # Find valid heights within internode ranges
      split_hts <- purrr::discard(splits, function(sp) {
        if (is.na(sp)) {
          TRUE
        } else {
          !(sp >= int1[1] && sp <= int1[2] || sp >= int2[1] && sp <= int2[2])
        }
      })

      # Choose split_int2: either equal-distance height or closest valid height
      suppressWarnings({
        split_int2 <- if (p_ht %in% split_hts) {
          p_ht
        } else {
          min(split_hts, na.rm = TRUE)
        }
      })

      # Fallback if no valid height is found
      if (is.infinite(split_int2)) {
        c <- ifelse(split_int1 >= bh, 1, 2)
        split_int2 <- get(split_int1_sect)[c]
      }
    }

    # Round split heights to four decimal places
    split_int1 <- round(split_int1, 2)
    split_int2 <- round(split_int2, 2)
  }

  # Get dbh if both dbh_class and nearest_int1 are available.
  if (!is.na(dbh_class) && !is.na(nearest_int1)) {
    if (dbh_class == "Split") {
      d_profile <- d_profile %>%
        dplyr::mutate(
          z_abs_int1 = abs(z - nearest_int1),
          z_abs_split_int1 = abs(z - split_int1),
          z_abs_split_int2 = abs(z - split_int2),
          dplyr::across(dplyr::where(is.numeric), round, 2)
        )
    }

    # Extract dbh from original diameter profile based on dbh_class.
    dbh <- switch(dbh_class,
      "DirectM" = d_profile %>%
        dplyr::filter(z == 1.4),
      "DirectSS" = d_profile %>%
        dplyr::slice(which.min(z_abs_int1)) %>%
        dplyr::mutate(z = dplyr::if_else(
          z > bh, z - seg_ht / 2, z + seg_ht / 2
        )),
      "DirectLS" = d_profile %>%
        dplyr::slice(which.min(z_abs_int1)) %>%
        dplyr::mutate(z = dplyr::if_else(
          z > bh, z - seg_ht / 2, z + seg_ht / 2
        )),
      "Split" = {
        # Retrieve closest diameters to split heights.
        dbh1 <- d_profile %>%
          dplyr::slice(which.min(z_abs_split_int1))
        dbh2 <- d_profile %>%
          dplyr::slice(which.min(z_abs_split_int2))

        # Combine and calculate mean diameter for Split class.
        dplyr::bind_rows(dbh1, dbh2) %>%
          dplyr::mutate(z = dplyr::if_else(
            z > bh, z - seg_ht / 2, z + seg_ht / 2
          )) %>%
          dplyr::summarise(dplyr::across(dplyr::where(is.numeric), mean))
      },
      stop("Unknown dbh_class")
    )

    # Final output structure.
    tree_dbh_var <- dbh %>%
      dplyr::mutate(dbh_class = dbh_class) %>%
      dplyr::select(dbh_class, dbh_height = z, diameter)
  } else {
    # Return NA values if dbh_class or nearest_int1 is not available.
    tree_dbh_var <- dplyr::tibble(
      dbh_class = NA,
      dbh_height = NA,
      diameter = NA
    )
  }

  return(tree_dbh_var)
}
