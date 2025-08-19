#' Stem detection from LAS individual tree point cloud
#'
#' This function detects and classifies stem points in an individual tree
#' LAS file, allowing restriction of the height range based on user input.
#' It uses the Hough transformation method to identify and classify the stem
#' cylinder in the point cloud.
#'
#' @param las_path Character. The file path to the LAS file.
#' @param stem_only Logical. If \code{TRUE}, returns only stem points if
#'                  \code{FALSE}, returns all points classified as Stem and
#'                  non-Stem.
#'                  Defaults to \code{TRUE}.
#' @param min_ht Numeric. The minimum height threshold in meters for filtering
#'               stem points. Points below this height will be excluded.
#' @param max_ht Numeric or character. The maximum height limit in meters for
#'               filtering stem points, either as "zmax" to use the file's
#'               maximum Z value or a specific numeric height.
#'               Default is \code{"zmax"}.
#' @param all_attr Logical. If \code{TRUE}, reads all attributes from the LAS
#'                 file and keep them through the processing if \code{FALSE},
#'                 reads and retains only XYZ coordinates.
#'                 Defaults to \code{TRUE}.
#' @param ht_step Numeric. Height step interval (in meters) used in the Hough
#'                transform for detecting stem segments. A smaller value
#'                increases accuracy but may increase computation time.
#'                Default is 0.1.
#' @param max_diameter Numeric. Maximum allowable stem diameter in meters
#'                     for classification. Points representing stems larger than
#'                     this diameter will be excluded.
#'                     Default is 0.6.
#' @param ht_base Numeric vector of length 2. Specifies the height range
#'                (in meters) considered as the base of the tree for stem
#'                detection.
#'                Default is \code{c(1, 2)}.
#' @param pixel_res Numeric. Pixel resolution (in meters) for grid-based stem
#'                  classification.
#'                  Default is 0.02 meters.
#' @param min_den Numeric. Minimum density of points (per square meter) required
#'                for an area to be considered part of the stem.
#'                Default is 0.2.
#' @param min_numvotes Integer. Minimum number of votes in the Hough transform
#'                     required to confirm stem points. Lower values may
#'                     increase false positives.
#'                     Default is 3.
#' @param defined_crs Numeric. Coordinate reference system (CRS) for the output
#'                    LAS object.
#'                    Default is EPSG:2193.
#'
#' @return A LAS object containing detected stem and non-stem (if specified)
#'         points that meet the specified height and classification criteria.
#'
#' @details
#' This function performs stem detection on a LAS file by:
#' \itemize{
#'   \item Reads the LAS file and assigns the specified CRS.
#'   \item Checks for missing values in the X, Y, and Z columns.
#'   \item Filters points based on the height range defined by \code{min_ht}
#'         and \code{max_ht}.
#'   \item Classify stem points using the Hough transformation, restricting
#'         results based on parameters like \code{max_diameter}, \code{ht_step},
#'         and \code{min_numvotes}.
#'   \item Returns a LAS object containing only the detected stem points
#'         meeting these criteria.
#' }
#'
#' @importFrom lidR readLAS las_check st_crs filter_poi filter_duplicates
#' @importFrom TreeLS stemPoints stm.hough
#' @importFrom dplyr filter mutate_if distinct
#' @importFrom utils globalVariables
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' stem_points <- detect_stem(
#'   las_path = "path/to/your/file.laz",
#'   stem_only = TRUE
#' )
#' lidR::plot(stem_points)
#' }
#'
#' @export
#'
#' @name detect_stem
#'
utils::globalVariables(c("%>%", "Z", "X", "Y", "Stem"))

detect_stem <- function(
    las_path,
    stem_only = TRUE,
    min_ht = 0.3,
    max_ht = "zmax",
    all_attr = TRUE,
    ht_step = 0.1,
    max_diameter = 0.6,
    ht_base = c(1, 2),
    pixel_res = 0.02,
    min_den = 0.2,
    min_num_votes = 3,
    defined_crs = 2193) {
  # Validate inputs.
  validate_path(las_path)
  validate_numeric(min_ht, min = 0, max = 100)
  validate_logical(stem_only)
  validate_logical(all_attr)
  validate_numeric(ht_step, min = 0.0001, max = 5)
  validate_numeric(max_diameter, min = 0.1, max = 5)
  validate_vector(ht_base, 2, "ht_base")
  validate_numeric(pixel_res, min = 0.0001, max = 5)
  validate_numeric(min_den, min = 0.0001, max = 10)
  validate_integer(min_num_votes)
  validate_crs(defined_crs)

  # Load LAS file with error handling.
  las <- tryCatch(
    {
      lidR::readLAS(las_path, select = if (all_attr) "*" else "xyz")
    },
    error = function(e) {
      stop("Error reading LAS file: ", e$message)
    }
  )

  # Check and validate max_ht with the point cloud.
  max_ht <- ifelse(
    max_ht == "zmax" || as.numeric(max_ht) > max(las$Z),
    round(max(las$Z), 4),
    as.numeric(max_ht)
  )

  validate_numeric(max_ht, min = 0.01, max = 100.01)

  # Ensure max_ht is greater than min_ht.
  if (min_ht >= max_ht) {
    stop("Invalid height range: min_ht must be less than max_ht.")
  }

  # Filter duplicate LAS points and points within the specified Z height range.
  las <- tryCatch(
    {
      las <- lidR::filter_duplicates(las)
      las <- lidR::filter_poi(las, Z >= min_ht & Z <= max_ht)
      las
    },
    error = function(e) {
      message("Error filtering LAS file by Z values: ", e$message)
      return(NULL)
    }
  )

  # Check validity of the LAS file.
  las_check_result <- lidR::las_check(las, print = FALSE)
  if (length(las_check_result$warnings) > 0 || length(
    las_check_result$errors
  ) > 0) {
    warning(
      "The input LAS file may have file las specification issues: ",
      paste(c(las_check_result$warnings, las_check_result$errors),
        collapse = "; "
      )
    )
  }

  # Remove rows with NA in X, Y, or Z only, and round numeric columns.
  las@data <- las@data %>%
    dplyr::filter(!is.na(X) & !is.na(Y) & !is.na(Z)) %>%
    dplyr::mutate_if(is.numeric, round, 4) %>%
    dplyr::distinct()

  # Assign CRS to the LAS file
  tryCatch(
    {
      lidR::projection(las) <- defined_crs
    },
    error = function(e) {
      stop("Error assigning CRS to LAS file: ", e$message)
    }
  )

  # Classify stem points using TreeLS.
  stem_points_class <- tryCatch(
    {
      TreeLS::stemPoints(
        las, TreeLS::stm.hough(
          h_step = ht_step,
          max_d = max_diameter,
          h_base = ht_base,
          pixel_size = pixel_res,
          min_density = min_den,
          min_votes = min_num_votes
        )
      )
    },
    error = function(e) {
      stop("Error in TreeLS stem detection: ", e$message)
    }
  )

  if (stem_only) {
    # Filter for stem points above specified min_ht.
    stem_points_class <- tryCatch(
      {
        lidR::filter_poi(
          stem_points_class,
          Stem == TRUE & Z >= min_ht & Z <= max_ht
        )
      },
      error = function(e) {
        stop("Error filtering stem points: ", e$message)
      }
    )
  }

  # Adjust column types if all attributes were read.
  if (all_attr) {
    all_cols <- c(
      "Intensity", "ReturnNumber", "NumberOfReturns",
      "ScanDirectionFlag", "EdgeOfFlightline", "Classification",
      "ScanAngleRank", "UserData", "PointSourceID"
    )
    available_cols <- intersect(names(stem_points_class), all_cols)

    stem_points_class@data <- stem_points_class@data %>%
      dplyr::mutate(dplyr::across(dplyr::any_of(available_cols), as.integer))
  }

  # Convert logical 'Stem' attribute to numeric
  stem_points_class@data$Stem <- as.numeric(stem_points_class@data$Stem)

  # Add the 'Stem' attribute to the LAS header
  stem_points_class <- add_lasattribute(stem_points_class, name = "Stem", desc = "Stem_points")

  # Return the filtered stem points.
  return(stem_points_class)
}
