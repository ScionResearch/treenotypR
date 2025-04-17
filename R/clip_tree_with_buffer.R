#' Clip individual tree segments to circular buffer from stem
#'
#' This function clips individual tree segments from a LAS object using a
#' circular buffer centering at stem.
#'
#' @param stem_class_las A LAS object of an individual tree with stem points
#'                       detected and classified.
#' @param buffer Numeric. The buffer distance from stem (in meters) added to the
#'               radius of the clipped segments.
#'               Default is 0.5.
#' @param min_ht Numeric. Minimum stem height for segmentation in meters.
#'               Points below this height will be excluded.
#'               Default is 0.2.
#' @param max_ht Numeric or character. Maximum height for segmentation, either
#'               as "zmax" (to use the maximum height in the LAS file) or a
#'               specific numeric value.
#'               Default is \code{"zmax"}.
#' @param defined_crs Numeric. Coordinate reference system (CRS) for output
#'                    layers.
#'                    Default is EPSG:2193.
#'
#' @return A LAS object containing the clipped tree segment.
#'
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Filters the LAS data to remove duplicated and NAs.
#'   \item Identifies stem points using Hough transformation.
#'   \item Loops through height segments to fit circular shapes using RANSAC and
#'         clips points within a buffer distance.
#'   \item Returns a LAS object with all clipped segments and removes duplicate
#'         points.
#' }
#' @importFrom dplyr bind_rows
#' @importFrom tibble tibble
#' @importFrom lidR filter_poi clip_circle filter_duplicates st_crs LAS
#' @importFrom utils tail
#' @importFrom stats na.omit
#' @importFrom utils globalVariables
#' @importFrom magrittr %>%
#'
#' @examples
#' \dontrun{
#' stem_class_las <- detect_stem("path/to/lasfile.las", stem_only = FALSE)
#' clipped_las <- clip_with_buffer(stem_class_las, buffer = 0.3)
#' lidR::plot(clipped_las)
#' }
#'
#' @export
#'
#' @name clip_with_buffer
#'
utils::globalVariables(c("Z", "Stem"))

clip_with_buffer <- function(stem_class_las,
                             buffer = 0.5,
                             min_ht = 0.3,
                             max_ht = "zmax",
                             seg_ht = 0.2,
                             n_points = 25,
                             quantiles = c(0.1, 0.9),
                             defined_crs = 2193) {
  # Validate inputs.
  validate_las(stem_class_las, required_cols = c("X", "Y", "Z", "Stem"))
  validate_numeric(buffer, min = 0.1, max = 50)
  validate_crs(defined_crs)
  validate_numeric(min_ht, min = 0.01, max = 100)

  # Define z range for segmenting
  z_min <- min_ht
  zmax <- ifelse(max_ht != "zmax", as.numeric(max_ht),
    round(max(stem_class_las$Z), 4)
  )

  validate_numeric(zmax, min = 0.01, max = 100)

  # Validate height range.
  if (min_ht >= zmax) {
    stop("Invalid Input: height range: min_ht must be less than max_ht.")
  }

  # Attempt to filter the LAS object by Z values, with error handling.
  stem_class_las <- tryCatch(
    {
      lidR::filter_poi(stem_class_las, Z <= zmax & Z >= z_min)
    },
    error = function(e) {
      message("Error filtering stem_class_las by Z values: ", e$message)
      return(NULL)
    }
  )

  # Create an empty object to store clipped segments.
  clipped_las_df <- tibble::tibble()
  diameters <- c()
  xs <- c()
  ys <- c()

  # Loop through height segments
  while (z_min < zmax) {
    z_max <- z_min + seg_ht

    # Attempt to filter by Z values first
    segment <- tryCatch(
      {
        lidR::filter_poi(stem_class_las, Z > z_min & Z < z_max)
      },
      error = function(e) {
        message("Error filtering by Z values: ", e$message)
        return(NULL) # Return NULL on error
      }
    )

    # Only proceed if the Z-filtering was successful
    if (!is.null(segment)) {
      # Now filter for points where Stem == TRUE
      segment_stem <- tryCatch(
        {
          lidR::filter_poi(segment, Stem == TRUE)
        },
        error = function(e) {
          message("Error isolating stem segment: ", e$message)
          return(NULL) # Return NULL on error
        }
      )
    } else {
      next
    }

    # Calculate diameter profile for the segment.
    if (nrow(segment_stem) > 1) {
      profile <- segment_diameter(segment_stem,
        n_points = n_points,
        quantiles = quantiles
      )
    } else {
      profile <- tibble::tibble(x = NA, y = NA, diameter = NA)
    }

    # Convert profile values to numeric once
    profile_x <- as.numeric(profile$x)
    profile_y <- as.numeric(profile$y)
    profile_diameter <- as.numeric(profile$diameter / 200 + buffer)

    # Handle missing values by combining with existing lists
    xs <- stats::na.omit(c(xs, profile_x))
    ys <- stats::na.omit(c(ys, profile_y))
    diameters <- stats::na.omit(c(diameters, profile_diameter))

    # Set values from the last elements of the lists if profile value is NA
    x <- ifelse(is.na(profile_x), utils::tail(xs, 1), profile_x)
    y <- ifelse(is.na(profile_y), utils::tail(ys, 1), profile_y)
    r <- ifelse(is.na(profile_diameter), utils::tail(diameters, 1),
      profile_diameter
    )

    # Attempt clipping the section, handling errors.
    clipped_segment <- tryCatch(
      {
        lidR::clip_circle(segment, x, y, r)
      },
      error = function(e) {
        message("Error in clipping segment with clip_circle: ", e$message)
        return(NULL)
      }
    )

    # Append clipped segment data to result LAS
    clipped_las_df <- dplyr::bind_rows(clipped_las_df, clipped_segment@data)

    # Increment the height range for next segment.
    z_min <- z_max
  }

  # Convert dataframe to a LAS object.
  clipped_las_df <- clipped_las_df %>%
    dplyr::mutate(dplyr::across(
      dplyr::any_of(
        c(
          "Intensity", "ReturnNumber", "NumberOfReturns", "ScanDirectionFlag",
          "EdgeOfFlightline", "Classification", "ScanAngleRank", "UserData",
          "PointSourceID"
        )
      ),
      as.integer
    ))

  clipped_las <- suppressMessages(lidR::LAS(clipped_las_df))

  # Attempt to filter duplicates
  clipped_las <- tryCatch(
    {
      clipped_las <- lidR::filter_duplicates(clipped_las)

      lidR::projection(clipped_las) <- defined_crs

      clipped_las
    },
    error = function(e) {
      message("Error in filtering duplicates or setting CRS: ", e$message)
      return(NULL)
    }
  )

  return(clipped_las)
}
