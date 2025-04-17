#' Validate various input types
#'
#' Validate directory path
#'
#' Checks if a path is a character string, if it exists, or if it can be
#' created.
#'
#' @param path Character. The path to validate.
#' @param create Logical. If \code{TRUE}, attempts to create the directory if
#'               given path does not exist
#'               Default is FALSE.
#'
#' @importFrom utils globalVariables
#'
#' @return NULL if valid, otherwise throws an error.
#' @export
#'
#' @name validate_path
#'
utils::globalVariables("name")

validate_path <- function(path, create = FALSE) {
  name <- if (is.symbol(substitute(path))) {
    as.character(substitute(path))
  } else {
    path
  }

  if (is.null(path) || !is.character(path)) {
    stop("Invalid input: path must be a character string.")
  }

  if (!dir.exists(path) && !file.exists(path)) {
    if (create) {
      dir.create(path, showWarnings = FALSE, recursive = TRUE)
      if (!dir.exists(path)) {
        stop(paste0(
          "Invalid input: Unable to create the specified directory
          or file path '",
          name, "'."
        ))
      }
    } else {
      stop(paste0(
        "Invalid input: The specified directory or file path '",
        name, "' does not exist."
      ))
    }
  }
  invisible(NULL)
}

#' Validate numeric input
#'
#' Checks if a value is numeric and within a specified range.
#'
#' @param x Numeric. The value to validate.
#' @param min Numeric. Optional minimum value.
#' @param max Numeric. Optional maximum value.
#' @return NULL if valid, otherwise throws an error.
#' @export
#'
validate_numeric <- function(x, min = -Inf, max = Inf) {
  name <- as.character(substitute(x))

  if (is.null(x) || !is.numeric(x)) {
    stop(paste0("Invalid input: ", name, " value must be numeric."))
  }

  if (x < min || x > max) {
    stop(paste0(
      "Invalid input: ", name, " value must be between ", min, " and ",
      max, "."
    ))
  }

  invisible(NULL)
}

#' Validate Coordinate Reference System (CRS) EPSG code
#'
#' Checks if a CRS code is numeric and within the EPSG code range.
#'
#' @param defined_crs Numeric. The EPSG code to validate.
#' @return NULL if valid, otherwise throws an error.
#' @export
#'
validate_crs <- function(defined_crs) {
  name <- as.character(substitute(defined_crs))

  if (is.null(defined_crs) || !is.numeric(defined_crs) || defined_crs < 1000 ||
    defined_crs > 99999) {
    stop(paste0("Invalid input: ", name, " must be a numeric EPSG code
                between 1000 and 99999."))
  }

  invisible(NULL)
}

#' Validate integer input
#'
#' Checks if a value is an integer and, optionally, if it is positive.
#'
#' @param x Numeric. The value to validate.
#' @param positive Logical. If TRUE, ensures the integer is positive
#'                 Default is TRUE.
#' @return NULL if valid, otherwise throws an error.
#' @export
#'
validate_integer <- function(x, positive = TRUE) {
  name <- as.character(substitute(x))

  if (is.null(x) || !is.numeric(x) || x != as.integer(x)) {
    stop(paste0("Invalid input: ", name, " value must be an integer."))
  }

  if (positive && x <= 0) {
    stop(paste0("Invalid input: ", name, " value must be a positive integer."))
  }

  invisible(NULL)
}

#' Validate character input
#'
#' Checks if a value is a non-null, non-empty character string.
#'
#' @param x Character. The value to validate.
#' @return NULL if valid, otherwise throws an error.
#' @export
#'
validate_character <- function(x) {
  name <- as.character(substitute(x))

  if (is.null(x) || !is.character(x) || x == "") {
    stop(paste("Invalid input:", name, "must be a non-null, non-empty
               character string."))
  }

  invisible(NULL)
}

#' Validate LAS object
#'
#' Checks if the provided object is a valid `LAS` object from the `lidR`
#' package, is not empty, and contains necessary attributes.
#'
#' @param las An object expected to be of class `LAS` from the `lidR` package.
#' @param required_cols A vector containing all required column names.
#'                      Default is \code{c("X", "Y", "Z")}
#' @return NULL if valid, otherwise throws an error.
#' @export
validate_las <- function(las, required_cols = c("X", "Y", "Z")) {
  name <- as.character(substitute(las))

  if (class(las)[1] != "LAS") {
    stop(paste0("Invalid input: ", name, " must be a `LAS` object."))
  }

  if (nrow(las@data) == 0) {
    stop(paste0("Invalid input: ", name, " LAS object is empty and
                contains no points."))
  }

  missing_columns <- setdiff(required_cols, names(las@data))

  if (length(missing_columns) > 0) {
    stop(paste0(
      "Invalid input: ", name,
      " LAS object is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    ))
  }

  invisible(NULL)
}

#' Validate logical input
#'
#' Checks if the provided input is a valid boolean/logical (TRUE or FALSE) value.
#'
#' @param x The input value to be validated.
#' @return NULL if valid, otherwise throws an error.
#' @export
validate_logical <- function(x) {
  name <- as.character(substitute(x))

  if (is.null(x) || !is.logical(x)) {
    stop(paste(
      "Invalid input:", name,
      "must be a logical value (TRUE or FALSE)."
    ))
  }

  invisible(NULL)
}

#' Validate vector input
#'
#' Checks if the provided input is a vector of the specified length.
#'
#' @param x The input value to be validated, expected to be a vector.
#' @param length The required length of the vector.
#' @param name (optional) The name of the variable, used for error messages.
#' @return NULL if valid, otherwise throws an error.
#' @export
validate_vector <- function(x, length, name) {
  if (is.null(x) || !is.vector(x)) {
    stop(paste("Invalid input:", name, "must be a vector."))
  }

  if (length(x) != length) {
    stop(paste(
      "Invalid input:", name, "must be a vector of length",
      length, "."
    ))
  }

  invisible(NULL)
}

#' Validate DataFrame input with specific columns
#'
#' Checks if a given DataFrame contains the required columns.
#'
#' @param df A data frame that is expected to contain certain columns.
#' @param required_cols A character vector of column names that the
#'                      DataFrame must have.
#' @return NULL if valid, otherwise throws an error.
#' @export
validate_dataframe <- function(df, required_cols) {
  if (!("data.frame" %in% class(df))) {
    stop(paste0("Invalid input: ", name, " must be a `data.frame` object."))
  }

  # Check if the required columns exist in the DataFrame
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop(paste0(
      "Invalid input: The DataFrame is missing the
                following required columns: ",
      paste(missing_cols, collapse = ", ")
    ))
  }

  invisible(NULL)
}
