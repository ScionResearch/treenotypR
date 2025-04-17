# Purpose:
## To create a canopy height model (CHM).
## Run this script to create a CHM from a normalised point cloud.

# Inputs:
## Normalised LiDAR or SfM point cloud (laz)

# Outputs:
## A CHM raster (tif).

# Author:
## Sadeepa Jayathunga (sadeepa.jayathunga@scionresearch.com)©

# Last updated:
Sys.time()
"2024-11-11 20:04:52 NZDT"

################################################################################
# Clear global environment and set up options.
rm(list = ls())
options(warn = 0)
par(mar = c(1, 1, 1, 0))

################################################################################
# Requirements.
## Requirement1: Install (only if required) and load the packages.
pkg_check <- function(pkgs) {
  missing_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, repos = "https://cran.rstudio.com/")
  }
  sapply(pkgs, library, character.only = TRUE)
}

pkgs <- c("rstudioapi", "future", "lidR", "terra", "raster", "dplyr")

pkg_check(pkgs)

################################################################################
## Requirement2: Set the working directory.
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

## Requirement3: Set up multi-processing.
future::plan(future::multisession, workers = 20)
lidR::set_lidr_threads(20)

################################################################################
## Requirement4: Set input path and parameters.
## Add site name. This will be used when assigning output names.
site_name <- "Test"

## Set input path of raw data.Can be a single LAS file or a folder with tiles.
## Must have a file extension if the input is a single file.
in_path_las <- file.path(getwd(), "Test_data", "Norm_ULS_P1.laz")

## Set CHM resolution in meters.
chm_res <- 0.25

## Set coordinate reference system (CRS).
defined_crs <- 2193

## Set moving window size.
## Options: 3x3 pixels (default), 5x5 pixels or 7x7 pixels.
ws <- matrix(1, 3, 3)

################################################################################
# Functions.
## Function1: Fill no data pixels.
fill.na <- function(x, i = 5) {
  if (is.na(x)[i]) {
    return(mean(x, na.rm = TRUE))
  } else {
    return(x[i])
  }
}

################################################################################
## Create output folder and set output path.
## Create main output directory, set out path and create sub directories to
## store intermediate results.
out_path <- file.path(getwd(), "CHM")

if (!dir.exists(out_path)) {
  dir.create(out_path)
}

################################################################################
# Main script.

## 1.Load input las data.

## Load input las file as a catalog and plot it with chunk pattern preview.
## Purpose: catalog format is used to avoid keeping tiles in the computer memory
## and allowing user to read files only when necessary.
## to identify any inconsistencies.
ctg <- lidR::readLAScatalog(in_path_las)
lidR::las_check(ctg)
lidR::plot(ctg, chunk = TRUE)

## Set coordinate reference system (CRS).
lidR::st_crs(ctg) <- defined_crs

################################################################################
## 2. Create a canopy height model (CHM) raster.

## Save CHM tiles to a directory for future reference.
## Purpose: To save output to a disk to avoid overloading the computer memory
## and allowing the user to read files only when necessary.
dir.create(file.path(out_path, "CHM_Tiles"))

lidR::opt_output_files(ctg) <- paste0(
  out_path,
  "/Tiles/CHM_{XLEFT}_{YBOTTOM}"
)

## Filter out noise points from the denoised point cloud.
lidR::opt_filter(ctg) <- "-drop_Z_below 0 "

## Create a merged CHM raster.
## Default CHM resolution is 0.25m but must be altered to a higher value if
## the density of the point cloud is significantly lower (e.g., <10 points/sqm).
chm_algo <- lidR::pitfree(
  thresholds = c(0, 10, 20),
  max_edge = c(0, 1.5),
  subcircle = 0.1,
  highest = TRUE
)

## Perform rasterisation.
chm_tiles <- lidR::rasterize_canopy(ctg,
  res = chm_res,
  algorithm = chm_algo
)

## Post-process CHM raster.
chm <- terra::focal(chm_tiles, ws, fun = fill.na) %>%
  raster::raster()

## Set coordinate system.
raster::crs(chm) <- defined_crs

## Save merged CHMs to a directory.
terra::writeRaster(chm,
  filename = paste0(
    out_path, "/Raw_CHM_", site_name, ".tif"
  ),
  overwrite = TRUE
)

## Remove unnecessary objects
rm(chm_tiles)

################################################################################
# End of Script.
