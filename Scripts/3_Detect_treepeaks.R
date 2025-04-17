# Purpose:
## To detect treepeaks.
## Run this script to detect treepeaks on both the CHM and the point cloud.

# Inputs:
## Normalised LiDAR or SfM point cloud (laz); CHM (tif)

# Outputs:
## Treepeak vectors (shp).

# Author:
## Sadeepa Jayathunga (sadeepa.jayathunga@scionresearch.com)©

# Last updated:
Sys.time()
"2024-11-11 20:34:49 NZDT"

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

pkgs <- c("rstudioapi", "future", "lidR", "terra", "raster", "sf", "magrittr")

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
in_path_las <- file.path(getwd(), "Test_data", "Norm_MLS_P1.laz")

## Set input path of CHM.
in_path_chm <- file.path(getwd(), "Test_data", "Raw_CHM_Test.tif")

## Set minimum distance between two points in meters.
min_dist <- 3

## Set minimum height of detected treepeaks.
min_ht <- 2

## Set thining resolution.
thin_res <- 0.25

## Set coordinate reference system (CRS).
defined_crs <- 2193

################################################################################
## Create output folder and set output path.
## Create main output directory, set out path and create sub directories to
## store intermediate results.
out_path <- file.path(getwd(), "Treepeaks")

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
## 2.Thin point cloud to a manageable size.

## Save thinned tiles to a directory.
## Purpose: To save output to a disk to avoid overloading the computer memory
## and allowing the user to read files only when necessary.
dir.create(file.path(out_path, "Thinned_Tiles"))

lidR::opt_output_files(ctg) <- paste0(
  out_path,
  "/Thinned_Tiles/Thinned_{XLEFT}_{YBOTTOM}"
)

## Perform point cloud thinning.
thinned_ctg <- lidR::decimate_points(ctg, lidR::highest(res = thin_res))

################################################################################
## 3. Detect treepeaks on the point cloud.

## Save tree peaks detected on the point cloud to a separate directory.
## Purpose: To save output to a disk to avoid overloading the computer memory
## and allowing the user to read files only when necessary.
dir.create(file.path(out_path, "Treepeaks_LAS_Tiles"))

lidR::opt_output_files(thinned_ctg) <- paste0(
  out_path,
  "/Treepeaks_LAS_Tiles/Treepeaks_LAS_{XLEFT}_{YBOTTOM}"
)

## Detect tree peaks on the point cloud.
treepeaks_las <- lidR::locate_trees(
  thinned_ctg,
  lidR::lmf(
    ws = min_dist,
    hmin = min_ht,
    shape = c("circular", "square")
  )
)

## Combine all detected tree peaks to one object.
treepeaks_las_merged <- lapply(treepeaks_las, st_read) %>%
  dplyr::bind_rows()

treepeaks_las_merged <- sf::st_transform(treepeaks_las_merged, 2193)

## Save detected tree peaks.
sf::write_sf(treepeaks_las_merged, paste0(
  out_path,
  "/Treepeaks_LAS_", site_name, ".shp"
))

################################################################################
## 4. Detect treepeaks on the CHM.

## Read and prep CHM.
chm <- raster::raster(in_path_chm)
raster::crs(chm) <- defined_crs

## Detect tree peaks on the CHM.
treepeaks_chm <- lidR::locate_trees(
  chm,
  lidR::lmf(
    ws = 3,
    hmin = 2,
    shape = c("circular", "square")
  )
) %>%
  sf::st_zm()

## Set coordinate system.
treepeaks_chm <- sf::st_transform(treepeaks_chm, 2193)

## Save detected tree peaks.
sf::write_sf(treepeaks_chm, paste0(
  out_path,
  "/Treepeaks_CHM_", site_name, ".shp"
))

################################################################################
# End of Script.
