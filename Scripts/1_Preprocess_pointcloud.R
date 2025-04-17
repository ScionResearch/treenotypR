# Purpose:
## To pre-process point cloud data.
## Run this script first to de-noise the point clouds, identify ground and
## ground normalise tree heights as per the following steps.

# Inputs:
## Raw LiDAR or SfM point cloud (laz); if SfM point cloud a DTM raster (tif)

# Steps:
## 1. Pre-process point clouds
## 1.1 Tile,
## 1.2. De-noise,
## 1.3. Ground classification,
## 1.4. Create a digital terrain model (DTM),
## 1.5. Height normalisation,
## 1.6. Create a canopy height model (CHM),

# Outputs:
## Final output will be a ground normalised las dataset.
## Products of all intermediate steps will be saved to the sub directories
## created within the Processed_pointclouds directory.

# Author:
## Sadeepa Jayathunga (sadeepa.jayathunga@scionresearch.com)©

# Last updated:
Sys.time()
"2024-08-28 11:55:38 NZST"

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

pkgs <- c("rstudioapi", "future", "lidR", "terra", "raster", "magrittr")

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

## Specify point cloud type: LiDAR; SfM.
## LiDAR data can be point cloud captured using airborne laser scanning (ALS),
## unmanned aerial vehicle laser scanning (ULS), mobile laser scanning (MLS),
## terrestrial laser scanning (TLS).
## MLS and TLS are usually ultra-high density and needs to be decimated to a
## manageable size to optimise computer resources.
## ALS, ULS, SfM, MLS and TLS works differently when capturing vegetation and
## have  their pros and cons. We highly recommend users get a proper
## understanding of these techniques and data sets.
data_type <- "LiDAR"

## Specify processing steps to be performed.
tile <- FALSE
denoise <- TRUE
ground_classify <- TRUE
create_dtm <- TRUE
height_normalise <- TRUE

## Set input path of raw data.Can be a single LAS file or a folder with tiles.
## Must have a file extension if the input is a single file.
in_path_las <- file.path(getwd(), "Test_data", "NonNorm_MLS_P1.laz")

## Set input path of DTM (only needed for SfM point clouds).
## SfM point clouds often fails to reconstruct the terrain in forest areas with
## closed canopy cover.
in_path_dtm <- file.path(getwd(), "Test_data", "Test_DTM.tif")

## Set coordinate reference system (CRS).
defined_crs <- 2193

################################################################################
# Functions.
## Function1: Create directories to store results.
dirs <- function(path, name) {
  if (!dir.exists(file.path(path, name))) {
    dir.create(file.path(path, name), showWarnings = FALSE) %>%
      invisible()
  } else {
    warning("Directory exists. Files may be overwritten.")
  }
}

## Function2: Plot ground classified point cloud for inspection.
plot_classified_pc <- function(filename, color_by, split = NA) {
  las_all <- lidR::readLAS(filename)
  lidR::plot(las_all, color = color_by)

  if (!is.na(split)) {
    las_split <- filter_poi(las_all, eval(parse(text = split)))
    lidR::plot(las_split, color = color_by)
    rm(las_split)
  }

  rm(las_all)
}

## Function3: Fill no data pixels.
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
dirs(getwd(), "Processed_pointclouds")

out_path <- file.path(getwd(), "Processed_pointclouds")

directory_names <- c("1_Tiles", "2_Denoised", "3_GroundC", "4_DTM", "5_NormHt")

lapply(directory_names, function(dir_name) dirs(out_path, dir_name))

################################################################################
# Main script.
## 1. Pre-processing point clouds.

## 1.1 Load input las file as a catalog and plot it with chunk pattern preview.
## Purpose: catalog format is used to avoid keeping tiles in the computer memory
## and allowing user to read files only when necessary.
## to identify any inconsistencies.
ctg <- lidR::readLAScatalog(in_path_las)
lidR::las_check(ctg)
lidR::plot(ctg, chunk = TRUE)

## Set coordinate reference system (CRS).
lidR::st_crs(ctg) <- defined_crs

################################################################################
## 1.2 Split the las file into a set of tiles/chunks.

if (tile) {
  ## Purpose: Tiling is done to efficiently allocate computer cores/memory for
  ## for processing. This also improves the processing speed.
  ## If needed, user can alter the tile size and tile buffer.
  lidR::opt_chunk_size(ctg) <- 50
  lidR::opt_chunk_buffer(ctg) <- 0

  ## Save new tiles to a directory.
  ## Purpose: To save output to a disk to avoid overloading the computer memory
  ## and allowing the user to read files only when necessary.
  ## Also set output format to .laz (compressed .las) to save disk space.
  lidR::opt_output_files(ctg) <- paste0(
    out_path,
    "/1_Tiles/Tile_{XLEFT}_{YBOTTOM}"
  )

  lidR::opt_laz_compression(ctg) <- TRUE

  ## Perform tiling with a preview of the new chunk pattern.
  new_ctg <- catalog_retile(ctg)

  ## Set catalog path for the next step.
  ctg <- new_ctg
  rm(new_ctg)
}

################################################################################
## 1.3 De-noise point clouds.

if (denoise) {
  ## Save de-noised tiles to a directory.
  ## Purpose: To save output to a disk to avoid overloading the computer memory
  ## and allowing the user to read files only when necessary.
  lidR::opt_output_files(ctg) <- paste0(
    out_path,
    "/2_Denoised/Denoised_{XLEFT}_{YBOTTOM}"
  )
  
  lidR::opt_laz_compression(ctg) <- TRUE

  ## De-noising algorithm: Statistical outlier removal.
  ## User can specify the number of neighbours (k) and maximum distance(m).
  denoise_algo <- lidR::sor(k = 15, m = 7)

  ## Perform de-noising.
  denoised_ctg <- lidR::classify_noise(ctg,
    algorithm = denoise_algo
  )

  ## Set catalog path for the next step.
  ctg <- denoised_ctg
  rm(denoised_ctg)
}

################################################################################
## 1.4 Classify ground points.

if (ground_classify & data_type == "LiDAR") {
  ## Save ground classified tiles to a directory.
  ## Purpose: To save output to a disk to avoid overloading the computer memory
  ## and allowing the user to read files only when necessary.
  lidR::opt_output_files(ctg) <- paste0(
    out_path,
    "/3_GroundC/GroundC_{XLEFT}_{YBOTTOM}"
  )
  
  lidR::opt_laz_compression(ctg) <- TRUE

  ## Filter out noise points from the denoised point cloud.
  lidR::opt_filter(ctg) <- "-drop_classification 18"

  ## Classification algorithm: Two options available.
  ## Purpose: Choosing an appropriate algorithm and identifying optimal
  ## parameters are not trivial tasks. Both these tasks often require several
  ## trial runs.
  ## We have included two algorithms below but we cannot say which one is better
  ## for a certain forest stand without inspecting the point clouds first. Also
  ## we cannot say which parameters would best suit the terrain of that
  ## particular forest stand. It's likely that users may need to dynamically
  ## adjust default parameters based on the characteristics of the local
  ## terrain.

  ## Option 1: Multiscale Curvature Classification algorithm used two
  ## parameters: a scalar parameter (s) and a curvature threshold (t).
  ground_algo1 <- lidR::mcc(s = 1.5, t = 0.3)

  ## Option 2: Cloth Simulation Filter algorithm uses six parameters.
  ground_algo2 <- lidR::csf(
    sloop_smooth = FALSE,
    class_threshold = 0.5,
    cloth_resolution = 0.5,
    rigidness = 1L,
    iterations = 500L,
    time_step = 0.65
  )

  ## Perform ground point classification using the chosen algorithm.
  groundC_ctg <- lidR::classify_ground(ctg,
    algorithm = ground_algo2,
    last_returns = TRUE
  )

  ## Once the ground points are classified, the user can plot a few random tiles
  ## and do visual assessment of the classification accuracy on R 3D plot window.
  ## In addition, user can load tiles onto CloudCompare software package and
  ## carry out the inspection. Cloudcompare has a very user-friendly interface
  ## and thus allows the user to conduct more interactive visualisations and
  ## assessments.

  ## Run following line of code as many times as necessary to choose and plot
  ## random tiles from the catalog. Two plots will be created for each chosen
  ## tile: one with all the classes and one only with ground class.
  plot_classified_pc(sample(groundC_ctg$filename, 1),
    color_by = "Classification",
    split = "Classification == 2L"
  )

  ## If re-running ground classification (with a different algorithm or refined
  ## parameters) is required, following line of code can be used to empty the
  ## 2_GroundC sub directory by removing previous results.
  # unlink(paste0(out_path, "/3_GroundC"), recursive = TRUE)

  ## Set catalog path for the next step.
  ctg <- groundC_ctg
  rm(groundC_ctg)
}

################################################################################
## 1.5 Create a digital terrain model (DTM) from the classified ground points.
## Default DTM resolution is 1m but must be altered to a higher value if the
## density of the point cloud is significantly lower (e.g., <10 points/sqm).

if (create_dtm) {
  ## Triangulation and extrapolation algorithm: k-nearest neighbour (KNN)
  ## approach with an inverse-distance weighting (IDW) is used as the
  ## extrapolation method.
  ## Parameters in the triangulation and extrapolation algorithm includes the
  ## number of k-nearest neighbours (k), power for inverse-distance weighting
  ## (p) and maximum search radius (rmax).
  dtm_algo <- lidR::tin(extrapolate = lidR::knnidw(k = 10, p = 2, rmax = 10))

  ## Save DTM tiles to a sub-directory for future reference.
  ## Purpose: To save output to a disk to avoid overloading the computer memory
  ## and allowing the user to read files only when necessary.
  dirs(file.path(getwd(), "Preprocessed_pointclouds", "4_DTM"), "Tiles")

  lidR::opt_output_files(ctg) <- paste0(
    out_path,
    "/4_DTM/Tiles/DTM_{XLEFT}_{YBOTTOM}"
  )

  ## Create a merged DTM raster.
  dtm_tiles <- lidR::rasterize_terrain(ctg,
    res = 1,
    algorithm = dtm_algo,
    use_class = 2L,
    pkg = "terra"
  )

  ## Post-process the DTM raster.
  ws <- matrix(1, 3, 3)
  dtm <- terra::focal(terra::focal(dtm_tiles, ws, fun = fill.na),
    w = ws,
    fun = mean,
    na.rm = TRUE
  ) %>%
    raster::raster()

  ## When using UAV-SfM point clouds, one might not get sufficient ground points
  ## create an accurate DTM. This is particularly true in complex terrains and
  ## in trials with closed canopies. In such, cases an accurate DTM created
  ## using LiDAR data is necessary. This line of code can be used to import an
  ## external DTM into the process.

  ## Set coordinate reference system (CRS).
  raster::crs(dtm) <- defined_crs

  ## Save merged DTM raster to a directory.
  terra::writeRaster(dtm,
    filename = paste0(out_path, "/4_DTM/DTM_", site_name, ".tif"),
    overwrite = TRUE
  )

  ## Remove unnecessary objects.
  rm(dtm_tiles)
}

################################################################################
## 1.6 Ground normalise the point cloud to estimate above ground heights.

if (height_normalise) {
  ## Save normalised tiles to a directory.
  ## Purpose: To save output to a disk to avoid overloading the computer memory
  ## and allowing the user to read files only when necessary.
  lidR::opt_output_files(ctg) <- paste0(
    out_path,
    "/5_NormHt/NormHt_{XLEFT}_{YBOTTOM}"
  )
  
  lidR::opt_laz_compression(ctg) <- TRUE

  ## Use appropriate dtm depending on the data_type.
  if (data_type == "SfM") {
    dtm <- raster::raster(in_path_dtm)
    raster::crs(dtm) <- defined_crs
  }

  ## Perform height normalisation.
  normHt_ctg <- lidR::normalize_height(ctg, dtm)

  ## Run following line of code as many times as necessary to choose and plot
  ## random tiles from the catalog.
  plot_classified_pc(sample(normHt_ctg$filename, 1),
    color_by = "Z",
    split = NA
  )

  ## Set catalog path for the next step.
  ctg <- normHt_ctg
  rm(normHt_ctg)
}

################################################################################
# End of Script.
