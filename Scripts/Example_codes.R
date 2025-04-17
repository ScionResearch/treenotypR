# Example code snippets for testing the package.

# Last updated:
Sys.time()
"2024-11-19 19:47:13 NZDT"

# If you encounter any errors or have questions pls contact
# We would like to hear from you - any feedback would be appreciated.
# Sadeepa Jayathunga(sadeepa.jayathunga@scionresearch.com)

################################################################################
# Clear environment.
rm(list = ls())

# Set library path.
# .libPaths("C:\\Users\\JayathuS\\AppData\\Local\\R\\win-library\\4.3")
# .libPaths()

################################################################################
# Set directory.
# dir <- "D:\\Projects\\ITD_Char\\Test_data"
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

# Install and load the library.
library_path <- file.path(getwd(), "treenotypR_0.1.0.tar.gz")
# detach("package:treenotypR", unload = TRUE)
devtools::install_local(library_path)
library("treenotypR")

# Set input data paths
path_MLS_ITD <- file.path(getwd(), "Test_data", "P1_1_MLS.laz")
path_MLS_plot <- file.path(getwd(), "Test_data", "Norm_MLS_P1.laz")

path_ULS_ITD <- file.path(getwd(), "Test_data", "P1_1_ULS.laz")
path_ULS_plot <- file.path(getwd(), "Test_data", "Norm_ULS_P1.laz")

path_stem_map <- file.path(getwd(), "Test_data", "Plot_stem_map.shp")

path_MLS_chm <- file.path(getwd(), "Test_data", "Raw_CHM_MLS_P1.tif")
path_ULS_chm <- file.path(getwd(), "Test_data", "Raw_CHM_ULS_P1.tif")

# Set output path.
out_path <- file.path(getwd(), "Test_Outputs")
dir.create(out_path)

################################################################################
# Process MLS point cloud.
## Decimate MLS
decimated_MLS <- decimate_LAS(path_MLS_plot, out_path)
lidR::plot(lidR::readLAS(decimated_MLS))

## Delineate crowns.
crowns_MLS <- delineate_crowns(path_MLS_chm, path_stem_map)
out_path_crowns_MLS <- file.path(out_path, "crowns_MLS.shp")
sf::write_sf(crowns_MLS, out_path_crowns_MLS)

ggplot2::ggplot(data = crowns_MLS) +
  ggplot2::geom_sf() +
  ggplot2::geom_sf_text(ggplot2::aes(label = treeID),
                        size = 3, color = "blue") +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "MLS crowns with tree IDs")

## Create density rasters.
raster_density_MLS <- create_density_raster(path_MLS_plot)
raster::plot(raster_density_MLS$raster, main = "MLS point density raster")
point_density_MLS <-  raster_density_MLS$mean_density +
  raster_density_MLS$sd_density
print(point_density_MLS)

## Homogenise the point cloud.
homogenised_MLS <- homogenise_LAS(path_MLS_plot, out_path, point_density_MLS)
lidR::plot(lidR::readLAS(homogenised_MLS))

## Clip individual tree point clouds.
ITs_MLS <- clip_IT_segments(path_MLS_plot, out_path_crowns_MLS, out_path)

## Test detect_stem and plot
stem_points <- detect_stem(path_MLS_ITD, stem_only = TRUE, all_attr = FALSE)
lidR::plot(stem_points)

## Get a segment and save.
segment <- lidR::filter_poi(stem_points, Z >= 1.3 & Z <= 1.5)
lidR::writeLAS(segment, file.path(out_path, "segment.laz"))
segment <- lidR::readLAS(file.path(out_path, "segment.laz"))

## Get segment diameter.
diameter <- segment_diameter(segment, n_points = 30, quantiles = c(0.05, 0.95))
print(diameter)

## Get diameter at 1.4m height from a LAS PC.
diameter_at_fixed_ht <- fixedht_diameter_LAS(stem_points, fixed_ht = 1.4,
                                             seg_ht = 0.25)
print(diameter_at_fixed_ht)

## Create a diameter profile.
d_profile <- stem_diameter_profile(stem_points, seg_ht = 0.2)
print(d_profile)

## Get diameter at 1.4m height from diameter profile.
diameter_at_fixed_ht <- fixedht_diameter_DF(d_profile, fixed_ht = 1.4)
print(diameter_at_fixed_ht)

## Create a node profile.
n_profile <- stem_node_profile(d_profile, tol = 1, base = 4,
                               split = c(FALSE, "Node"))
print(n_profile)

## Get variable dbh height.
dbh_ht <- get_dbh_ht(n_profile, bh = 1.4, min_ht = 0.2, max_ht = 2.9,
                     tol1 = 0.1, tol2 = 0.3)
print(dbh_ht)

## Get dbh at variable height.
dbh <- variable_ht_dbh(d_profile, dbh_ht)
print(dbh)

## Get stem volume.
tree_volume <- stem_volume(d_profile)
print(tree_volume)

## Get tree metrics MLS.
metrics_MLS_PC <- get_tree_metrics_TLS(path_MLS_ITD, height = TRUE,
                                       crown = TRUE, min_ht = 2)
print(metrics_MLS_PC)

# Get tree metrics ULS.
metrics_MLS_CHM <- get_tree_metrics_CHM(path_MLS_chm, out_path_crowns_MLS)
print(metrics_MLS_CHM)

# Clip with buffer from tree stem.
stem_class_las <- detect_stem(path_MLS_ITD, stem_only = FALSE, all_attr = TRUE)
clipped_las <- clip_with_buffer(stem_class_las, buffer = 0.5)
lidR::plot(clipped_las)

################################################################################
# Process ULS point cloud.
# Delineate crowns.
crowns_ULS <- delineate_crowns(path_ULS_chm, path_stem_map)
out_path_crowns_ULS <- file.path(out_path, "crowns_ULS.shp")
sf::write_sf(crowns_ULS, out_path_crowns_ULS)

ggplot2::ggplot(data = crowns_ULS) +
  ggplot2::geom_sf() +
  ggplot2::geom_sf_text(ggplot2::aes(label = treeID),
                        size = 3, color = "blue") +
  ggplot2::theme_bw() +
  ggplot2::labs(title = "ULS crowns with tree IDs")

# Create density rasters.
raster_density_ULS <- create_density_raster(path_ULS_plot)
raster::plot(raster_density_ULS$raster, main = "ULS point density raster")
point_density_ULS <-  raster_density_ULS$mean_density +
  raster_density_ULS$sd_density
print(point_density_ULS)

# Homogenise the point cloud.
homogenised_ULS <- homogenise_LAS(path_ULS_plot, out_path, point_density_ULS)
lidR::plot(lidR::readLAS(homogenised_ULS))

# Clip individual tree point clouds.
ITs_ULS <- clip_IT_segments(path_ULS_plot, out_path_crowns_ULS, out_path)

# Get tree metrics from the ULS point cloud.
tree_metrics_ULS_PC <- tibble::tibble()
for (tree in ITs_ULS){
  print(tree)
  tree_metrics_ULS_PC <- get_tree_metrics_ALS(tree, min_ht = 2,
                                         height = TRUE,
                                         crown = TRUE) %>%
    dplyr::select(treeID = file_id, dplyr::everything()) %>%
    bind_rows(., tree_metrics_ULS_PC)
}
print(tree_metrics_ULS_PC)
write.csv(tree_metrics_ULS_PC, file.path(getwd(), "Test_Outputs",
                                         "ULS_pointcloud_metrics_Test.csv"),
          row.names = FALSE)

# Get tree metrics from the ULS CHM.
tree_metrics_ULS_CHM <- get_tree_metrics_CHM(path_ULS_chm, out_path_crowns_ULS)
print(tree_metrics_ULS_CHM)
write.csv(tree_metrics_ULS_CHM, file.path(getwd(), "Test_Outputs",
                                         "ULS_pointcloud_metrics_Test.csv"),
          row.names = FALSE)

################################################################################
# End of Script
