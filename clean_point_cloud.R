# Load packages
library(rgl)
library(geometry)
library(alphashape3d)
library(dplyr)
library(Rvcg)
library(shapes)
library(Morpho)
library(RANN)

# 1. Data Import and Visualization ---------------------------------------------

ply_dir <- "~/tmp/dummdata_armasuisse"
ply_files <- list.files(ply_dir, pattern = "\\.ply$", full.names = TRUE)

file <- ply_files[3]
mesh <- vcgPlyRead(file)
mesh <- vcgUniformRemesh(mesh, voxelSize = NULL, offset = 0) # Downsample Shapes (Reduce Number of Vertices)

vertices <- mesh$vb[1:3, ] |> t()  # Transpose to get Nx3 matrix

# Convert to data frame for easier handling
points_df <- data.frame(X = vertices[,1], Y = vertices[,2], Z = vertices[,3])

# Visualisierung
x_limits <- c(-300, 300)
y_limits <- c(-300, 300)
z_limits <- c(-300, 300)

rgl::plot3d(points_df, col = "blue", size = 3,
       xlim = x_limits, ylim = y_limits, zlim = z_limits,
       xlab = "X", ylab = "Y", zlab = "Z")

# 1. Detect the plane (using RANSAC) and remove it -----------------------------
source("~/Nextcloud/project-fab-forschung/Publikationen/ISB/Koska_2025/R/detect_plane.R")

clean_point_cloud <- remove_plane(points_df)

rgl::plot3d(clean_point_cloud, col = "black", size = 3,
            xlim = x_limits, ylim = y_limits, zlim = z_limits,
            xlab = "X", ylab = "Y", zlab = "Z")


# Nächster Schritt: Den Fuss bei einer bestimmten Höhe abschneiden

# Danach: Ausreisser entfernen, downsamplen, Punkte gleichmäßig verteilen



# 2. Preprocessing: Noise and Outlier Removal ----------------------------------
# It's essential to clean the point cloud by removing noise and outliers to ensure accurate hull extraction.

# a. Statistical Outlier Removal
# While R doesn't have a direct equivalent to PCL's Statistical Outlier Removal (SOR), you can implement a similar approach using k-nearest neighbors.

# Function to perform Statistical Outlier Removal
remove_outliers <- function(data, k = 20, sigma = 1.0) {
  # Find k nearest neighbors
  nn <- nn2(data = data[, c("X", "Y", "Z")], query = data[, c("X", "Y", "Z")], k = k + 1)

  # Exclude the point itself by removing the first neighbor (distance = 0)
  distances <- rowMeans(nn$nn.dists[, 2:(k + 1)])

  # Compute mean and standard deviation of distances
  mean_dist <- mean(distances)
  sd_dist <- sd(distances)

  # Keep points within mean + sigma * sd
  filtered_data <- data %>%
    mutate(avg_dist = distances) %>%
    filter(avg_dist < (mean_dist + sigma * sd_dist)) %>%
    select(-avg_dist)

  return(filtered_data)
}

# Apply outlier removal
points_filtered <- remove_outliers(points_df, k = 10, sigma = 0.0001)

# Visualize filtered point cloud
rgl::plot3d(points_filtered, col = "blue", size = 3,
            xlim = x_limits, ylim = y_limits, zlim = z_limits,
            xlab = "X", ylab = "Y", zlab = "Z")

# plot3d(points_filtered$X, points_filtered$Y, points_filtered$Z, size = 3, col = "black",
#        main = "Filtered Point Cloud")

# # c. Voxel Grid Downsampling (Optional)
# # If further downsampling is needed to reduce computational load:
# 
# # Voxel Grid Downsampling
# voxel_size <- 0.005  # Adjust based on your data's scale
# 
# # Compute voxel indices
# points_filtered <- points_filtered %>%
#   mutate(vx = floor(X / voxel_size),
#          vy = floor(Y / voxel_size),
#          vz = floor(Z / voxel_size))
# 
# # Aggregate points within each voxel by taking the centroid
# points_downsampled <- points_filtered %>%
#   group_by(vx, vy, vz) %>%
#   summarize(X = mean(X), Y = mean(Y), Z = mean(Z)) %>%
#   ungroup()
# 
# # Visualize downsampled point cloud
# plot3d(points_downsampled$X, points_downsampled$Y, points_downsampled$Z, 
#        size = 1, col = "red", main = "Downsampled Point Cloud")


# 4. Surface Reconstruction to Extract the Outer Hull---------------------------

# Create an alpha shape
alpha_value <- 0.02  # Adjust based on data density and desired hull tightness
ashape <- ashape3d(as.matrix(points_df), alpha = alpha_value)
ashape <- ashape3d(as.matrix(points_filtered), alpha = alpha_value)

# Visualize the alpha shape
# plot(ashape, col = "cyan", alpha = 0.3, main = "Alpha Shape Hull")


# 5. Filtering to Retain Only the Outer Hull -----------------------------------

# Extract the hull mesh from the alpha shape
hull_mesh <- as.mesh3d(ashape)

# Extract vertices of the hull mesh
hull_vertices <- hull_mesh$vb[1:3, ] %>% t()

# Convert to data frame
hull_points <- data.frame(X = hull_vertices[,1], Y = hull_vertices[,2], Z = hull_vertices[,3])

# Visualize the hull points
rgl::plot3d(hull_points, col = "black", size = 3,
            xlim = x_limits, ylim = y_limits, zlim = z_limits,
            xlab = "X", ylab = "Y", zlab = "Z")



