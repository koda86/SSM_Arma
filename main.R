# TODO: Entfernung der Ebenen verbessern

# Load packages
library(Rvcg)
library(reticulate)
library(Morpho)
library(rgl)
library(geometry)
library(alphashape3d)
library(dplyr)
library(shapes)
library(RANN)

# Selbstgeschriebene Funktionen
source("~/SSM_Arma/detect_plane_mod.R") # Detect the plane (using RANSAC) and remove it
# source("~/SSM_Arma/detect_plane.R") # Alte Version ohne "Hüllenerkennung"
# source("~/SSM_Arma/remove_plane.R")
source("~/SSM_Arma/reconstruct_surface_reticulate.R") # Surface reconstruction

# 1. Data import and visualization ---------------------------------------------

# setwd("/run/user/1000/gvfs/smb-share:server=ex2100c.local,share=projektarchiv/Armasuisse/2014-12-11_Backup_AFS_Footscans/Armasuisse3D/konvertiert")
# ply_files <- list.files()

ply_dir <- "~/tmp/dummdata_armasuisse"
ply_files <- list.files(ply_dir, pattern = "\\.ply$", full.names = TRUE)

file <- ply_files[19] # Beispielfiles: 3 - mit Ebene, 10 - ohne Ebene im Fuß, 18 - viele Ebenen

mesh <- Rvcg::vcgPlyRead(file)

v1 <- c(0, 150, 0) # Cutoff threshold
norm <- c(0, 1, 0) # Normal plane
mesh_cut <- Morpho::cutMeshPlane(mesh, v1, normal = norm, keep.upper = FALSE)

vertices <- mesh_cut$vb[1:3, ] |> t() # Transpose to get Nx3 matrix

# Convert to data frame for easier handling
points_df <- data.frame(X = vertices[, 1], Y = vertices[, 2], Z = vertices[, 3])

# Visualisierung
x_limits <- c(-300, 300)
y_limits <- c(-300, 300)
z_limits <- c(-300, 300)

rgl::plot3d(points_df, col = "blue", size = 3,
            xlim = x_limits, ylim = y_limits, zlim = z_limits,
            xlab = "X", ylab = "Y", zlab = "Z")

# Detektionsschema schreiben um festzustellen, ob überhaupt eine Ebene im Fuß existiert
# Grundidee: Füße mit zusätzlichen Ebenen haben (deutlich) mehr Punkte ...
number_points <- max(dim(mesh_cut$vb))
cat(number_points)

threshold_planes_existance <- 8300
if (number_points > threshold_planes_existance) {
  # Hier ne Schleife reinbauen die drüberläuft bis die 3 dünnen Ebenen weg sind. Anschließend: Outlier entfernen!
  # clean_point_cloud <- remove_plane(points_df)
  points_df <- remove_plane(points_df,
                            num_iterations = 1000, # Number of iterations
                            distance_threshold = 5, # Adjust based on your data's scale
                            inlier_ratio_threshold = 0.5 # Minimum ratio of inliers to accept a plane
  )
  clean_point_cloud <- points_df
  
  rgl::plot3d(clean_point_cloud, col = "black", size = 3,
              xlim = x_limits, ylim = y_limits, zlim = z_limits,
              xlab = "X", ylab = "Y", zlab = "Z")
  
  
} else {
  
  reconstructed_mesh3d_object <- reconstruct_poisson(clean_point_cloud)
}

# Visualize the original and reconstructed meshes
open3d()
shade3d(reconstructed_mesh3d_object, color = "magenta", alpha = 0.7)

#  Wenn nicht, das Mesh direkt resamplen
# Downsample (Reduce Number of Vertices)
mesh_cut_downsampled <- Rvcg::vcgUniformRemesh(mesh_cut,
                                               voxelSize = 5,
                                               offset = 0,
                                               mergeClost = TRUE
)

vertices_downsampled <- mesh_cut_downsampled$vb[1:3, ] |> t()
points_df_downsampled <- data.frame(X = vertices_downsampled[, 1], Y = vertices_downsampled[, 2], Z = vertices_downsampled[, 3])

rgl::plot3d(points_df_downsampled, col = "blue", size = 3,
            xlim = x_limits, ylim = y_limits, zlim = z_limits,
            xlab = "X", ylab = "Y", zlab = "Z")




# pcd <- rgl::as.mesh3d(clean_point_cloud)
# resample_mesh <- vcgUniformRemesh(mesh, voxelSize = 0.06)

# Advanced front surface reconstruction
# package ‘RCGAL’ is not available for this version of R
# https://www.r-bloggers.com/2022/01/surface-reconstruction-with-rcgal/
# remotes::install_github(
#   "stla/RCGAL", dependencies = TRUE, build_opts = "--no-multiarch"
# )
library(RCGAL)
afs_mesh <- AFSreconstruction(SolidMobiusStrip_cloud)

open3d(windowRect = c(50, 50, 562, 562))
view3d(0, -50, zoom = 0.75)
shade3d(afs_mesh, color = "darkred")


# Export für Weiterverarbeitung in MeshLab -------------------------------------
write_ply <- function(df, filename) {
  n <- nrow(df)
  header <- c(
    "ply",
    "format ascii 1.0",
    paste("element vertex", n),
    "property float x",
    "property float y",
    "property float z",
    "end_header"
  )
  
  # Write header
  writeLines(header, con = filename)
  
  # Write vertex data
  write.table(df, file = filename, append = TRUE, sep = " ", 
              row.names = FALSE, col.names = FALSE)
}

# Replace 'clean_point_cloud' with your actual data frame name
write_ply(clean_point_cloud, "~/tmp/clean_point_cloud.ply")

#  Reimport reconstructed foot
# Replace with your actual file path
ply_file <- "~/tmp/reconstructed_foot.ply"

# Read the PLY file
mesh <- vcgImport(ply_file, updateNormals = TRUE)

vertices <- mesh$vb[1:3, ] |> t() # Transpose to get Nx3 matrix

# Convert to data frame for easier handling
points_df <- data.frame(X = vertices[, 1], Y = vertices[, 2], Z = vertices[, 3])

rgl::plot3d(points_df, col = "blue", size = 3,
            xlim = x_limits, ylim = y_limits, zlim = z_limits,
            xlab = "X", ylab = "Y", zlab = "Z")




# Downsample (Reduce Number of Vertices)
mesh_cut_downsampled <- Rvcg::vcgUniformRemesh(mesh_cut,
                               voxelSize = 30,
                               offset = 0,
                               # mergeClost = TRUE
)

#  Alternativer Downsamplingansatz: **Downsampling the Mesh to a Specific Number of Vertices**
# desired_vertices <- 1000  # Adjust this number as needed
mesh_cut_downsampled <- Rvcg::vcgQEdecim(mesh_cut, percent = 0.8)

vertices <- mesh_cut_downsampled$vb[1:3, ] |> t()  # Transpose to get Nx3 matrix

# Convert to data frame for easier handling
points_df <- data.frame(X = vertices[, 1], Y = vertices[, 2], Z = vertices[, 3])

rgl::plot3d(points_df, col = "blue", size = 3,
            xlim = x_limits, ylim = y_limits, zlim = z_limits,
            xlab = "X", ylab = "Y", zlab = "Z")


# 2. Preprocessing: Noise and Outlier Removal ----------------------------------
# It's essential to clean the point cloud by removing noise and outliers to ensure accurate hull extraction.

# a. Statistical Outlier Removal
# While R doesn't have a direct equivalent to PCL's Statistical Outlier Removal (SOR), you can implement a similar approach using k-nearest neighbors.

# Function to perform Statistical Outlier Removal
remove_outliers <- function(data, k = 20, sigma = 1.0) {
  # Find k nearest neighbors
  nn <- RANN::nn2(data = data[, c("X", "Y", "Z")], query = data[, c("X", "Y", "Z")], k = k + 1)
  
  # Exclude the point itself by removing the first neighbor (distance = 0)
  distances <- rowMeans(nn$nn.dists[, 2:(k + 1)])
  
  # Compute mean and standard deviation of distances
  mean_dist <- mean(distances)
  sd_dist <- sd(distances)
  
  # Keep points within mean + sigma * sd
  filtered_data <- data %>%
    mutate(avg_dist = distances) %>%
    filter(avg_dist < (mean_dist + sigma * sd_dist)) %>%
    dplyr::select(-avg_dist)
  
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



# vertical_axis <- 2
# 
# z_coords <- mesh$vb[vertical_axis, ] # Extract y-coordinates of all vertices
# z_cutoff <- 150 # Alternative: Replace with your desired value
# keep_vertices <- z_coords <= z_cutoff # Logical vector for vertices to keep
# 
# # Extract the vertex indices to keep
# keep_indices <- which(keep_vertices)
# 
# # Subset the vertices
# new_vertices <- mesh$vb[, keep_indices]
# 
# # Extract the original faces (each row represents a triangle)
# original_faces <- mesh$it
# 
# # Identify faces where all three vertices are kept
# faces_to_keep <- apply(original_faces, 1, function(face) all(keep_vertices[face]))
# 
# # Subset the faces
# new_faces <- original_faces[faces_to_keep, ]
# 
# # Create a mapping from old vertex indices to new indices
# # Initialize all mappings to zero
# old_to_new <- integer(length(keep_vertices))
# # Assign new indices to kept vertices
# old_to_new[keep_indices] <- seq_along(keep_indices)
# 
# # Update the face indices to match the new vertex subset
# new_faces <- old_to_new[new_faces]
# 
# # Create a new mesh object
# new_mesh <- mesh
# 
# # Update the vertices
# new_mesh$vb <- new_vertices
# 
# # Update the faces
# new_mesh$it <- new_faces
# 
# # Optional: Recompute normals and clean the mesh
# new_mesh <- vcgClean(new_mesh) # , sel = "mesh"
# 
# 
# # Downsample (Reduce Number of Vertices)
# mesh <- Rvcg::vcgUniformRemesh(new_mesh,
#                                voxelSize = NULL,
#                                offset = 0,
#                                # mergeClost = TRUE
# )
# 
# 
# 
# 
# 
# vertices <- new_mesh$vb[1:3, ] |> t()  # Transpose to get Nx3 matrix
# 
# # Convert to data frame for easier handling
# points_df <- data.frame(X = vertices[, 1], Y = vertices[, 2], Z = vertices[, 3])
# 
# # Visualisierung
# x_limits <- c(-300, 300)
# y_limits <- c(-300, 300)
# z_limits <- c(-300, 300)
# 
# rgl::plot3d(points_df, col = "blue", size = 3,
#             xlim = x_limits, ylim = y_limits, zlim = z_limits,
#             xlab = "X", ylab = "Y", zlab = "Z")





# Danach: Ausreisser entfernen, downsamplen, Punkte gleichmäßig verteilen






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



