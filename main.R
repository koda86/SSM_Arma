# Load packages
library(Rvcg)
library(Morpho)
library(rgl)
library(geometry)
library(alphashape3d)
library(dplyr)
library(shapes)
library(RANN)

# 1. Data import and visualization ---------------------------------------------

ply_dir <- "~/tmp/dummdata_armasuisse"
ply_files <- list.files(ply_dir, pattern = "\\.ply$", full.names = TRUE)

file <- ply_files[3]

mesh <- Rvcg::vcgPlyRead(file)

v1 <- c(0,150,0) #Laenge wo abgeschnitten wird
norm <- c(0,1,0) #Normalenvktor der Ebene
mesh_cut <- Morpho::cutMeshPlane(mesh, v1, normal = norm, keep.upper = FALSE)

vertices <- mesh_cut$vb[1:3, ] |> t()  # Transpose to get Nx3 matrix

# Convert to data frame for easier handling
points_df <- data.frame(X = vertices[, 1], Y = vertices[, 2], Z = vertices[, 3])

# Visualisierung
x_limits <- c(-300, 300)
y_limits <- c(-300, 300)
z_limits <- c(-300, 300)

rgl::plot3d(points_df, col = "blue", size = 3,
            xlim = x_limits, ylim = y_limits, zlim = z_limits,
            xlab = "X", ylab = "Y", zlab = "Z")

# Downsample (Reduce Number of Vertices)
mesh_cut_downsampled <- Rvcg::vcgUniformRemesh(mesh_cut,
                               voxelSize = NULL,
                               offset = 0,
                               # mergeClost = TRUE
)

vertices <- mesh_cut_downsampled$vb[1:3, ] |> t()  # Transpose to get Nx3 matrix

# Convert to data frame for easier handling
points_df <- data.frame(X = vertices[, 1], Y = vertices[, 2], Z = vertices[, 3])

rgl::plot3d(points_df, col = "blue", size = 3,
            xlim = x_limits, ylim = y_limits, zlim = z_limits,
            xlab = "X", ylab = "Y", zlab = "Z")

# 1. Detect the plane (using RANSAC) and remove it -----------------------------
source("~/SSM_Arma/detect_plane.R")

clean_point_cloud <- remove_plane(points_df)

rgl::plot3d(clean_point_cloud, col = "black", size = 3,
            xlim = x_limits, ylim = y_limits, zlim = z_limits,
            xlab = "X", ylab = "Y", zlab = "Z")


# # Cut off the upper part of the foot to ensure all scans have the same height --
# # Vertical axis in these scans is the y axis!
# y_cutoff <- 150  # Adjust this value as needed
# 
# # Create vertices for a large plane at y = y_cutoff
# # The plane is a square spanning from -1000 to 1000 in x and z for ample coverage
# vertices <- matrix(c(
#   -300, y_cutoff, -300,  # Vertex 1
#   300, y_cutoff, -30,  # Vertex 2
#   300, y_cutoff,  300,  # Vertex 3
#   -300, y_cutoff,  300   # Vertex 4
# ), nrow = 3, byrow = FALSE)
# 
# vertices <- matrix(c(
#   -300, y_cutoff, -300,  # Vertex 1
#   -300, y_cutoff, 300,  # Vertex 2
#   -300, y_cutoff,  300,  # Vertex 3
#   -300, y_cutoff,  300   # Vertex 4
# ), nrow = 3, byrow = FALSE)
# 
# # Define indices for two triangular faces forming the plane
# indices <- matrix(c(
#   1, 2, 3,  # Triangle 1
#   1, 3, 4   # Triangle 2
# ), nrow = 3, byrow = FALSE)
# 
# # Create the mesh3d object representing the clipping plane
# clipping_plane <- tmesh3d(
#   vertices = vertices,
#   indices = indices,
#   homogeneous = FALSE,
#   material = list(color = "red")
# )
# 
# # Verify the class
# class(clipping_plane)  # Should return "mesh3d"
# 
# # Define the offset to position the plane at y_cutoff
# # The plane equation is: 0*x + 1*y + 0*z + d = 0 => y + d = 0 => d = -y_cutoff
# # cut_mesh <- vcgClost(clipping_plane, mesh)
# 
# cut_mesh <- vcgClost(
#   x = clipping_plane,
#   mesh = mesh,
#   sign = TRUE,         # Determines the side to keep based on the plane normal
#   keep = "below"       # Options: "below" or "above"
# )
# 
# 
# # After cutting, it's good practice to recompute normals and clean the mesh to ensure its integrity.
# cut_mesh <- vcgUpdateNormals(cut_mesh)
# cut_mesh <- vcgClean(cut_mesh)
# 
# vertices <- cut_mesh$vb[1:3, ] |> t()  # Transpose to get Nx3 matrix
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


vertical_axis <- 2

z_coords <- mesh$vb[vertical_axis, ] # Extract y-coordinates of all vertices
z_cutoff <- 150 # Alternative: Replace with your desired value
keep_vertices <- z_coords <= z_cutoff # Logical vector for vertices to keep

# Extract the vertex indices to keep
keep_indices <- which(keep_vertices)

# Subset the vertices
new_vertices <- mesh$vb[, keep_indices]

# Extract the original faces (each row represents a triangle)
original_faces <- mesh$it

# Identify faces where all three vertices are kept
faces_to_keep <- apply(original_faces, 1, function(face) all(keep_vertices[face]))

# Subset the faces
new_faces <- original_faces[faces_to_keep, ]

# Create a mapping from old vertex indices to new indices
# Initialize all mappings to zero
old_to_new <- integer(length(keep_vertices))
# Assign new indices to kept vertices
old_to_new[keep_indices] <- seq_along(keep_indices)

# Update the face indices to match the new vertex subset
new_faces <- old_to_new[new_faces]

# Create a new mesh object
new_mesh <- mesh

# Update the vertices
new_mesh$vb <- new_vertices

# Update the faces
new_mesh$it <- new_faces

# Optional: Recompute normals and clean the mesh
new_mesh <- vcgClean(new_mesh) # , sel = "mesh"


# Downsample (Reduce Number of Vertices)
mesh <- Rvcg::vcgUniformRemesh(new_mesh,
                               voxelSize = NULL,
                               offset = 0,
                               # mergeClost = TRUE
)





vertices <- new_mesh$vb[1:3, ] |> t()  # Transpose to get Nx3 matrix

# Convert to data frame for easier handling
points_df <- data.frame(X = vertices[, 1], Y = vertices[, 2], Z = vertices[, 3])

# Visualisierung
x_limits <- c(-300, 300)
y_limits <- c(-300, 300)
z_limits <- c(-300, 300)

rgl::plot3d(points_df, col = "blue", size = 3,
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



