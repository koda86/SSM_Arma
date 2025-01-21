remove_plane <- function(point_cloud, num_iterations = 1000, distance_threshold = 5, 
                         inlier_ratio_threshold = 0.5, alpha = 2){
  
  # **Install and Load Required Packages**
  required_packages <- c("geometry", "alphashape3d", "rgl")
  
  for(pkg in required_packages){
    if(!requireNamespace(pkg, quietly = TRUE)){
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  # **Data Validation and Preprocessing**
  
  # Check if required columns exist and are numeric
  required_cols <- c("X", "Y", "Z")
  if(!all(required_cols %in% colnames(point_cloud))){
    stop("Point cloud must contain columns named 'X', 'Y', and 'Z'.")
  }
  
  # Ensure columns are numeric
  if(!all(sapply(point_cloud[, required_cols], is.numeric))){
    stop("Columns 'X', 'Y', and 'Z' must be numeric.")
  }
  
  # Check for missing values
  if(any(is.na(point_cloud[, required_cols]))){
    stop("Point cloud contains NA values in 'X', 'Y', or 'Z' columns.")
  }
  
  # Check if enough points to compute concave hull
  if(nrow(point_cloud) < 4){
    stop("Not enough points to compute a concave hull.")
  }
  
  # **Remove Duplicate Points**
  
  # Convert to matrix for ashape3d
  point_matrix <- as.matrix(point_cloud[, required_cols])
  
  # Remove duplicate points
  point_matrix_unique <- unique(point_matrix)
  
  # Check the number of points before and after removing duplicates
  initial_count <- nrow(point_matrix)
  final_count <- nrow(point_matrix_unique)
  cat("Removed", initial_count - final_count, "duplicate points.\n")
  
  # Ensure there are still enough points after removing duplicates
  if(final_count < 4){
    stop("Not enough unique points to compute a concave hull after removing duplicates.")
  }
  
  # **Define Helper Functions**
  
  # Cross product of two vectors
  cross_product <- function(u, v) {
    return(c(
      u[2] * v[3] - u[3] * v[2],
      u[3] * v[1] - u[1] * v[3],
      u[1] * v[2] - u[2] * v[1]
    ))
  }
  
  # Plane fitting using three points
  fit_plane <- function(p1, p2, p3) {
    v1 <- p2 - p1
    v2 <- p3 - p1
    normal <- cross_product(v1, v2)
    
    if (all(abs(normal) < 1e-6)) {
      return(NULL)  # Points are collinear
    }
    
    normal <- normal / sqrt(sum(normal^2))  # Normalize
    d <- -sum(normal * p1)
    return(list(normal = normal, d = d))
  }
  
  # Compute distances from points to a plane
  point_plane_distance <- function(points, plane) {
    # points is a matrix with columns X, Y, Z
    normal <- plane$normal
    d <- plane$d
    distances <- abs(points %*% normal + d) / sqrt(sum(normal^2))
    return(as.vector(distances))
  }
  
  # **Define Vertical Axis and Angle Threshold**
  
  vertical_axis <- c(1, 0, 0)  # Adjust based on your coordinate system
  max_angle_degrees <- 10
  max_angle_radians <- max_angle_degrees * (pi / 180)
  
  # **Initialize Variables for RANSAC**
  
  best_plane <- NULL
  max_inliers <- 0
  best_inlier_indices <- NULL
  
  set.seed(123)  # For reproducibility
  
  ### **Compute Concave Hull Using Alpha Shape**
  
  ashape <- tryCatch({
    ashape3d(point_matrix_unique, alpha = alpha)
  }, error = function(e){
    stop("Error in computing concave hull: ", e$message)
  })
  
  # Extract hull point indices
  concave_hull_indices_unique <- unique(ashape$ia)
  
  if(length(concave_hull_indices_unique) < 3){
    stop("Concave hull could not be computed with the given alpha. Consider adjusting the alpha parameter.")
  }
  
  # Extract hull facets (triangles)
  concave_facets <- ashape$ih
  
  if(is.null(concave_facets)){
    stop("No facets found in the concave hull. The alpha parameter might be too low.")
  }
  
  # Map hull indices back to original unique point matrix
  hull_points_matrix <- point_matrix_unique[concave_hull_indices_unique, ]
  
  # Define non-hull points
  non_hull_indices_unique <- setdiff(1:nrow(point_matrix_unique), concave_hull_indices_unique)
  
  # Check if there are non-hull points to process
  if(length(non_hull_indices_unique) < 3){
    stop("Not enough non-hull points to perform plane fitting.")
  }
  
  ### **RANSAC Loop to Detect the Best Plane**
  
  for (i in 1:num_iterations) {
    # Randomly select 3 unique non-hull points
    sample_indices <- sample(non_hull_indices_unique, 3, replace = FALSE)
    p1 <- point_matrix_unique[sample_indices[1], ]
    p2 <- point_matrix_unique[sample_indices[2], ]
    p3 <- point_matrix_unique[sample_indices[3], ]
    
    # Fit plane
    plane <- fit_plane(p1, p2, p3)
    
    # Continue to next iteration if plane fitting failed (collinear points)
    if (is.null(plane)) {
      next
    }
    
    # Calculate distances of all non-hull points to the plane
    distances <- point_plane_distance(point_matrix_unique[non_hull_indices_unique, ], plane)
    
    # Identify inliers among non-hull points
    inlier_indices_non_hull <- which(distances < distance_threshold)
    num_inliers <- length(inlier_indices_non_hull)
    
    # Compute angle between plane normal and vertical axis
    dot_product <- sum(plane$normal * vertical_axis)
    angle <- acos(abs(dot_product))  # Absolute to account for normal direction
    
    # Check if angle is within the threshold
    if (angle > max_angle_radians) {
      next  # Skip this plane as it's not aligned with vertical axis
    }
    
    # Update best plane if current one has more inliers
    if (num_inliers > max_inliers) {
      max_inliers <- num_inliers
      best_plane <- plane
      # Map inlier indices back to unique point matrix indices
      best_inlier_indices_unique <- non_hull_indices_unique[inlier_indices_non_hull]
      best_inlier_indices <- best_inlier_indices_unique
    }
    
    # Early exit if a good enough plane is found
    if ((max_inliers / length(non_hull_indices_unique)) > inlier_ratio_threshold) {
      break
    }
  }
  
  # **Check if a Plane Was Found**
  
  if (!is.null(best_plane)) {
    cat("Plane detected with", max_inliers, "inliers out of", nrow(point_matrix_unique), "points.\n")
    
    ### **Extract Inlier and Outlier Points Among Non-Hull Points**
    inliers_unique <- point_matrix_unique[best_inlier_indices, ]
    outliers_unique <- point_matrix_unique[setdiff(non_hull_indices_unique, best_inlier_indices), ]
    
    ### **Integrate Hull Points Back**
    hull_points_unique <- hull_points_matrix  # Already a matrix
    
    ### **Combine Outliers with Hull Points to Form the Cleaned Point Cloud**
    clean_point_matrix <- rbind(outliers_unique, hull_points_unique)
    
    # Convert back to data frame for consistency
    clean_point_cloud <- as.data.frame(clean_point_matrix)
    colnames(clean_point_cloud) <- c("X", "Y", "Z")
    
    ### **Control Plot: Visualize Concave Hull Points**
    # Open a new 3D plotting window
    open3d()
    
    # Plot all cleaned points in light grey
    plot3d(clean_point_cloud$X, clean_point_cloud$Y, clean_point_cloud$Z, 
           col = "lightgrey", size = 1, alpha = 0.5, 
           xlab = "X", ylab = "Y", zlab = "Z",
           main = "Cleaned Point Cloud with Concave Hull")
    
    # Highlight concave hull points in red
    points3d(clean_point_cloud$X[(nrow(outliers_unique) + 1):nrow(clean_point_cloud)], 
             clean_point_cloud$Y[(nrow(outliers_unique) + 1):nrow(clean_point_cloud)], 
             clean_point_cloud$Z[(nrow(outliers_unique) + 1):nrow(clean_point_cloud)], 
             col = "red", size = 3)
    
    # Draw the concave hull facets as a mesh in semi-transparent blue
    for(i in 1:nrow(concave_facets)){
      idx <- concave_facets[i, ]
      if(length(idx) != 3){
        warning(paste("Facet", i, "does not have 3 points. Skipping this facet."))
        next
      }
      # Adjust indices to match clean_point_matrix
      triangle <- as.matrix(clean_point_matrix[idx, c("X", "Y", "Z")])
      triangles3d(triangle, color = "blue", alpha = 0.3, lit = FALSE)
    }
    
    # Add a legend for clarity
    legend3d("topright", legend = c("Cleaned Points", "Concave Hull Points", "Concave Hull Facets"), 
             pch = c(16, 16, NA), 
             col = c("lightgrey", "red", "blue"), 
             cex = 1.2, 
             inset = c(0.02, 0.02),
             fill = c(NA, NA, "blue"))
    
    ### **Optionally, Visualize the Detected Plane**
    # Define grid range based on point cloud bounds
    x_range <- range(clean_point_matrix[, "X"])
    z_range <- range(clean_point_matrix[, "Z"])  # Assuming Y is vertical
    
    grid_x <- seq(x_range[1], x_range[2], length.out = 10)
    grid_z <- seq(z_range[1], z_range[2], length.out = 10)
    
    # Create a meshgrid for X and Z
    mesh_grid <- expand.grid(X = grid_x, Z = grid_z)
    
    # Calculate Y based on the plane equation: ax + by + cz + d = 0 => y = (-ax - cz - d)/b
    if (abs(best_plane$normal[2]) > 1e-6) {  # Avoid division by zero
      mesh_grid$Y <- (-best_plane$normal[1] * mesh_grid$X - best_plane$normal[3] * mesh_grid$Z - best_plane$d) / best_plane$normal[2]
      
      # Plot the plane as a surface
      surface3d(mesh_grid$X, mesh_grid$Z, mesh_grid$Y, color = "green", alpha = 0.5, front = "lines")
    } else {
      warning("Plane normal's Y component is too small; cannot plot the plane.")
    }
    
    # **Return the Cleaned Point Cloud**
    return(clean_point_cloud)
    
  } else {
    cat("No plane detected aligned with the vertical axis.\n")
    return(as.data.frame(point_matrix_unique))  # Return unique point cloud as data frame
  }
}
