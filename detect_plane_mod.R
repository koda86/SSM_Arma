remove_plane <- function(point_cloud, num_iterations, distance_threshold, inlier_ratio_threshold){
  
  # Load required libraries inside the function
  if(!requireNamespace("geometry", quietly = TRUE)) {
    install.packages("geometry")
  }
  library(geometry)
  
  # Define a function to compute the cross product of two 3D vectors
  cross_product <- function(u, v) {
    return(c(
      u[2] * v[3] - u[3] * v[2],
      u[3] * v[1] - u[1] * v[3],
      u[1] * v[2] - u[2] * v[1]
    ))
  }
  
  # Plane fitting function using the manual cross product
  fit_plane <- function(p1, p2, p3) {
    # Create vectors lying on the plane
    v1 <- p2 - p1
    v2 <- p3 - p1
    
    # Compute the normal vector using the cross product
    normal <- cross_product(v1, v2)
    
    # Check if the points are collinear (i.e., cross product is zero vector)
    if (all(abs(normal) < 1e-6)) {
      return(NULL)  # Collinear points; no unique plane
    }
    
    # Normalize the normal vector
    normal <- normal / sqrt(sum(normal^2))
    
    # Calculate the plane offset 'd'
    d <- -sum(normal * p1)
    
    return(list(normal = normal, d = d))
  }
  
  # Function to compute distances from points to a plane
  point_plane_distance <- function(points, plane) {
    # Ensure 'points' has columns X, Y, Z
    if (!all(c("X", "Y", "Z") %in% colnames(points))) {
      stop("points must have columns named 'X', 'Y', and 'Z'")
    }
    
    # Ensure plane has 'normal' and 'd'
    if (!("normal" %in% names(plane) && "d" %in% names(plane))) {
      stop("plane must be a list with 'normal' and 'd'")
    }
    
    # Ensure plane$normal is a numeric vector of length 3
    if (!is.numeric(plane$normal) || length(plane$normal) != 3) {
      stop("plane$normal must be a numeric vector of length 3")
    }
    
    # Convert points to matrix
    points_mat <- as.matrix(points[, c("X", "Y", "Z")])
    
    # Extract normal and d
    normal <- plane$normal
    d <- plane$d
    
    # Compute distances
    numerators <- abs(points_mat %*% normal + d)
    denominator <- sqrt(sum(normal^2))  # Should be 1 since normal is normalized
    distances <- numerators / denominator
    
    return(as.vector(distances))
  }
  
  # Define the vertical axis (Y-axis in this example)
  vertical_axis <- c(1, 0, 0)  # Adjust if your vertical axis is different
  
  # Define the maximum angle (in degrees) between plane normal and vertical axis
  max_angle_degrees <- 10
  max_angle_radians <- max_angle_degrees * (pi / 180)
  
  # Initialize variables to store the best plane
  best_plane <- NULL
  max_inliers <- 0
  best_inlier_indices <- NULL
  
  set.seed(123)  # For reproducibility
  
  ### Identify Hull Points (added 20.01.2025)
  
  # Compute the convex hull of the point cloud
  hull <- convhulln(point_cloud[, c("X", "Y", "Z")], options = "FA")
  
  # Extract hull point indices correctly
  hull_indices <- unique(as.vector(hull$hull))
  
  # Define non-hull points
  non_hull_indices <- setdiff(1:nrow(point_cloud), hull_indices)
  
  # Check if there are non-hull points to process
  if(length(non_hull_indices) < 3){
    stop("Not enough non-hull points to perform plane fitting.")
  }
  
  # RANSAC Loop to detect the best plane
  for (i in 1:num_iterations) {
    # Randomly select 3 unique non-hull points
    sample_indices <- sample(non_hull_indices, 3, replace = FALSE)
    p1 <- as.numeric(point_cloud[sample_indices[1], c("X", "Y", "Z")])
    p2 <- as.numeric(point_cloud[sample_indices[2], c("X", "Y", "Z")])
    p3 <- as.numeric(point_cloud[sample_indices[3], c("X", "Y", "Z")])
    
    # Fit plane
    plane <- fit_plane(p1, p2, p3)
    
    # Continue to next iteration if plane fitting failed (collinear points)
    if (is.null(plane)) {
      next
    }
    
    # Calculate distances of all non-hull points to the plane
    distances <- point_plane_distance(point_cloud[non_hull_indices, ], plane)
    
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
      # Map inlier indices back to original point cloud indices
      best_inlier_indices <- non_hull_indices[inlier_indices_non_hull]
    }
    
    # Early exit if a good enough plane is found
    if ((max_inliers / length(non_hull_indices)) > inlier_ratio_threshold) {
      break
    }
  }
  
  # Check if a plane was found
  if (!is.null(best_plane)) {
    cat("Plane detected with", max_inliers, "inliers out of", nrow(point_cloud), "points.\n")
    
    # Extract inlier and outlier points among non-hull points
    inliers <- point_cloud[best_inlier_indices, ]
    outliers_non_hull <- point_cloud[setdiff(non_hull_indices, best_inlier_indices), ]
    
    # **Integrate Hull Points Back**
    hull_points <- point_cloud[hull_indices, ]
    
    # Combine outliers with hull points to form the cleaned point cloud
    clean_point_cloud <- rbind(outliers_non_hull, hull_points)
    
    # Optionally, visualize the detected plane and hull
    # Uncomment the following lines to visualize
    # library(rgl)
    # plot3d(outliers_non_hull$X, outliers_non_hull$Y, outliers_non_hull$Z, col = "blue", size = 1,
    #        main = "Plane Detection with RANSAC (Aligned with Vertical Axis)",
    #        xlab = "X", ylab = "Y", zlab = "Z")
    # points3d(inliers$X, inliers$Y, inliers$Z, col = "red", size = 3)
    # points3d(hull_points$X, hull_points$Y, hull_points$Z, col = "green", size = 3)
    
    # Optionally, visualize the plane
    # Create a grid of points on the plane for visualization
    plane_normal <- best_plane$normal
    plane_d <- best_plane$d
    
    # Define grid range based on point cloud bounds
    x_range <- range(point_cloud$X)
    z_range <- range(point_cloud$Z)  # Assuming Y is vertical
    
    grid_x <- seq(x_range[1], x_range[2], length.out = 10)
    grid_z <- seq(z_range[1], z_range[2], length.out = 10)
    
    # Create a meshgrid for X and Z
    mesh_grid <- expand.grid(X = grid_x, Z = grid_z)
    
    # Calculate Y based on the plane equation: ax + by + cz + d = 0 => y = (-ax - cz - d)/b
    if (abs(plane_normal[2]) > 1e-6) {  # Avoid division by zero
      mesh_grid$Y <- (-plane_normal[1] * mesh_grid$X - plane_normal[3] * mesh_grid$Z - plane_d) / plane_normal[2]
      
      # Plot the plane
      # Uncomment the following line to visualize the plane
      # surface3d(mesh_grid$X, mesh_grid$Z, mesh_grid$Y, color = "green", alpha = 0.5, front = "lines")
    } else {
      warning("Plane normal's Y component is too small; cannot plot the plane.")
    }
    
    # Return the cleaned point cloud
    return(clean_point_cloud)
    
  } else {
    cat("No plane detected aligned with the vertical axis.\n")
    return(point_cloud)  # Return original point cloud if no plane is detected
  }
}
