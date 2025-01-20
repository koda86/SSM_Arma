remove_plane <- function(point_cloud){
  
  # Define a function to compute the cross product of two 3D vectors
  cross_product <- function(u, v) {
    return(c(
      u[2] * v[3] - u[3] * v[2],
      u[3] * v[1] - u[1] * v[3],
      u[1] * v[2] - u[2] * v[1]
    ))
  }
  
  # Corrected plane fitting function using the manual cross product
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
  
  point_plane_distance <- function(points, plane) {
    # Ensure 'points' is a data frame or matrix with columns X, Y, Z
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
    
    # Ensure 'points_mat' is N x 3
    if (ncol(points_mat) != 3) {
      stop("points matrix must have exactly 3 columns")
    }
    
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
  
  num_iterations <- 1000          # Number of iterations
  distance_threshold <- 5      # Adjust based on your data's scale
  inlier_ratio_threshold <- 0.5    # Minimum ratio of inliers to accept a plane
  
  # Initialize variables to store the best plane
  best_plane <- NULL
  max_inliers <- 0
  best_inlier_indices <- NULL
  
  set.seed(123)  # For reproducibility
  
  for (i in 1:num_iterations) {
    # Randomly select 3 unique points
    sample_indices <- sample(1:nrow(point_cloud), 3, replace = FALSE)
    p1 <- as.numeric(point_cloud[sample_indices[1], c("X", "Y", "Z")])
    p2 <- as.numeric(point_cloud[sample_indices[2], c("X", "Y", "Z")])
    p3 <- as.numeric(point_cloud[sample_indices[3], c("X", "Y", "Z")])
    
    # Fit plane
    plane <- fit_plane(p1, p2, p3)
    
    # Continue to next iteration if plane fitting failed (collinear points)
    if (is.null(plane)) {
      next
    }
    
    # Calculate distances of all points to the plane
    distances <- point_plane_distance(point_cloud, plane)
    
    # Identify inliers
    inlier_indices <- which(distances < distance_threshold)
    num_inliers <- length(inlier_indices)
    
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
      best_inlier_indices <- inlier_indices
    }
    
    # Early exit if a good enough plane is found
    if ((max_inliers / nrow(point_cloud)) > inlier_ratio_threshold) {
      break
    }
  }
  
  # Check if a plane was found
  if (!is.null(best_plane)) {
    cat("Plane detected with", max_inliers, "inliers out of", nrow(point_cloud), "points.\n")
    
    # Extract inlier and outlier points
    inliers <- point_cloud[best_inlier_indices, ]
    outliers <- point_cloud[-best_inlier_indices, ]
    
    # Visualizing the detected plane
    # plot3d(outliers$X, outliers$Y, outliers$Z, col = "blue", size = 1, 
    #        main = "Plane Detection with RANSAC (Aligned with Vertical Axis)", 
    #        xlab = "X", ylab = "Y", zlab = "Z")
    # points3d(inliers$X, inliers$Y, inliers$Z, col = "red", size = 3)
    
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
    mesh_grid$Y <- (-plane_normal[1] * mesh_grid$X - plane_normal[3] * mesh_grid$Z - plane_d) / plane_normal[2]
    
    # Plot the plane
    # surface3d(mesh_grid$X, mesh_grid$Z, mesh_grid$Y, color = "green", alpha = 0.5, front = "lines")
    
    # Remove the detected plane from the point cloud
    clean_point_cloud <- outliers
    
    # Visualize the cleaned point cloud
    # rgl::plot3d(clean_point_cloud, col = "black", size = 3,
    #             xlim = x_limits, ylim = y_limits, zlim = z_limits,
    #             xlab = "X", ylab = "Y", zlab = "Z")
    
  } else {
    cat("No plane detected aligned with the vertical axis.\n")
  }
  
  return(clean_point_cloud)
}











