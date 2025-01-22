# Function to identify planes using a ray-csting logic
# The idea is that a ray is cast through the foot along its anterior-posterior direction (slightly above the ground).
# Along this line, there should only be two intersections with the foot. Additional
# intersections represent addition planes (that need to be removed in the rest of the script).

library(rgl)

identify_planes <- function(point_cloud, distance_threshold) {
  
  # Define the ray
  # x - anterior-posteriore Achse
  # y - vertikale Achse
  ray_start <- c(-100, 50, -40)   # Start slightly below the foot
  ray_end <- c(500, 50, -40)      # End slightly above the foot
  
  # Ray-casting logic
  count_intersections <- function(ray_start, ray_end, point_cloud) {
    
    points <- as.matrix(point_cloud)
    
    # Direction of the ray
    ray_dir <- ray_end - ray_start
    
    # Compute parameter t for each point using the ray equation: ray_start + t * ray_dir
    t_values <- (points - matrix(ray_start, nrow = nrow(points), ncol = 3, byrow = TRUE)) %*% ray_dir
    t_values <- t_values / sum(ray_dir^2) # Normalize by ray_dir length squared
    
    # Points on the ray are valid if t is positive and within the ray segment (0 to 1)
    valid_points <- points[t_values >= 0 & t_values <= 1, ]
    
    # Count intersections by finding points very close to the ray
    distances <- apply(valid_points, 1, function(point) {
      v <- point - ray_start
      cross_prod <- c(
        ray_dir[2] * v[3] - ray_dir[3] * v[2], # x component
        ray_dir[3] * v[1] - ray_dir[1] * v[3], # y component
        ray_dir[1] * v[2] - ray_dir[2] * v[1]  # z component
      )
      sqrt(sum(cross_prod^2)) / sqrt(sum(ray_dir^2)) # Perpendicular distance
    })
    
    # Points within a small threshold distance are considered intersections
    sum(distances < distance_threshold)
  }
  
  # Count intersections
  intersections <- count_intersections(ray_start, ray_end, point_cloud)
  # cat("Number of intersections:", intersections, "\n")
  
  # Highlight intersected points (for debugging/visualization)
  if (intersections > 0) {
    intersected_points <- point_cloud[apply(point_cloud, 1, function(p) {
      # Vector from ray_start to point p
      v <- p - ray_start
      ray_dir <- ray_end - ray_start
      ray_dir <- ray_dir / sqrt(sum(ray_dir^2)) # Normalize ray direction
      
      # Compute the cross product for perpendicular distance
      cross_prod <- c(
        ray_dir[2] * v[3] - ray_dir[3] * v[2], # x component
        ray_dir[3] * v[1] - ray_dir[1] * v[3], # y component
        ray_dir[1] * v[2] - ray_dir[2] * v[1]  # z component
      )
      dist_to_ray <- sqrt(sum(cross_prod^2)) / sqrt(sum(ray_dir^2)) # Perpendicular distance
      
      # Check if the point is within the threshold
      return(dist_to_ray < distance_threshold)
    }), ]
    
    # # Visualize the point cloud, ray, and intersected points
    # plot3d(point_cloud, col = "blue", size = 3, xlab = "X", ylab = "Y", zlab = "Z")
    # lines3d(rbind(ray_start, ray_end), col = "red", lwd = 3) # Ray
    # points3d(intersected_points, col = "green", size = 15) # Intersected points
  }
  
  return(intersected_points)
}



