# library(reticulate)
# library(Rvcg)
# library(rgl)

reconstruct_poisson <- function(point_cloud) {
  
  # Create a virtual environment
  use_condaenv("r-open3d-env", required = TRUE)
  
  # Install Open3D in the "r-open3d-env" Conda environment
  # conda_install(
  #   envname = "r-open3d-env",
  #   packages = "open3d",
  #   pip = FALSE  # Use Conda for installation
  # )
  
  # Import Open3D and check its version
  o3d <- import("open3d")
  # open3d_version <- o3d$`__version__`
  # cat("Open3D version:", open3d_version, "\n")
  
  # Falls das Objekt nicht als mesh3d vorliegt erst umwandeln
  point_cloud_mesh3d <- mesh3d(point_cloud)
  # In dem Fall fehlen die normals, die zunächst noch bestimmt werden müssen
  point_cloud_mesh3d <- vcgUpdateNormals(point_cloud_mesh3d) # Compute normals for the mesh3d object
  
  # Extract vertices and normals
  vertices <- t(point_cloud_mesh3d$vb[1:3, ])
  normals <- t(point_cloud_mesh3d$normals[1:3, ])
  
  # Create Open3D PointCloud object
  pcd <- o3d$geometry$PointCloud()
  pcd$points <- o3d$utility$Vector3dVector(vertices)
  pcd$normals <- o3d$utility$Vector3dVector(normals)
  
  # Perform Poisson Surface Reconstruction
  # Explicitly cast the depth parameter to an integer
  depth_int <- as.integer(9)
  
  result <- o3d$geometry$TriangleMesh$create_from_point_cloud_poisson(
    pcd,
    depth = depth_int
  )
  
  # Extract the reconstructed mesh
  poisson_mesh <- result[[1]]
  
  # Convert the reconstructed mesh back to R's mesh3d object
  numpy <- import("numpy")
  vertices_np <- numpy$asarray(poisson_mesh$vertices)
  vertices_recon <- py_to_r(vertices_np)
  
  # Convert Open3D Vector3iVector to a NumPy array
  triangles_np <- numpy$asarray(poisson_mesh$triangles)
  triangles_recon <- py_to_r(triangles_np) # Convert NumPy array to R matrix
  triangles_recon <- triangles_recon + 1 # Add 1 to convert from 0-based to 1-based indexing
  
  # Create mesh3d object in R
  mesh3d_recon <- tmesh3d(
    vertices = t(vertices_recon),
    indices = t(triangles_recon),
    homogeneous = FALSE
  )
  
  # Visualize the original and reconstructed meshes
  # open3d()
  # shade3d(mesh3d_recon, color = "magenta", alpha = 0.7)
  
  return(mesh3d_recon)
}



