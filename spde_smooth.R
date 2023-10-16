smooth.construct.spde.test <- function(data,coords=NULL,mesh_in, knots){
  if(any(class(data) == "sf") & is.null(coords)){
    x <- st_coordinates(data)
  } else{
    if(!is.null(coords)){
      # observation locations
      x <- matrix(0, nr = length(data[[1]]), nc = 2) 
      x[,1] <- data[[coords[[1]]]]
      x[,2] <- data[[coords[[2]]]]
    } else{rlang::abort("Need to provide sf object or coords column names")}
  }
  # setup accept user mesh
  if (!any(class(mesh_in)=="inla.mesh" ) ) stop("xt must be NULL or an inla.mesh object")
  mesh <- mesh_in
  
  # browser()
  # model matrix: projects parameters to observation locations on mesh 
  X <- as.matrix(fm_basis(mesh, x))
  # compute finite element matrices used as smoothing penalty matrices 
  inlamats <- fm_fem(mesh)
  S <- list()
  S[[1]] <- as.matrix(inlamats$c1)
  S[[2]] <- 2 * as.matrix(inlamats$g1)
  S[[3]] <- as.matrix(inlamats$g2)
  # L is a matrix with a column for each smoothing parameter (tau, kappa) 
  # and a row for each smoothing matrix (c1, g1, g2). 
  # The (i,j)^th entry of L contains the power that smoothing parameter i 
  # is computed to before being multiplied by smoothing matrix j. 
  # E.g. If (1, 2) has value 4, then smoothing parameter 2 (kappa) is taken
  # to the power 4 before being multiplied by smoothing matrix 1 (c1): i.e. kappa^4*c1
  # All of these computations for each element of L are then summed to create a single
  # smoothing matrix. 
  L <- matrix(c(2,2,2,4,2,0), ncol = 2)
  # Rank is the basis dimension, it is repeated for each smoothing matrix 
  rank <- rep(knots,3)
  # As kappa > 0, the null space of the Matern SPDE is empty 
  null.space.dim <- 0 
  
  df <- ncol(X)     # maximum DoF (if unconstrained)
  # Give object a class
  object <- list(A = X, S = S, L = L,
                 rank = rank, 
                 null.space.dim = null.space.dim,
                 df = df)
  class(object) <- "spde.smooth" 
  return(object)
}