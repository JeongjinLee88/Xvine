# If format = FALSE, then 'edges' is of the same form as the result of VineEdges. 
# So, it is a list of length (d-1) of matrices, where the first two rows of each matrix represent conditioned variables,
# and the following rows representing conditioning variables.
# If format = TRUE, then 'edges' is a list of length (d-1), where each list element is a list of two matrices, one representing the
# conditioned variables (transposed), the other representing the conditioning variables (transposed).

# Of course, the structure matrix is not unique. When choosing between two indices, we always pick the largest one
# Is this a good convention?
XVineEdges <- function(M, Mmod = NULL, Mpar = NULL, d = ncol(M), format = FALSE){ 
  if((is.null(Mmod) & !is.null(Mpar)) | (!is.null(Mmod) & is.null(Mpar))){
    return("If you specify Mpar, you need to specify Mmod, and vice versa")
  } else{
    lmax <- d
    test <- which(M == 0, arr.ind = T)
    indics <- test[test[,2]>test[,1],1]
    if(length(indics) != 0){lmax <- min(indics)} 
    edges <- models <- pars <- vector('list', length = lmax-1)
    edges[[1]] <- sapply(c(2:d), function(j) sort(M[c(1,j),j]))
    if(lmax > 2){
      for(i in 2:(lmax-1)){
        edges[[i]] <- sapply(c((i+1):d), function(j) c(sort(M[c(i,j),j]),sort(M[(i-1):1,j])))
      }
    }
    if(!is.null(Mmod)){
      models[[1]] <- Mmod[1,2:d]
      pars[[1]] <- Mpar[1,2:d]
      if(lmax > 2){
        for(i in 2:(lmax-1)){
          models[[i]] <- Mmod[i,(i+1):d]
          pars[[i]] <- Mpar[i,(i+1):d]
        }
      }
    }
    if(format){
      edgestemp <- vector('list', length = lmax-1)
      edgestemp[[1]] <- list('C' = t(edges[[1]]), 'D' = NULL)
      if(lmax > 2){
        for(i in 2:(lmax-1)){
          edgestemp[[i]] <- list('C' = t(edges[[i]][1:2,,drop=F]), 'D' = t(edges[[i]][3:(i+1),,drop=F]))
        }
      }
      edges <- edgestemp
    }
    return(list('edges' = edges, 'models' = models, 'pars' = pars))
  }
}