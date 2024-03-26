#' Reproduce the given X-vine matrix specification depending on each conditioning variable
#'
#' @description
#' Given the initial X-vine matrix specification of (structure matrix, family matrix, parameter matrix),
#' `XVineSpec()` generates a list of redefined X-vine specifications for each conditioning variable.
#' Preserving the vine structure, the function reconstructs structure matrices where the first diagonal element
#' , denoted as \eqn{m_{11}\in [d]}, corresponds to the conditioning variable.
#' Additionally, `XVineSpec()` includes the associated parameter and family matrices based on the given structure matrices.
#' 
#' @details
#' To use the exact simulation algorithm for multivariate Pareto distribution in Engelke and Hitz (2020) and Kiriliouk et al. (2023),
#' we condition on each variable to calculate extremal functions.
#' For each conditioning variable, the index of the conditioning variable is placed on the first diagonal element in the structure matrix.
#' Specifically, in a \eqn{d \times d} structure matrix, the element \eqn{m_{jj}\in [d]} corresponds to the index of the conditioning variable \eqn{X_j} for \eqn{j=1,\ldots,d}.
#' All reproduced structure matrices maintain the same vine tree structure.
#' The corresponding parameter matrices contain the parameter values for the bivariate (tail) copulas.
#' In family matrices, the first row specifies bivariate tail copulas and the remaining rows specify bivariate copulas.
#' 
#' 
#' @param M A \eqn{d \times d} upper triangular structure matrix.
#' @param Mmod A \eqn{d \times d} strict upper triangular family matrix for identifying
#'  bivariate (tail) copulas.
#' @param Mpar A \eqn{d \times d} strict upper triangular parameter matrix associated with the family matrix.
#'
#' @return 
#' A list of arrays includes:
#' * xmat: array of dimension \eqn{d \times d \times d} containing rearranged X-vine structure matrices.
#' * fmat: array of dimension \eqn{d \times d \times d} containing family matrices associated with the structure matrices.
#' * pmat: array of dimension \eqn{d \times d \times d} containing parameter matrices associated with the structure matrices.
#' * xedge: A list of edges corresponding to the initial structure matrix.
#' @export
#' @references Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
#' 
#' Engelke, S., & Hitz, A. S. (2020). Graphical models for extremes. Journal of the Royal Statistical Society Series B: Statistical Methodology, 82(4), 871-932.
#' 
#'
#' @seealso \code{VineEdges()}
#' @examples
#' ##  A 6 x 6 structure matrix
#' M <- cbind(c(1,0,0,0,0,0), c(1,4,0,0,0,0), c(4,1,6,0,0,0),
#'           c(4,1,6,3,0,0), c(6,4,1,3,5,0), c(3,4,1,6,5,2))
#' ##  A 6 x 6 family matrix
#' Mmod <- cbind(rep(0,6),c("neglog", rep(0,5)),
#'               c("hr", "clayton",  rep(0,4)),
#'               c("log", "gumbel", "joe", rep(0,3)),
#'               c("dir", "frank", "clayton", "gumbel", 0, 0),
#'               c("neglog", "joe", "frank", "clayton", "gumbel", 0))
#' ##  A 6 x 6 parameter matrix
#' Mpar <- cbind(rep(0,6),c(0.25, rep(0,5)), c(0.3, 0.7,  rep(0,4)),
#'               c(0.5,0.7,0.1, rep(0,3)), c(0.2,0.8,0.3,0.55, 0, 0),
#'               c(0.1,0.9,0.75,0.3,0.45, 0))
#' ##  Permute the given X-vine specification
#' PermXVineSpec <- XVineSpec(M,Mmod,Mpar)
#' PermXVineSpec[[1]]
#' cbind(PermXVineSpec[[1]][,,5],PermXVineSpec[[2]][,,5])
#' ##  Compare with cbind(M,Mmod): looks ok!
#' ##  Truncated X-vine specification at tree level 3
#' M <- cbind(c(1,0,0,0,0,0), c(1,4,0,0,0,0), c(4,1,6,0,0,0),
#'            c(4,1,6,3,0,0), c(6,4,1,0,5,0), c(3,4,1,0,0,2))
#' Mmod <- cbind(rep(0,6),c("neglog", rep(0,5)), c("hr", "clayton",  rep(0,4)), 
#'            c("log", "gumbel", "joe", rep(0,3)), c("dir", "frank", "clayton", 0, 0, 0),
#'            c("neglog", "joe", "frank", 0,0, 0))
#' Mpar <- cbind(rep(0,6),c(0.25, rep(0,5)), c(0.3, 0.7,  rep(0,4)), 
#'               c(0.5,0.7,0.1, rep(0,3)), c(0.2,0.8,0.3,0, 0, 0),
#'               c(0.1,0.9,0.75,0,0, 0))
#' PermTrucXVine <- XVineSpec(M,Mmod,Mpar)
#' PermTrucXVine[[1]]
#' cbind(PermTrucXVine[[1]][,,5],PermTrucXVine[[2]][,,5]) #compare with cbind(M,Mmod): looks ok!
XVineSpec <- function(M, Mmod, Mpar){
  d <- nrow(M)
  vedges <- VineEdges(M,Mmod,Mpar)
  
  kinit <- M[1,1]
  Dtotal <- c(1:d)[-kinit] 
  mlist <- mlistmod <- mlistpar <- vector('list', length = d)
  mlist[[kinit]] <- M
  mlistmod[[kinit]] <- Mmod
  mlistpar[[kinit]] <- Mpar
  for(k in Dtotal){
    matr <- matrmod <- matrpar <- matrix(0,nrow=d,ncol=d)
    edgtemp <- vedges$edges
    l <- length(edgtemp)
    matr[1,1] <- k
    for(i in d:2){
      inds <- c(edgtemp[[l]][1:2,])
      temp <- unique(as.numeric(names(which(table(inds) == 1)))) 
      if(i > 2 & l >= 2){
        start <- max(temp[which(!(temp %in% c(k,unique(c(edgtemp[[l]][3:(l+1),])))))])
      } else{
        start <- max(temp[which(temp != k)])
      }
      ind <- sapply(edgtemp, function(mat) which(mat[1:2,,drop=F] == start, arr.ind = T)[1,2])
      temp <- unlist(sapply(c(1:l), function(j) edgtemp[[j]][1:2,ind[j]]))
      diagn <- temp[temp!=start]
      if(length(diagn) != (i-1)){
        diagn <- c(diagn,rep(0,(i-1)-length(diagn)))
      }
      matr[1:(i-1),i] <- diagn
      edgtemp <- lapply(edgtemp, function(mat) { 
        mat <- mat[,apply(mat, 2, function(col) all(col[1:2] != start)), drop = F]
      })
      matr[i,i] <- start 
      if(dim(edgtemp[[length(edgtemp)]])[2] == 0){
        edgtemp <- edgtemp[-length(edgtemp)]
        l <- l-1
      }
    }
    
    edgnew <- VineEdges(matr,d = d)$edges 
    for(i in 1:length(edgnew)){
      ind <- apply(edgnew[[i]], 2, function(j) which(apply(vedges$edges[[i]], 2, function(x) all(x== j))))
      matrpar[i,(i+1):d] <- vedges$pars[[i]][ind]
      matrmod[i,(i+1):d] <- vedges$models[[i]][ind]
    }
    mlist[[k]] <- matr
    mlistmod[[k]] <- matrmod
    mlistpar[[k]] <- matrpar
  }
  
  xmat=array(unlist(mlist),dim = c(d,d,d)) # structure mtx
  fmat=array(unlist(mlistmod),dim = c(d,d,d)) # family mtx
  pmat=array(unlist(mlistpar),dim = c(d,d,d)) # parameter mtx
  
  #return(list('matrices' = mlist, 'matricesPar' = mlistpar, 'matricesMod' = mlistmod, 'vedges'=vedges$edges))
  return(list('xmat' = xmat, 'fmat' = fmat, 'pmat' = pmat, 'xedge'=vedges$edges))
}





