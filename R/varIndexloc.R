#' Store the order of variable indices \eqn{(1,\ldots,d)} placed on the diagonals of a structure matrix
#' 
#' @description
#' `varIndexloc()` stores the order of variable indices \eqn{(1,\ldots,d)} placed on the diagonals of a \eqn{d \times d} structure matrix. 
#' The sequence of diagonal elements \eqn{(m_{11},\ldots,m_{dd})} corresponds to the one of the indices of variables.
#' For instance, if \eqn{(m_{11},m_{22},m_{33},m_{44},m_{55})=(5,4,2,3,1)}, then it indicates the order of variable indices \eqn{(5,3,4,2,1)}.
#' The stored order will be used to rearrange the data matrix generated from [XVineSim()] in ascending order.
#' 
#' 
#' @param Diag A numeric vector of diagonal elements of a structure matrix.
#'
#' @return A numeric vector returning the corresponding order of diagonal elements.
#' @export
#'
#' @examples
#' ##  A structure matrix
#' StrMtx=matrix(c(5,5,4,2,2,
#' 0,4,5,4,3,
#' 0,0,2,5,4,
#' 0,0,0,3,5,
#' 0,0,0,0,1),5,5,byrow=TRUE)
#' ##  Save diagonal elements
#' Delements=diag(StrMtx)
#' ##  Return the corresponding order of the diagonal elements
#' varIndexloc(Delements)
varIndexloc <- function(Diag){
  d=length(Diag)
  order.new=rep(NA,d)
  for(i in 1:d){
    order.new[i] <- which(Diag==i)
  }
  return(order.new)
}

