#' Store the positions of variable indices in a structure matrix
#' 
#' @description
#' `storeVarInd()` stores the positions of variable indices (\eqn{1,\ldots,d})
#' in a structure matrix.
#' The output vector will be used to rearrange the order of generated variables
#'  from [XVineSim()], which may not necessarily align with \eqn{1,\ldots,d}.
#' 
#' @param Diag Numeric \eqn{d \times d} structure matrix.
#'
#' @return A numeric vector containing the locations of variable indices.
#' @export
#'
#' @examples
#' ##  
storeVarInd <- function(Diag){
  d=length(Diag)
  order.new=rep(NA,d)
  for(i in 1:d){
    order.new[i] <- which(Diag==i)
  }
  return(order.new)
}
