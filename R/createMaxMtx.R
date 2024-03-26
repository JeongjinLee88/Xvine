#' Create a maximum matrix from a structure matrix with diagonal elements in ascending order
#'
#' @description
#' `createMaxMtx()` creates a max-matrix from a structure matrix where diagonal elements are put in increasing order.
#' This max-matrix is used in sequential parameter estimation
#'  to determine the appropriate argument for conditional pair-copulas.
#' For more details, refer to Chapter 6 in Czado (2019).
#' @param Matrix Numeric \eqn{d \times d} upper triangular matrix that specifies
#' the regular vine structure.
#'
#' @return Numeric \eqn{d \times d} upper triangular max-matrix.
#' @export
#' @references Czado, C. (2019). Analyzing dependent data with vine copulas.
#'  Lecture Notes in Statistics, Springer, 222.
#' @examples
#' ##  Create a 5 x 5 upper triangular structure matrix
#' StrMtx <- matrix(c(1, 1, 2, 2, 4,
#' 0, 2, 1, 3, 2,
#' 0, 0, 3, 1, 3,
#' 0, 0, 0, 4, 1,
#' 0, 0, 0, 0, 5),5,byrow = TRUE)
#' ##  Derive the corresponding 5 x 5 upper triangular max-matrix
#' MaxMtx <- createMaxMtx(StrMtx)
createMaxMtx <- function(Matrix) 
{
  if (dim(Matrix)[1] != dim(Matrix)[2]) 
    stop("Structure matrix has to be quadratic.")
  MaxMat <- Matrix
  n <- nrow(MaxMat)
  for (j in n:2) {
    # j=5, 4, 3, 2
    # i=5,4,3,2 / 4,3,2 / 3,2 / 2  
    for (i in j:2) {
      MaxMat[i, j] <- max(MaxMat[i:1, j])
    }
  }
  return(MaxMat)
}




