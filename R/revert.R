#' Convert a lower (or upper) vine structure matrix
#' 
#' @description
#' `revert()` transforms a lower (or upper) triangular vine structure matrix
#'  to the corresponding upper (or lower) triangular matrix.
#' 
#' @param m A numeric \eqn{d\times d} lower (or upper) triangular matrix specifying
#' the regular vine structure.
#'
#' @return A numeric \eqn{d\times d} upper triangular matrix if the input matrix is lower
#' triangular; otherwise, it is a lower triangular matrix.
#' @export
#'
#' @examples 
#' ##  A 5 x 5 lower triangular structure matrix
#' LowerStrMtx <- matrix(c(5, 2, 3, 1, 4,
#'  0, 2, 3, 4, 1,
#'  0, 0, 3, 4, 1,
#'  0, 0, 0, 4, 1,
#'  0, 0, 0, 0, 1),5,5)
#'  
#' ##  Convert to the upper triangular matrix
#' revert(LowerStrMtx)
revert <- function (m) 
{
  if (length(dim(m)) == 2) {
    return(m[nrow(m):1, ncol(m):1, drop = F])
  }
  else {
    return(m[nrow(m):1, ncol(m):1, , drop = F])
  }
}

