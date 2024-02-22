#' Permute the diagonal elements of a structure matrix
#'
#' @description
#' `ReorderStrMtx()` reorders the indices of variables to ensure that the
#'  diagonal elements of a structure matrix are arranged in ascending order.
#' 
#' @param Matrix Numeric \eqn{d \times d} structure matrix
#' @param oldOrder Logical; If `NULL` (the default), the function uses the order
#' of diagonal elements in the structure matrix. 
#' 
#' @return Numeric \eqn{d \times d} structure matrix with diagonal elements
#' arranged in ascending order.
#' @export
#'
#' @examples
#' ##  A 5 x 5 structure matrix where diagonal elements are not arranged in ascending order.
#' StrMtx <- cbind(c(2,0,0,0,0),c(2,4,0,0,0),c(4,2,5,0,0),c(2,4,5,3,0),c(2,3,4,5,1))
#' 
#' ##  A 5 x 5 structure matrix with diagonal elements arranged in increasing order.
#' ReorderStrMtx(StrMtx)
ReorderStrMtx <- function (Matrix, oldOrder = NULL)
{
  if (length(oldOrder) == 0) {
    oldOrder <- diag(Matrix)
  }
  O <- apply(t(1:nrow(Matrix)), 2, "==", Matrix)
  for (i in 1:nrow(Matrix)) {
    Matrix[O[, oldOrder[i]]] <- i
  }
  return(Matrix)
}