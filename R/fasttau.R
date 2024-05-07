fastTau <- function (x, y, weights = NA) 
{
  if (any(is.na(weights))) {
    m <- length(x)
    n <- length(y)
    if (m == 0 || n == 0) 
      stop("both 'x' and 'y' must be non-empty")
    if (m != n) 
      stop("'x' and 'y' must have the same length.")
    out <- .C("ktau", x = as.double(x), y = as.double(y), 
              N = as.integer(n), tau = as.double(0), S = as.double(0), 
              D = as.double(0), T = as.integer(0), U = as.integer(0), 
              V = as.integer(0), PACKAGE = "VineCopula")
    ktau <- out$tau
  }
  else {
    ktau <- TauMatrix(matrix(c(x, y), length(x), 2), weights)[2, 
                                                              1]
  }
  return(ktau)
}