#' An empirical chi matrix
#'
#' @description
#' The function `ChiMatrixMC` takes a inverted-Pareto data matrix and returns an empirical chi matrix regardless of rank transformations.
#' 
#' @param Dat An \eqn{n\times d} data matrix from the inverted Pareto distribution. 
#' @param quan Numeric; a lower threshold quantile
#' 
#' @return A \eqn{d\times d} strict upper triangular matrix of empirical chi matrix.
#' @export
#'
#' @examples 
#' StrMtx <- matrix(c(1, 1, 2, 2, 4,
#' 0, 2, 1, 3, 2,
#' 0, 0, 3, 1, 3,
#' 0, 0, 0, 4, 1,
#' 0, 0, 0, 0, 5),5,byrow = TRUE)
#' ParMtx <- matrix(c(0, 1.5, 2, 2.5, 2,
#'                 0, 0, 2, 2.5, 0.7,
#'                 0, 0, 0, 0.4, -0.3,
#'                 0, 0, 0, 0, 0.1,
#'                 0, 0, 0, 0, 0),5,byrow = TRUE)
#' FamMtx <- matrix(c(0, 1, 2, 3, 4,
#'                    0, 0, 3, 4, 1,
#'                    0, 0, 0, 3, 1,
#'                    0, 0, 0, 0, 1,
#'                    0, 0, 0, 0, 0),5,byrow = TRUE)
#' # X-vine specification
#' XVS=XVineSpec(M = StrMtx, Mmod = FamMtx, Mpar = ParMtx)
#' # Pareto random samples
#' Dat_P=ParetoSim(n = 5000, XVS = XVS) # Pareto scale
#' ChiMatrixMC(1/Dat_P)
ChiMatrixMC <- function(Dat, quan=NULL)
{
  if (!is.matrix(Dat)) {
    stop("The data should be a matrix")
  }
  if (ncol(Dat) <= 1) {
    stop("The data should be a matrix with at least two columns.")
  }
  if (!is.null(quan)) {
    Dat <- ParetoTransRank(data = Dat, u_quan = quan, scaleType = "U")
    #data <- data2mpareto(data, p)
  }
  d <- ncol(Dat)
  ChiMatrixMC <- matrix(0, d, d)
  
  for(i in 1:d){ # i=1, j=1,2,3,4,5 / i=2, j=2,3,4,5 ...
    for(j in i:d){ 
      Z=Dat[,c(i,j)] # 12
      Z1.2=Z[Z[,2]<1,] #(X1,X2) given X2<1
      Z2.1=Z[Z[,1]<1,c(2,1)] #(X2,X1) given X1<1
      
      chi1 <- sum(Z1.2[,1] < 1)/nrow(Z1.2) # X1<1 | X2<1
      chi2 <- sum(Z2.1[,1] < 1)/nrow(Z2.1) # X2<1 | X1<1
      chi <- mean(c(chi1,chi2))
      ChiMatrixMC[i,j]=chi 
    }
  }
  ChiMatrixMC=ChiMatrixMC+t(ChiMatrixMC)-diag(diag(ChiMatrixMC))
  return(ChiMatrixMC)
}
