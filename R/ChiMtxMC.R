#' Pairwise tail dependence matrix via Monte Carlo simulations
#'
#' @description
#' `ChiMtxMC` takes an inverted-Pareto data matrix as input and computes an empirical pairwise tail dependence measure via Monte Carlo simulations.
#' The resulting matrix is upper strict and computes the tail dependence coefficient between variables numerically.
#'  Each element in the chi-matrix consists of the (lower) tail dependence coefficient \eqn{\chi_{L,e}=R(1,1)} for \eqn{e\in T_i}, \eqn{i=1,\ldots,d-1}.
#' 
#' @param Dat An \eqn{n\times d} data matrix from multivariate inverted-Pareto distribution.
#' If the original data are used, then the data transformation is performed to obtain pseudo-observations for the inverted-Pareto distribution
#'  See \code{ParetoTransRank()}.
#' @param Quan A numeric quantile, set a lower enough quantile \eqn{u_{quan}\in(0,1)} to convert data into pseudo-observations for the inverted-Pareto distribution.
#' Default is `Quan=NULL`.
#' 
#' @return A \eqn{d\times d} strict upper triangular matrix of empirical pairwise chi-matrix.
#' @export
#' 
#' @seealso [ParetoTransRank()]
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
#' Dat_P=ParetoSim(n = 2000, XVS = XVS) # Pareto scale
#' # Pairwise chi-matrix
#' ChiMtxMC(Dat_P)
ChiMtxMC <- function(Dat, Quan=NULL)
{
  if (!is.matrix(Dat)) {
    stop("The data should be a matrix")
  }
  if (ncol(Dat) <= 1) {
    stop("The data should be a matrix with at least two columns.")
  }
  if (!is.null(Quan)) {
    Dat <- ParetoTransRank(data = Dat, u_quan = Quan, scaleType = "U")
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
