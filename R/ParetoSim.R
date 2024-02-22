#' Simulation from multivariate Pareto distribution
#'
#' @description
#' `ParetoSim()` generates samples from multivariate Pareto distribution by rejection sampling.
#' Drawing uniform samples from \eqn{k \in {1,\ldots,d}}, 
#' the \eqn{k}th conditioning index 'k' is passed to the function [XVineSim()]
#'  and then calculate its extremal function.
#'  For more details on the X-vine simulation algorithm, refer to Kiriliouk, A., Lee, J., & Segers, J. (2023).
#' For more details on the exact simulation algorithm, refer to Engelke, S., & Hitz, A. S. (2020).
#' @param n Integer; the number of d-dimensional observations to generate
#' @param XVS A list of three components: 
#' * permuted structure matrices.
#' * permuted family matrices.
#' * permuted parameter matrices.
#'
#' @return A \eqn{n\times d} data matrix of multivariate Pareto samples.
#' @export
#' @references Kiriliouk, A., Lee, J., & Segers, J. (2023). X-Vine Models for Multivariate Extremes. arXiv preprint arXiv:2312.15205.
#' 
#' Engelke, S., & Hitz, A. S. (2020). Graphical models for extremes. Journal of the Royal Statistical Society Series B: Statistical Methodology, 82(4), 871-932.
#' 
#' @examples
#' ##  A 5-dim X-vine model
#' StrMtx <- matrix(c(1, 1, 2, 2, 4,
#'                   0, 2, 1, 3, 2,
#'                   0, 0, 3, 1, 3,
#'                   0, 0, 0, 4, 1,
#'                   0, 0, 0, 0, 5),5,byrow = TRUE)
#'                   
#' ParMtx <- matrix(c(0, 1.5, 2, 2.5, 2,
#'                    0, 0, 2, 2.5, 0.7,
#'                    0, 0, 0, 0.4, -0.3,
#'                    0, 0, 0, 0, 0.1,
#'                    0, 0, 0, 0, 0),5,byrow = TRUE)
#'                    
#' FamMtx <- matrix(c(0, 1, 2, 3, 4,
#'                    0, 0, 3, 4, 1,
#'                    0, 0, 0, 3, 1,
#'                    0, 0, 0, 0, 1,
#'                    0, 0, 0, 0, 0),5,byrow = TRUE)
#'                    
#' ##  X-Vine specification
#' XVS=XVineSpec(M = StrMtx, Mmod = FamMtx, Mpar = ParMtx)
#' ##  Generate samples from multivariate Pareto distribution
#' Dat_P=ParetoSim(n = 5000, XVS = XVS)
ParetoSim <- function(n, XVS){
  
  d <- dim(XVS$xmat)[1]
  out <- numeric(0)
  ite <- 0
  
  while (ite < n) {
    shift <- sample(1:d, n, replace = TRUE)
    for (k in 1:d) {
      n.k <- sum(shift == k)
      if (n.k > 0) {
        Xk <- XVineSim(N = n.k, XVS=XVS, k=k)
        Uk <- Xk[,k]/Xk # extremal fts of Yk
        proc <- Uk/rowSums(Uk)/(1 - stats::runif(NROW(Uk)))
        idx.sim <- which(apply(proc, 1, max) > 1)
        out <- rbind(out, proc[idx.sim, ])
        ite <- NROW(out)
      }
    }
  }
  Psamp=out[sample(1:NROW(out), n, replace = FALSE), ]
  return(Psamp)
}

