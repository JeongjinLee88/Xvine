#' Specify the selection criterion for regular vine trees
#' 
#' @description
#' `TreeCrit()` chooses the selection criterion for regular vine trees. Possible criterion for the first tree
#' are the lower tail dependence measure called `chi` and the `variogram` for the Husler-Reiss distribution.
#' Other criteria for subsequent trees are Kendall's tau, Spearsman rho, information critera such as AIC, BIC, and adjusted AIC. 
#' This function is built upon the function [set_treecrit()] in `VineCopula` package.
#' Users can add additional tree criteria to the function.
#' 
#' 
#' @param treecrit A function that includes three arguments (u1,u2,weights). The first two arguments are for bivariate parametric models and the last argument 'weights' can be 
#' used for weighting missing values.
#' @param famset A numeric vector; the list of families for pair-copula models.
#'
#' @return The defined function `treecrit` used to calculate edge weights.
#' @export
#'
#' @examples
#' TreeCrit('chi')
TreeCrit <- function(treecrit, famset){
  if (is.function(treecrit)) {
    w <- try(treecrit(u1 = runif(10), u2 = runif(10), weights = rep(1, 
                                                                    10)), silent = TRUE)
    if (inherits(w, "try-error")) 
      stop("treecrit must be of the form 'function(u1, u2, weights)'")
    if (!is.numeric(w) || length(w) > 1) 
      stop("treecrit does not return a numeric scalar")
  }
  else if (all(treecrit == "chi")) {
    treecrit <- function(x1, x2, weights) {
      complete.i <- which(!is.na(x1 + x2))
      if (length(complete.i) < 10) {
        chi <- 0
      }
      else {
        complete.freq <- mean(!is.na(x1 + x2))
        Z=cbind(x1,x2)
        Z1.2=Z[Z[,2]<1,] #(X1,X2) given X2<1
        Z2.1=Z[Z[,1]<1,c(2,1)] #(X2,X1) given X1<1
        #chi <- sum(Z1.2[,1] < 1)/nrow(Z1.2)  # under the rank trans, chi1=chi2
        chi1 <- sum(Z1.2[,1] < 1)/nrow(Z1.2)
        chi2 <- sum(Z2.1[,1] < 1)/nrow(Z2.1)
        chi <- mean(c(chi1,chi2))
        chi * sqrt(complete.freq)
      }
    }
  }
  else if (all(treecrit == "variogram")) {
    treecrit <- function(x1, x2, weights) {
      Z=cbind(x1,x2)
      -emp_vario(1/Z)[1,2]
    }
  }
  else if (all(treecrit == "tau")) {
    treecrit <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 10) {
        tau <- 0
      }
      else {
        complete.freq <- mean(!is.na(u1 + u2))
        tau <- abs(VineCopula:::fasttau(u1[complete.i], u2[complete.i], 
                           weights))
        tau * sqrt(complete.freq)
      }
    }
  }
  else if (all(treecrit == "rho")) {
    treecrit <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 10) {
        tau <- 0
      }
      else {
        complete.freq <- mean(!is.na(u1 + u2))
        rho <- abs(cor(u1, u2, method = "spearman", use = "complete.obs"))
        rho * sqrt(complete.freq)
      }
    }
  }
  else if (all(treecrit == "AIC")) {
    treecrit <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 2) {
        AIC <- 0
      }
      else {
        AIC <- -suppressWarnings(BiCopSelect(u1[complete.i], 
                                             u2[complete.i], familyset = famset, weights = weights)$AIC)
      }
      AIC
    }
  }
  else if (all(treecrit == "BIC")) {
    treecrit <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 2) {
        BIC <- 0
      }
      else {
        BIC <- -suppressWarnings(BiCopSelect(u1[complete.i], 
                                             u2[complete.i], familyset = famset, weights = weights)$BIC)
      }
      BIC
    }
  }
  else if (all(treecrit == "cAIC")) {
    treecrit <- function(u1, u2, weights) {
      complete.i <- which(!is.na(u1 + u2))
      if (length(complete.i) < 2) {
        cAIC <- 0
      }
      else {
        fit <- suppressWarnings(BiCopSelect(u1[complete.i], 
                                            u2[complete.i], familyset = famset, weights = weights))
        n <- fit$nobs
        p <- fit$npars
        cAIC <- -(fit$AIC + 2 * p * (p + 1)/(n - p - 
                                               1))
      }
      cAIC
    }
  }
  else {
    txt1 <- "treecrit must be one of \"tau\", \"rho\", \"AIC\", \"BIC\", \"cAIC\""
    txt2 <- "or a function like function(u1, u2, weights) ... returning a numeric value."
    stop(paste(txt1, txt2))
  }
  treecrit
}
