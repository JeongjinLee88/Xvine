#' Redefine pseudo-observations for the conditional tail copula above the second level tree
#'
#' @description
#' `Condfit.FirstTree()` redefine pseudo-observations for the conditional tail copula above the second level tree
#'  where the cardinality of the conditioning set is greater than 2.
#' 
#' @param condSet A numeric vector of the conditioning set where the cardinality is greater than 2.
#' @param data A \eqn{n\times d} inverted-Pareto sample data set.
#' @param VineTree A list of selected vine structure and fitted vine models.
#'
#' @return A list of pseudo-observations for the conditional tail copula:
#' * Pseudo1
#' * Pseudo2
#' @export
#'
Condfit.FirstTree <- function (condSet, data, VineTree) 
{
  Pseudo1=list()
  Pseudo2=list()
  d <- nrow(VineTree$E$nums)
  Pam  <- do.call(rbind,VineTree$E$Copula.param)[,1]
  Fam  <- VineTree$E$Copula.type
  #con <- MST$E$nums[1, ]
  #temp <- oldVineGraph$E$nums[con[1], ]
  #temp <- oldVineGraph$E$nums[con[2], ]
  #condSet <- temp
  #condSet=c(7,8,9);condSet=c(3,4,6)
  CondData <- data[apply(data[,c(condSet)],1,function(x)all(x<1)),]
  #CondData=data[,all(data[,1:29]<1)] # when no data points in CondData
  for (i in 1:d) { # i=10
    ind <- VineTree$E$nums[i, ]
    if(length(CondData)==0){
      Pseudo1[i] <- list(0)
      Pseudo2[i] <- list(0)
    }else{
      if(is.vector(CondData)){ # if CondData is a vector
        Z=CondData[c(ind)]
        if(Z[2]<=1){
          Z1.2=Z #lambda1.2 given X2 < 1
          CondOn.1 <- CondTC(x1 = Z1.2[1],x2 = Z1.2[2],par = Pam[i],family = Fam[i])
        }else{
          CondOn.1 <- 0
        }
        if(Z[1]<=1){
          Z2.1=rev(Z)
          CondOn.2 <- CondTC(x1 = Z2.1[1],x2 = Z2.1[2],par = Pam[i],family = Fam[i])
        }else{ 
          CondOn.2 <- 0
        }
      }else{ # if CondData is a matrix
        Z=CondData[,c(ind)]
        Z1.2=Z[Z[,2]<=1,] #lambda1.2 given X2 < 1
        if(length(Z1.2)==0){
          CondOn.1 <- 0
        }
        if(is.vector(Z1.2)){
          CondOn.1 <- CondTC(x1 = Z1.2[1],x2 = Z1.2[2],par = Pam[i],family = Fam[i])
        }
        if(is.matrix(Z1.2)){
          CondOn.1 <- CondTC(x1 = Z1.2[,1],x2 = Z1.2[,2],par = Pam[i],family = Fam[i])
        }
        Z2.1=Z[Z[,1]<=1,c(2,1)] #lambda2.1 given X1 < 1
        if(length(Z2.1)==0){
          CondOn.2 <- 0
        }
        if(is.vector(Z2.1)){
          CondOn.2 <- CondTC(x1 = Z2.1[1],x2 = Z2.1[2],par = Pam[i],family = Fam[i])
        }
        if(is.matrix(Z2.1)){
          CondOn.2 <- CondTC(x1 = Z2.1[,1],x2 = Z2.1[,2],par = Pam[i],family = Fam[i])
        }
      }
      #min(nrow(Z1.2),nrow(Z2.1)) < 10 #truncate the tree
      Pseudo1[i] <- list(CondOn.1)
      Pseudo2[i] <- list(CondOn.2)
    }
  }
  return(list("Pseudo1"=Pseudo1,"Pseudo2"=Pseudo2))
}



