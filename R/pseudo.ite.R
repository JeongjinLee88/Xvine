pseudo.ite <- function(condSet, treeMinus, data, MST, VineTree){
  Vine_temp=lapply(1:treeMinus,function(i)NULL)
  #Vine_temp[[1]] all bivariate exponent densities in T1
  #Vine_temp[[2]] all pair-copulas in T2 (single conditioning variable)
  #Vine_temp[[3]] all pair-copulas in T3 (two conditioning variables)
  for(j in 2:treeMinus){
    if(j==2){
      d=nrow(MST[[j]]$E$nums)
      for(i in 1:d){ #i=1
        con <- MST[[j]]$E$nums[i, ]  #MST[[2]]
        PairEdge <- MST[[1]]$E$nums[con, ] #VineTree[[1]]
        common=ifelse(test = any(PairEdge[1,]==PairEdge[2,1]), PairEdge[2,1], PairEdge[2,2])
        Output=Condfit.FirstTree(condSet = condSet,data = data,VineTree = VineTree[[1]])  
        Vine_temp[[1]]$CondData.1=Output$Pseudo1
        Vine_temp[[1]]$CondData.2=Output$Pseudo2
        if (PairEdge[1, 1] == common) {
          zr1 <- Vine_temp[[1]]$CondData.2[[con[1]]]
          n1 <- VineTree[[1]]$E$Copula.CondName.2[con[1]]
        }else {
          zr1 <- Vine_temp[[1]]$CondData.1[[con[1]]]
          n1 <- VineTree[[1]]$E$Copula.CondName.1[con[1]]
        }
        if (PairEdge[2, 1] == common) {
          zr2 <- Vine_temp[[1]]$CondData.2[[con[2]]]
          n2 <- VineTree[[1]]$E$Copula.CondName.2[con[2]]
        }else {
          zr2 <- Vine_temp[[1]]$CondData.1[[con[2]]]
          n2 <- VineTree[[1]]$E$Copula.CondName.1[con[2]]
        }
        CoPam  <- do.call(rbind,VineTree[[j]]$E$Copula.param)[,1] # from VineTree[[2]]
        CoFam  <- VineTree[[j]]$E$Copula.type
        #CondOn.1 <- CondCop(u1 = zr1,u2 = zr2,par = CoPam[i],family = CoFam[i])
        #CondOn.2 <- CondCop(u1 = zr2,u2 = zr1,par = CoPam[i],family = CoFam[i])
        if(!length(zr1)==length(zr2)){
          CondOn.1 <- 0
          CondOn.2 <- 0
        }
        if(length(zr1)==length(zr2)){
          CondOn.1 <- BiCopHfunc2(zr1, zr2, family = CoFam[i], par = CoPam[i], par2 = 0, check.pars = FALSE)
          CondOn.2 <- BiCopHfunc2(zr2, zr1, family = CoFam[i], par = CoPam[i], par2 = 0, check.pars = FALSE)
        }
        Vine_temp[[j]]$CondData.1[i]<- list(CondOn.1)
        Vine_temp[[j]]$CondData.2[i]<- list(CondOn.2)
      }
    }
    if(j > 2){ # should separate j=2 from j > 2
      d=nrow(MST[[j]]$E$nums)
      for(i in 1:d){ #i=12
        con <- MST[[j]]$E$nums[i, ]  #MST[[3]]
        PairEdge <- MST[[j-1]]$E$nums[con, ] #VineTree[[2]]
        common=ifelse(test = any(PairEdge[1,]==PairEdge[2,1]), PairEdge[2,1], PairEdge[2,2])
        if (PairEdge[1, 1] == common) {
          zr1 <- Vine_temp[[j-1]]$CondData.2[[con[1]]]
          n1 <- VineTree[[j-1]]$E$Copula.CondName.2[con[1]]
        }else {
          zr1 <- Vine_temp[[j-1]]$CondData.1[[con[1]]]
          n1 <- VineTree[[j-1]]$E$Copula.CondName.1[con[1]]
        }
        if (PairEdge[2, 1] == common) {
          zr2 <- Vine_temp[[j-1]]$CondData.2[[con[2]]]
          n2 <- VineTree[[j-1]]$E$Copula.CondName.2[con[2]]
        }else {
          zr2 <- Vine_temp[[j-1]]$CondData.1[[con[2]]]
          n2 <- VineTree[[j-1]]$E$Copula.CondName.1[con[2]]
        }
        CoPam  <- do.call(rbind,VineTree[[j]]$E$Copula.param)[,1] # from VineTree[[2]]
        CoFam  <- VineTree[[j]]$E$Copula.type
        if(!length(zr1)==length(zr2)){
          CondOn.1 <- 0
          CondOn.2 <- 0
        }
        if(length(zr1)==length(zr2)){
          CondOn.1 <- BiCopHfunc2(zr1, zr2, family = CoFam[i], par = CoPam[i], par2 = 0, check.pars = FALSE)
          CondOn.2 <- BiCopHfunc2(zr2, zr1, family = CoFam[i], par = CoPam[i], par2 = 0, check.pars = FALSE)
        }
        Vine_temp[[j]]$CondData.1[i]<- list(CondOn.1)
        Vine_temp[[j]]$CondData.2[i]<- list(CondOn.2)  
      }
    }
  }
  Vine_temp        
}