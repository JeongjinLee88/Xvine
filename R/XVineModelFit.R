#' Fit X-vine models
#' 
#' @description 
#' Fits X-Vine models to a d-dimensional inverted-Pareto data set. The function selects vine tree structures, using maximum spanning tree.
#' and select bivariate exponent measure classes and pair-copula classes and estimate parameters for the corresponding classes
#' Users can specify the type of edge weights for the first tree and subsequent trees separately and the class of bivariate exponent measures and pair-copulas
#' 
#' @param data a \eqn{n\times d} data matrix of multivariate Pareto samples.
#' @param tcfamset Numeric vector; the class of bivariate exponent measures.
#' @param pcfamset Numeric vector; the class of bivariate pair-copula models.
#' @param selectioncrit Character string; indicates the selection criteria for the class of bivariate exponent measures or pair-copulas.
#' @param trunclevel Numeric; indicates the specified truncation level.
#' @param progress Logical; whether the progress of selecting vine tree structures is printed.
#' @param treecritT1 a character string, indicating the tree criterion for the first tree.
#' @param treecritT2 a character string, indicating the tree criterion for subsequent trees.
#' @param weights Logical; whether weights should be assigned to observations when missing values exist.
#' @param cores Numeric; indicates the number of cores for parallel computing (optional).
#' @param Rank Logical; whether rank transformation is performed or not (\code{Rank=T}; default).
#' @param qt Numeric; a lower threshold for the rank transformation. It switches from Pareto scale to uniform scale.

#' @return a nested list object with maximum spanning trees and fitted vine trees.
#' @export
#'
#' @examples
XVineModelFit <- function(data, Rank=TRUE, Rank_chiU=TRUE, Rank_chiL=FALSE, Rank_chi3=TRUE, MST1_HR=FALSE, qt=0.2, N=8000
                           , XVS, tcfamset = c(1,2,3,4), pcfamset = c(0,1,3,4,5,6,13,14,16), selectioncrit = "AIC"
                           , BIC_graph=TRUE, Chi3_graph=FALSE, trunclevel = FALSE, progress = TRUE
                           , treecritT1 = "chi", treecritT2= 'tau', si=0.9, effsampsize=10, tau_threshold=0.05, se=FALSE, weights=NA, cores = 1)
{
  ##  Note that the 'XVineModelFit' function uses multivariate 'inverted' Pareto samples
  ##  If you directly use samples from the limiting distribution, they must be ones from Pareto distribution with Pareto margin.
  if(Rank){
    data <- ParetoTransRank(data = data, u_quan = qt, scaleType = "U") 
  }
  # Determine the truncation level if specified.
  d <- ifelse(test = trunclevel,trunclevel+1,ncol(data))
  d_col <- ncol(data)
  n <- nrow(data)
  # Determine the selection criteria for T_1 and T_i, i=2,...
  T1crit <- TreeCrit(treecrit = treecritT1)
  Ticrit <- TreeCrit(treecrit = treecritT2)
  #Ticrit <- VineCopula:::set_treecrit(treecrit = treecritT2)
  # Set NULL lists
  MST <- list()
  VineTree=list()
  # Fit the first MST (maximum spanning tree) and its extremal graph
  if(MST1_HR){
    MST[[1]] <- MST_HR(Dat_Pareto = data,quan = qt)
  }else{
    g <- VineCopula:::initializeFirstGraph(data, T1crit, weights)
    MST[[1]] <- VineCopula:::findMaxTree(g, mode = "RVine")
  }
  VineTree[[1]] <- fit.FirstTree(MST[[1]], data, tcfamset, si=si,
                                 selectioncrit = selectioncrit, cores = 1)
  if(d > 2){
    for (tree in 2:(d - 1)) { # tree=2,3,4
      g <- VineCopula:::buildNextGraph(VineTree[[tree-1]], weights, treecrit = Ticrit, 
                                       cores > 1, truncated = FALSE)
      MST[[tree]] <- VineCopula:::findMaxTree(g, mode = "RVine", truncated = FALSE)
      VineTree[[tree]] <- fit.SubTree(data = data, MST = MST, VineTree = VineTree, copfamset = pcfamset, tree = tree, si=si,
                                          selectioncrit, progress, effsampsize = effsampsize, tau_threshold = tau_threshold, weights = weights, 
                                          se = se, cores = cores)
    }  
  }
  
  mBIC_tree=do.call(c,lapply(1:(d-1),function(i)sum(VineTree[[i]]$E$mBIC)))
  mBIC_cum=cumsum(mBIC_tree[-1])
  TruncLevel_mBICmin=which.min(mBIC_cum)+1
  ThresholdLevel=which(sapply(1:(d-1),function(i)all(VineTree[[i]]$E$Copula.type==0)))
  if(identical(ThresholdLevel,integer(0))){
    TruncLevel_Ind=d_col-1
  }else{
    TruncLevel_Ind=min(which(sapply(1:(d-1),function(i)all(VineTree[[i]]$E$Copula.type==0))))
  }  
  TruncLevelStar=min(TruncLevel_mBICmin,TruncLevel_Ind)
  
  mBIC_g = ggplot() +
    geom_point(aes(x = 2:(d-1), y = mBIC_cum)) +
    geom_vline(xintercept = TruncLevelStar, linetype=3, size=1.5) +
    labs(title="",x=expression("Truncation Level (q)"), y = "mBIC values") +
    theme(
      axis.text.x = element_text(color = "black",face="bold", size=14),
      axis.text.y = element_text(face="bold", size = 14),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      panel.background = element_rect(fill = "white",
                                      colour = "white",
                                      linewidth = 0.5, linetype = "solid"),
      panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                      colour = "grey"), 
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                      colour = "grey"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    scale_x_continuous(limits = c(1,(d-1)))
  mBIC_g
  
  
  BIC_tree=do.call(c,lapply(1:(d-1),function(i)sum(VineTree[[i]]$E$BIC)))
  BIC_total=sum(BIC_tree[2:(d-1)])
  #BICPropSubTrees=BIC_total*Prop_BIC # total from T_2 to T_d-1.
  BIC_cum=cumsum(BIC_tree[-1])
  TruncLevel_BICmin=which.min(BIC_cum)+1
  #TruncLevel_BICprop=which.min(abs(cumsum(BIC_tree[-1])) <= abs(BICPropSubTrees))+1
  TruncLevelDelta=min(TruncLevel_BICmin,TruncLevel_Ind)
  
  if(BIC_graph){
    BIC_g = ggplot() +
      geom_point(aes(x = 2:(d-1), y = BIC_cum)) +
      geom_vline(xintercept = TruncLevelDelta, linetype=3, size=1.5) +
      geom_vline(xintercept = TruncLevel_BICmin, linetype=4, size=1.5) +
      labs(title="",x=expression("Truncation Level (q)"), y = "BIC values") +
      theme(
        axis.text.x = element_text(color = "black",face="bold", size=14),
        axis.text.y = element_text(face="bold", size = 14),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                        colour = "grey"), 
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "grey"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
      scale_x_continuous(limits = c(1,(d-1)))
    BIC_g
  }
  
  ##  Store conditioned sets and conditioning sets in maximum spanning trees 
  VineGraph <- list()
  ParList <- list()
  FamList <- list()
  for(i in 1:(d-1)){
    if(i==1){
      VineGraph[[1]] <- list('C'=do.call(rbind,MST[[1]]$E$conditionedSet),'D'=NULL)
      ParList[[1]]=do.call(cbind,VineTree[[1]]$E$Copula.param)[1,]
      FamList[[1]]=VineTree[[1]]$E$Copula.type
    }else{
      VineGraph[[i]] <- list('C'=do.call(rbind,MST[[i]]$E$conditionedSet),'D'=do.call(rbind,MST[[i]]$E$conditioningSet))
      ParList[[i]]=do.call(cbind,VineTree[[i]]$E$Copula.param)[1,]
      FamList[[i]]=VineTree[[i]]$E$Copula.type
    }
  }
  ##  Construct a structure matrix using edge information
  StrMtx_MST=edgesToM(VineGraph, format = T)
  
  dstar=ifelse(test = trunclevel,trunclevel-1,d_col-2) # if trunclevel=F, i ranges from 1:d_col-2 b/c no order in the deepest edge 
  Edges <- append(list(t(VineGraph[[1]][[1]])), lapply(VineGraph[-1], function(i) rbind(t(i[[1]]),t(i[[2]]))))
  Edges_str <- list()
  for(i in 1:dstar){
    if(i==1){
      Edges_str[[1]]=rbind(StrMtx_MST[1,-1],diag(StrMtx_MST[-1,-1]))  
    }else{
      Edges_str[[i]]=rbind(StrMtx_MST[1:i,-(1:i)],diag(StrMtx_MST[-(1:i),-(1:i)]))
    }
  }
  
  for(i in 1:dstar){
    # Find the column permutation
    Reorder <- apply(Edges_str[[i]], 2, function(col_orig) {
      match_result <- apply(Edges[[i]], 2, function(col_perm) setequal(col_orig,col_perm))
      which(match_result)
    })
    if(i==1){
      #VineGraph[[1]]=list('C'=VineGraph[[1]]$C[Reorder,],'D'=NULL)  
      ParList[[1]]=ParList[[1]][Reorder]
      FamList[[1]]=FamList[[1]][Reorder]
    }else{
      #VineGraph[[i]]=list('C'=VineGraph[[i]]$C[Reorder,],'D'=as.matrix(VineGraph[[i]]$D[Reorder,]))  
      ParList[[i]]=ParList[[i]][Reorder]
      FamList[[i]]=FamList[[i]][Reorder]
    }
  }
  
  ##  Specify a parameter matrix
  all.pairs <- combn(1:d_col, 2)
  Param=do.call(c,ParList)
  ParMtx_MST=matrix(0,d_col,d_col)
  for(i in 1:length(Param)){
    ParMtx_MST[all.pairs[1,i],all.pairs[2,i]]=Param[i]  
  }
  
  ##  Specify a family matrix
  Fam=do.call(c,FamList)
  FamMtx_MST=matrix(0,d_col,d_col)
  for(i in 1:length(Fam)){
    FamMtx_MST[all.pairs[1,i],all.pairs[2,i]]=Fam[i]  
  }
  
  ##  Define a XVineSpec() including permuted matrices
  XVS_spec=XVineSpec(M = StrMtx_MST, Mmod = FamMtx_MST, Mpar = ParMtx_MST)
  #for(i in 1:d){
  #  print(VineCopula::RVineMatrixCheck(XVS$xmat[,,i]))
  #}
  
  FittedDat_P=ParetoSim(n = N, XVS = XVS_spec) # Pareto scale (no need to relabel nodes b/c the ft 'XVineSim' rearrange them in ascending order)
  
  if(Rank_chiL){
    emp_chimat <- ChiMatrixMC(data) # already rank-based samples
    emp_chimat <- emp_chimat[upper.tri(emp_chimat)]
    XVine_chimat <- ChiMatrixMC(FittedDat_P,Quan = qt)
    XVine_chimat <- XVine_chimat[upper.tri(XVine_chimat)]
  }else{
    # Model-based chi
    Dat_P=ParetoSim(n = N, XVS = XVS)
    emp_chimat <- ChiMatrixMC(1/Dat_P) # samples from multivariate Pareto dist
    emp_chimat <- emp_chimat[upper.tri(emp_chimat)]
    XVine_chimat <- ChiMatrixMC(1/FittedDat_P)
    XVine_chimat <- XVine_chimat[upper.tri(XVine_chimat)]
  }
  ChiLPlot <- ggplot() +
    geom_point(aes(x = c(emp_chimat),y = c(XVine_chimat))) +
    geom_abline(slope = 1, intercept = 0) +
    labs(title="",x=expression(paste("Model"," ",chi)), y = expression(paste("Fitted"," ",chi))) +
    theme(
      axis.text.x = element_text(color = "black",face="bold", size=14),
      axis.text.y = element_text(face="bold", size =14),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      panel.background = element_rect(fill = "white",
                                      colour = "white",
                                      linewidth = 0.5, linetype = "solid"),
      panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                      colour = "grey"), 
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                      colour = "grey"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    scale_x_continuous(n.breaks = 5)
  
  ChiLPlot
  
  if(Rank_chiU){
    emp_chimat <- emp_chi(data,p = 1-qt)
    emp_chimat <- emp_chimat[upper.tri(emp_chimat)]
    XVine_chimat <- emp_chi(FittedDat_P,p = 1-qt)
    XVine_chimat <- XVine_chimat[upper.tri(XVine_chimat)]
  }else{
    emp_chimat <- emp_chi(ParetoSim(n = N, XVS = XVS))
    emp_chimat <- emp_chimat[upper.tri(emp_chimat)]
    XVine_chimat <- emp_chi(FittedDat_P)
    XVine_chimat <- XVine_chimat[upper.tri(XVine_chimat)]
  }
  #chiXVine <- Gamma2chi(emp_vario(data = FittedDat_P)) # unless it is a pure HR case, not use
  ChiUPlot <- ggplot() +
    geom_point(aes(x = c(emp_chimat),y = c(XVine_chimat))) +
    geom_abline(slope = 1, intercept = 0) +
    labs(title="",x=expression(paste("Model"," ",chi)), y = expression(paste("Fitted"," ",chi))) +
    theme(
      axis.text.x = element_text(color = "black",face="bold", size=14),
      axis.text.y = element_text(face="bold", size =14),
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      panel.background = element_rect(fill = "white",
                                      colour = "white",
                                      linewidth = 0.5, linetype = "solid"),
      panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                      colour = "grey"), 
      panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                      colour = "grey"),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
    scale_x_continuous(n.breaks = 5)
  
  ChiUPlot
  
  if(Chi3_graph){
    if(Rank_chi3){
      emp_chi3mat <- EmpChiArray(data)
      emp_chi3mat <- emp_chi3mat[upper.tri(emp_chi3mat[,,1])]
      XVine_chi3mat <- EmpChiArray(FittedDat_P,Quan = qt)
      XVine_chi3mat <- XVine_chi3mat[upper.tri(XVine_chi3mat[,,1])]
    }else{
      emp_chi3mat <- EmpChiArray(1/ParetoSim(n = N, XVS = XVS))
      emp_chi3mat <- emp_chi3mat[upper.tri(emp_chi3mat[,,1])]
      XVine_chi3mat <- EmpChiArray(1/FittedDat_P)
      XVine_chi3mat <- XVine_chi3mat[upper.tri(XVine_chi3mat[,,1])]
    }
    Chi3Plot <- ggplot() +
      geom_point(aes(x = c(emp_chi3mat),y = c(XVine_chi3mat))) +
      geom_abline(slope = 1, intercept = 0) +
      labs(title="",x=expression(paste("Empirical"," ",chi)), y = expression(paste("Fitted"," ",chi))) +
      theme(
        axis.text.x = element_text(color = "black",face="bold", size=14),
        axis.text.y = element_text(face="bold", size =14),
        axis.title.x = element_text(color="black", size=14, face="bold"),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linewidth = 0.5, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                        colour = "grey"), 
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "grey"),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
      scale_x_continuous(n.breaks = 5)
    Chi3Plot
  }
  
  
  return(list("MST"=MST,"VineTree"=VineTree,"mBIC_g"=mBIC_g,"BIC_g"=BIC_g,"XVS_spec"=XVS_spec,"emp_chimat"=emp_chimat,"XVine_chimat"=XVine_chimat,"ChiLPlot"=ChiLPlot,"ChiUPlot"=ChiUPlot,"TruncLevelStar"=TruncLevelStar,"TruncLevel_mBICmin"=TruncLevel_mBICmin))
}

