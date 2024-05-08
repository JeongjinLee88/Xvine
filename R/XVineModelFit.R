#' Fit X-vine models
#' 
#' @description 
#' Fits X-Vine models to a d-dimensional inverted-Pareto data set. The function selects vine tree structures, using maximum spanning tree.
#' and select bivariate exponent measure classes and pair-copula classes and estimate parameters for the corresponding classes
#' Users can specify the type of edge weights for the first tree and subsequent trees separately and the class of bivariate exponent measures and pair-copulas
#' Note that if you do not use samples from the limiting distribution, then the `Rank` transformation is required to create chi-plots or variogram-plots.
#' That is, if you fit the X-vine model to the real dataset, then the rank transformation is required throughout.
#' 
#' @param data A \eqn{n\times d} data matrix of multivariate Pareto samples.
#' @param Rank Logical; whether rank transformation is performed or not (\code{Rank=T}; default).
#' @param Rank_chiL Logical; whether rank-based samples are used or not for creating the chi-plot of the lower pairwise tail dependence measure.
#' If \code{Rank_chiL=FALSE}, then `ChiLPlot` compares model-based chi's with fitted chi's via Monte Carlo simulation.
#' @param Rank_chi3 Logical; whether rank-based samples are used or not for creating the chi-plot of the trivariate tail dependence measure.
#' If \code{Rank_chi3=FALSE}, then `ChiLPlot` compares model-based chi's with fitted chi's via Monte Carlo simulation.
#' @param Chi3_graph Logical; whether the chi-plot of the trivariate tail dependence measure is plotted or not.
#' @param Rank_Vario Logical; whether rank-based samples are used or not for the plot of the empirical pairwise variogram versus the fitted pairwise variogram. 
#' @param Vario_graph Logical; whether the plot of the empirical pairwise variogram versus the fitted pairwise variogram is created or not.
#' @param MST1_HR Logical; the minimum spanning tree for the Husler-Reiss model is plotted or not.
#' @param quan Numeric; a lower threshold for the rank transformation. It switches from Pareto scale to uniform scale.
#' @param N Numeric; the sample size to draw samples from the limiting distribution of multivariate Pareto distribution.
#' @param XVS A list consisting of three components: reconstructed structure matrices, family matrices, parameter matrices, see:[XVineSpec()].
#' @param tcfamset Numeric vector; the class of bivariate exponent measures.
#' @param pcfamset Numeric vector; the class of bivariate pair-copula models.
#' @param selectioncrit Character string; indicates the selection criteria for the class of bivariate exponent measures or pair-copulas.
#' @param BIC_graph Logical; whether the BIC graph over the tree level is plotted or not.
#' @param trunclevel Numeric; indicates the specified truncation level (\code{trunclevel=NULL}; default).
#' @param progress Logical; whether the progress of selecting vine tree structures is printed.
#' @param treecritT1 A character string, indicating the tree criterion for the first tree.
#' @param treecritT2 A character string, indicating the tree criterion for subsequent trees.
#' @param si Numeric; a tuning parameter for mBIC (\eqn{\si_0=0.9}; default).
#' @param effsampsize Numeric; the specified effective sample size for the independence copula (\eqn{n_{D_e}}<10; default).
#' @param tau_threshold Numeric; the specified Kendall's tau value for the independence copula (\eqn{\hat{\tau}_e}<0.05; default)
#' @param se Logical; whether standard errors for ML estimators are reported.
#' @param weights Logical; whether weights should be assigned to observations when missing values exist.
#' @param cores Numeric; indicates the number of cores for parallel computing (optional).
#'
#' @return A nested list object containing:
#' * MST: Maximum spanning trees
#' * VineTree: Fitted X-vine models
#' * mBIC_g: the plot of mBIC values across tree levels
#' * BIC_g: the plot of BIC values across tree levels
#' * XVS_spec: the X-vine specification from fitted X-vine models
#' * emp_chimat: the matrix of empirical pairwise chi's
#' * XVine_chimat: the matrix of fitted pairwise chi's
#' * ChiLPlot: the chi-plot of the lower pairwise tail dependence measure
#' * Chi3Plot: the chi-plot of the trivariate tail dependence measure
#' * VarioPlot: the plot of the pairwise empirical varigoram versus fitted variogram
#' * TruncLevelStar: the optimal truncation level determined as the minimum between the tree level with the lowest mBIC value
#'  and the tree level such that all pair-copulas are set to independence copula
#' * TruncLevel_mBICmin: the truncation level with the lowest mBIC value
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
#' Dat_P=ParetoSim(n = 2000, XVS = XVS) # Pareto scale
#' XVineFitOut=XVineModelFit(data = 1/Dat_P, N = 2000, XVS = XVS
#'                          , Rank = FALSE, Rank_chiL = FALSE, Rank_chi3 = FALSE
#'                          , Chi3_graph = TRUE, Rank_Vario = FALSE, Vario_graph = FALSE
#'                          , tcfamset = c(1:4), pcfamset = c(0,1,3,4,5,6,13,14,16)
#'                          , selectioncrit = "AIC", trunclevel = FALSE, progress = TRUE
#'                          , treecritT1 = "chi", treecritT2 = "tau"
#'                          , effsampsize = 10, tau_threshold = 0.05
#'                          , weights = NA, cores=1)
XVineModelFit <- function(data, Rank=TRUE, Rank_chiL=FALSE, Rank_Vario=FALSE, Rank_chi3=TRUE
                          , Chi3_graph=FALSE, Vario_graph=FALSE, MST1_HR=FALSE, quan=0.2, N=2000
                          , XVS, tcfamset = c(1,2,3,4), pcfamset = c(0,1,3,4,5,6,13,14,16), selectioncrit = "AIC"
                          , BIC_graph=TRUE, trunclevel = FALSE, progress = TRUE
                          , treecritT1 = "chi", treecritT2= 'tau', si=0.9, effsampsize=10, tau_threshold=0.05, se=FALSE, weights=NA, cores = 1)
{
  ##  Note that the 'XVineModelFit' function uses multivariate 'inverted' Pareto samples
  ##  If you directly use samples from the limiting distribution, they must be ones from Pareto distribution with Pareto margin.
  if(Rank){
    data <- ParetoTransRank(data = data, u_quan = quan, scaleType = "U") 
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
    MST[[1]] <- MST_HR(Dat_Pareto = data,quan = quan)
  }else{
    g <- InitializeFirstGraph(data, T1crit, weights)
    MST[[1]] <- findMST(g, mode = "RVine")
  }
  VineTree[[1]] <- fit.FirstTree(MST[[1]], data, tcfamset, si=si,
                                 selectioncrit = selectioncrit)
  if(d > 2){
    for (tree in 2:(d - 1)) { # tree=2,3,4
      g <- BuildNextGraph(VineTree[[tree-1]], weights, treecrit = Ticrit, truncated = FALSE)
      MST[[tree]] <- findMST(g, mode = "RVine", truncated = FALSE)
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
  if(!is.null(XVS)){
    Dat_Pa=ParetoSim(n = N, XVS = XVS)  
  }else{
    if(any(c(Rank_chiL,Rank_Vario,Rank_chi3)==FALSE))
    stop("Either 'XVS' must be specified or all ('Rank_chiL','Rank_Vario','Rank_chi3') must be 'TRUE'")
  }
  
  if(Rank_chiL){
    # Empirical chi's vs Fitted chi's
    emp_chimat <- ChiMtxMC(data) # already rank-based samples
    emp_chimat <- emp_chimat[upper.tri(emp_chimat)]
    XVine_chimat <- ChiMtxMC(FittedDat_P,quan = quan)
    XVine_chimat <- XVine_chimat[upper.tri(XVine_chimat)]
  }else{
    # Model-based chi's vs Fitted chi's via Monte Carlo simulation
    emp_chimat <- ChiMtxMC(1/Dat_Pa) # limiting samples from the specified inverted-MPD
    emp_chimat <- emp_chimat[upper.tri(emp_chimat)]
    XVine_chimat <- ChiMtxMC(1/FittedDat_P)
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
  
  if(Vario_graph){ # use samples on Pareto scale
    if(Rank_Vario){
      emp_Variomat <- emp_vario(data = data,p = 1-quan)
      emp_Variomat <- emp_Variomat[upper.tri(emp_Variomat)]
      XVine_Variomat <- emp_vario(FittedDat_P,p = 1-quan)
      XVine_Variomat <- XVine_Variomat[upper.tri(XVine_Variomat)]
    }else{
      emp_Variomat <- emp_vario(Dat_Pa)
      emp_Variomat <- emp_Variomat[upper.tri(emp_Variomat)]
      XVine_Variomat <- emp_chi(FittedDat_P)
      XVine_Variomat <- XVine_Variomat[upper.tri(XVine_Variomat)]
    }
    VarioPlot <- ggplot() +
      geom_point(aes(x = c(emp_Variomat),y = c(XVine_Variomat))) +
      geom_abline(slope = 1, intercept = 0) +
      labs(title="",x=expression(paste("Empirical"," ",Gamma)), y = expression(paste("Fitted"," ",Gamma))) +
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
    VarioPlot
  }else{
    VarioPlot=NULL
  }
  
  if(Chi3_graph){
    if(Rank_chi3){
      emp_chi3mat <- EmpChiArray(data)
      emp_chi3mat <- emp_chi3mat[upper.tri(emp_chi3mat[,,1])]
      XVine_chi3mat <- EmpChiArray(FittedDat_P,quan = quan)
      XVine_chi3mat <- XVine_chi3mat[upper.tri(XVine_chi3mat[,,1])]
    }else{
      emp_chi3mat <- EmpChiArray(1/Dat_Pa)
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
  }else{
    Chi3Plot=NULL
  }
  
  
  return(list("MST"=MST,"VineTree"=VineTree,"mBIC_g"=mBIC_g,"BIC_g"=BIC_g
              ,"XVS_spec"=XVS_spec,"emp_chimat"=emp_chimat,"XVine_chimat"=XVine_chimat
              ,"ChiLPlot"=ChiLPlot,"Chi3Plot"=Chi3Plot
              ,"VarioPlot"=VarioPlot,"TruncLevelStar"=TruncLevelStar
              ,"TruncLevel_mBICmin"=TruncLevel_mBICmin))
}

