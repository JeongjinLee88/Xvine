XVineModelFit_lightversion <- function(data, Rank=TRUE, MST1_HR=FALSE, quan=0.2, BIC_graph=TRUE
                          , tcfamset = c(1,2,3,4), pcfamset = c(0,1,3,4,5,6,13,14,16), selectioncrit = "AIC"
                          , trunclevel = FALSE, progress = TRUE
                          , treecritT1 = "chi", treecritT2= 'tau', si=0.9, effsampsize=10, tau_threshold=0.05, se=FALSE, weights=NA, cores = 1)
{
  ##  Note that the 'XVineModelFit' function uses multivariate 'inverted' Pareto samples
  ##  If you directly use samples from the limiting distribution, they must be ones from Pareto distribution with Pareto margin.
  if(Rank){
    Dat_U <- ParetoTransRank(data = data, u_quan = quan, scaleType = "U") 
  }else{
    Dat_U <- 1/data #switch to Uniform scale from Pareto scale
  }
  # Determine the truncation level if specified.
  d <- ifelse(test = trunclevel,trunclevel+1,ncol(Dat_U))
  d_col <- ncol(Dat_U)
  n <- nrow(Dat_U)
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
    g <- InitializeFirstGraph(Dat_U, T1crit, weights)
    MST[[1]] <- findMST(g, mode = "RVine")
  }
  VineTree[[1]] <- fit.FirstTree(MST[[1]], Dat_U, tcfamset, si=si,
                                 selectioncrit = selectioncrit)
  if(d > 2){
    for (tree in 2:(d - 1)) { # tree=2,3,4
      g <- BuildNextGraph(VineTree[[tree-1]], weights, treecrit = Ticrit, truncated = FALSE)
      MST[[tree]] <- findMST(g, mode = "RVine", truncated = FALSE)
      VineTree[[tree]] <- fit.SubTree(data = Dat_U, MST = MST, VineTree = VineTree, copfamset = pcfamset, tree = tree, si=si,
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
  
  return(list("MST"=MST,"VineTree"=VineTree,"mBIC_g"=mBIC_g,"BIC_g"=BIC_g
              ,"TruncLevelStar"=TruncLevelStar
              ,"TruncLevel_mBICmin"=TruncLevel_mBICmin))
}

