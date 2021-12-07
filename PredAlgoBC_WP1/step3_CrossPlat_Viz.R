

#!/usr/bin/Rscript

#=========================================================================
#=================================================================================


#         STEP3  vizualisation and beginning of metric calculation


#===========================================================================
#=====================================================================


setwd("tcga")

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(dendextend)
library(FactoMineR)

# prepare to call file  =>  we can't load all because it is heavy and too slow

fileNames <- list.files(full.names = T,pattern = "exprSet_",ignore.case=TRUE)
fileNames <- fileNames[grep("_array|_Rsq", invert = T,fileNames)]  # remove not transforemd
fileNames <- gsub("\\./|\\.rds","",fileNames)  # remove.rds to have nicer names
saveRDS(fileNames,"fileNames.rds")   

# should get :
# exprSet_noCPN
# exprSet_TDM
# exprSet_FSQN
# exprSet_RBE
# exprSet_Comb
# exprSet_MM
# exprSet_GQ
# exprSet_Zsc


sampleAnnot <- readRDS("sampleAnnot_noCPN.rds")
condCol<-"subtype"   
batchCol<-"batch" 

#=====================

# function

#===================



# adapted gPCA function from https://rdrr.io/cran/gPCA/man/gPCA.batchdetect.html
#=======================

# attention, initial function has been modified to have batch information in numeric format.
# If batch is not numeric,  it gives wrong results!!!!!



gPCA.batchdetect <-function (x, batch, filt = NULL, nperm = 1000, center = FALSE,
                             scaleY = FALSE, seed = NULL)
{
  if (!is.null(seed)) {
    set.seed(seed)
  }
  batch = as.numeric(as.factor(batch))
  permute <- matrix(NA, ncol = length(batch), nrow = 50000)
  for (j in 1:50000) {
    permute[j, ] <- sample(batch, replace = FALSE)
  }
  samp <- sample(1:dim(permute)[1], nperm, replace = FALSE)
  permute.samp <- permute[samp, ]
  if (center == FALSE) {
    x2 <- scale(x, center = T, scale = F)
  }
  else {
    x2 <- x
  }
  if (sum(is.na(x)) > 0) {
    missing <- readline(prompt = "Missing values detected. Continue with mean value imputation? (Note this may take a very long time, but it will automatically save in your working dir so you don't have to ever run it again.) [y/n] ")
    if (substr(missing, 1, 1) == "n") {
      stop("The PC cannot be calculated with missing values.")
    }
    else {
      x2.imp <- ifelse(is.na(x2), rowMeans(x2, na.rm = TRUE),
                       x2)
      save(x2.imp, "x2.imputed.RData")
    }
  }
  else {
    x2.imp <- x2
  }
  if (is.null(filt)) {
    data.imp <- x2.imp
  }
  else {
    sd <- apply(x2.imp, 2, sd)
    rank <- rank(sd)
    keep <- (1:length(sd))[rank %in% (length(rank) - filt +
                                        1):length(rank)]
    data.imp <- x2.imp[, keep]
  }
  n <- dim(data.imp)[1]
  p <- dim(data.imp)[2]
  b <- length(unique(batch))
  n
  p
  b
  if (length(batch) != n) {
    stop("Matrices do not conform: length(batch)!=n")
  }
  y <- matrix(nrow = length(batch), ncol = length(unique(batch)))
  for (j in 1:length(unique(batch))) {
    y[, j] <- ifelse(batch == j, 1, 0)
  }
  if (scaleY == FALSE) {
    y2 <- scale(y, center = T, scale = F)
  }
  else {
    ys <- matrix(nrow = length(batch), ncol = length(unique(batch)))
    nk <- apply(y, 2, sum)
    for (j in 1:length(unique(batch))) {
      ys[, j] <- ifelse(batch == j, 1/nk[j], 0)
    }
    y2 <- scale(ys, center = F, scale = F)
  }
  svd.x <- svd(data.imp)
  PC.u <- data.imp %*% svd.x$v
  var.x <- var(PC.u)
  varPCu1 <- diag(var.x)[1]/sum(diag(var.x))
  cumulative.var.u <- numeric()
  for (i in 1:dim(var.x)[1]) {
    cumulative.var.u[i] <- sum(diag(var.x)[1:i])/sum(diag(var.x))
  }
  svd.bat <- svd(t(y2) %*% data.imp)
  PC.g <- data.imp %*% svd.bat$v
  var.bat <- var(PC.g)
  varPCg1 <- diag(var.bat)[1]/sum(diag(var.bat))
  cumulative.var.g <- numeric()
  for (i in 1:dim(var.bat)[1]) {
    cumulative.var.g[i] <- sum(diag(var.bat)[1:i])/sum(diag(var.bat))
  }
  delta <- diag(var.bat)[1]/diag(var.x)[1]
  delta.p <- numeric()
  for (i in 1:nperm) {
    batch.p <- permute.samp[i, ]
    y <- ys <- matrix(nrow = length(batch.p), ncol = length(unique(batch.p)))
    for (j in 1:length(unique(batch.p))) {
      y[, j] <- ifelse(batch.p == j, 1, 0)
    }
    if (scaleY == FALSE) {
      y2 <- scale(y, center = T, scale = F)
    }
    else {
      nk <- apply(y, 2, sum)
      for (j in 1:length(unique(batch.p))) {
        ys[, j] <- ifelse(batch.p == j, 1/nk[j], 0)
      }
      y2 <- scale(ys, center = F, scale = F)
    }
    svd.bat.p <- svd(t(y2) %*% data.imp)
    var.bat.p <- var(data.imp %*% svd.bat.p$v)
    PC.g.p <- diag(var.bat.p)[1]/sum(diag(var.bat.p))
    delta.p[i] <- diag(var.bat.p)[1]/diag(var.x)[1]
  }
  p.val <- sum(delta < delta.p)/length(delta.p)
  p.val
  p.val <- ifelse(p.val == 0, "<0.001", round(p.val,
                                              3))
  out <- list(delta = delta, p.val = p.val, delta.p = delta.p,
              batch = batch, filt = filt, n = n, p = p, b = b, PCg = PC.g,
              PCu = PC.u, varPCu1 = varPCu1, varPCg1 = varPCg1, nperm = nperm,
              cumulative.var.u = cumulative.var.u, cumulative.var.g = cumulative.var.g)
}


# dendrogram plot 
#==================

#cuttree by pair

cutree_bypair <-function(dend, talk= F){
  require("dendextend")
  labels_node <- dend %>%get_nodes_attr("label")
  results <-list()
  for (i in 1:(length(labels_node)-1)){
    
    if ( if  (i==1){
      !is.na(labels_node[[i+1]])&!is.na(labels_node[[i]])
    } else{
      is.na(labels_node[[i-1]])&!is.na(labels_node[[i]])& !is.na(labels_node[[i+1]])
      
    }
    )     {
      if (talk==T){print(paste0("1er IF is true - loop ",i))}
      results[[i]] <- c(i,i)
      names(results[[i]] ) <- c(labels_node[[i]],labels_node[[i+1]])
    } else {
      results[[i]]  <- 0
      names(results[[i]] ) <- c(labels_node[[i]])
      if (talk==T){print(paste0("1er IF is false - loop ",i))}
    }
  }
  results <-unlist(results)
  results <-results[!is.na(names(results))]
  results <-data.frame(keyName=names(results), cluster=results, row.names=NULL) #to keep num num
  results <- results[results$cluster!=0,]  #new code to test...
  # results <-aggregate(.~keyName, results,  FUN = sum)   #...instead of this code but made fake positive with double 0
  return(results)
}

colDendPlot= function(col_choice,sampleAnnot,clusters, plotTitle,labels_cex){
  variable <- as.data.frame(col_choice)
  variable$ID<-row.names(sampleAnnot)
  dend <- as.dendrogram(clusters)
  ordered_names <- as.data.frame(dendextend::cutree(dend, 1, order_clusters_as_data = FALSE))
  ordered_names$ID<-rownames(ordered_names)
  dend_col<-merge(ordered_names,variable,by="ID", sort=FALSE)  #join est sur plyr
  dend_col<-dend_col$col_choice
  dend<-hang.dendrogram(dend, hang=0.1)
  dendCol<-rainbow(nlevels(dend_col))
  names(dendCol)<-levels(dend_col)
  colorsdend<-dendCol[dend_col]
  labels_colors(dend)<-colorsdend
  labels_cex(dend)<-labels_cex
  #hang.dendrogram(dend, hang = 0.1, cex=0.3,main=plotTitle)
  plot(dend, cex=0.5,main=plotTitle)
}


#================================================
#===============================


#  ALL anlaysis, dataset by dataset (to avoid to lad everything because it's very heavy)


#==========================
#================================================

corPltf_all <- list()     # correlation stat
goodmatch<-list()      # dendro stat
corcoef_3D <-list()   # PCA stat
delta_gPCA_all <-list()    # gPCA stat


for (i in fileNames){
  
  exprSet_toplot <- readRDS(paste0(i,".rds"))
  
  
  
  #=======================

  # correlation

  #======================



  nbSpl <- ncol(exprSet_toplot[[1]])/2

  corPltf2 <- list()
  for (j in names(exprSet_toplot)){
    corPltf <- list()
    for (k in (1:nbSpl)){
      corPltf[[k]]  <- cor(exprSet_toplot[[j]][,k],exprSet_toplot[[j]][,nbSpl+ k])

    }

    corPltf2[[j]] <-unlist(corPltf)
    print(paste0(j,  ": corr calculation done"))
  }

  corPltf_all[[i]] <-as.data.frame(do.call(cbind, corPltf2))
  names(corPltf_all[[i]])<-paste0(i,"_",names(corPltf_all[[i]]))
  rownames(corPltf_all[[i]]) <- colnames(exprSet_toplot[[1]])[(nbSpl+1):(nbSpl*2)]



  #=======================

  # dendrogramm

  #======================


  pdf(paste0("dend_aftCPN__",i,".pdf"),width = 20, height = 12)
  par(mfrow=c(2,1))

  nb_mismatch <- list()
  perc_goodmach  <- list()

  for (k in names(exprSet_toplot)){

    distance <- dist(t(exprSet_toplot[[k]]),method="euclidian")
    clusters <- hclust(distance)

    # calculate pair matching
    dend <- as.dendrogram(clusters)
    cutree_cell <-cutree_bypair(dend, talk=F)
    rownames(cutree_cell)<-cutree_cell$keyName
    cutree_cell$keyName <-gsub("_array|_arr","",cutree_cell$keyName)
    cutree_cell$dupli <-duplicated(cutree_cell) | duplicated(cutree_cell, fromLast = TRUE)
    nb_mismatch[[k]] <- (length ((cutree_cell$dupli)[cutree_cell$dupli == FALSE]))/2
    nb_goodmatch <- (length ((cutree_cell$dupli)[cutree_cell$dupli == TRUE]))/2
    perc_goodmach[[k]] <- round(nb_goodmatch*100/(nrow(sampleAnnot[[1]]) /2),digits = 1)


    colDendPlot= function(col_choice,DendsampleAnnot,clusters, plotTitle,labels_cex){
      variable <- as.data.frame(col_choice)
      variable$ID<-row.names(DendsampleAnnot)
      dend <- as.dendrogram(clusters)
      ordered_names <- as.data.frame(dendextend::cutree(dend, 1, order_clusters_as_data = FALSE))
      ordered_names$ID<-rownames(ordered_names)
      dend_col<-merge(ordered_names,variable,by="ID", sort=FALSE)  #join est sur plyr
      dend_col<-dend_col$col_choice
      dend<-hang.dendrogram(dend, hang=0.1)
      dendCol<-rainbow(nlevels(dend_col))
      names(dendCol)<-levels(dend_col)
      colorsdend<-dendCol[dend_col]
      labels_colors(dend)<-colorsdend
      labels_cex(dend)<-labels_cex
      #hang.dendrogram(dend, hang = 0.1, cex=0.3,main=plotTitle)
      plot(dend, cex=0.5,main=plotTitle)
    }


    # plot dendro

    for (j in c(batchCol,condCol) ){
      col_choice<- as.factor(sampleAnnot[[1]][,j])
      colDendPlot(col_choice,DendsampleAnnot=sampleAnnot[[1]],clusters,
                  plotTitle=paste0(k,"-",i," by ", j),labels_cex=0.5)
    }

    print(paste0(i,"-", k ,": dendrogram done"))

  } # k end

  goodmatch[[i]] <- data.frame(nb_mismatch=unlist(nb_mismatch),
                               perc_goodmach= unlist(perc_goodmach),
                               row.names=paste0(i,"_", names(nb_mismatch)) )

  dev.off()


  
  #====================
  
  # PCA
  
  #====================
  
  pdf(paste0("PCA_aftCPN__",i,".pdf"),width = 6, height = 5)
  
  corcoef <- list()
  for (k in names(exprSet_toplot)){
    
    allintable <-cbind(t(exprSet_toplot[[k]]),sampleAnnot[[k]][,c(batchCol,condCol)])
    pca<- FactoMineR:: PCA(allintable, scale.unit=TRUE, ncp=5 ,
                           quali.sup= (ncol(allintable)-1):ncol(allintable)  , graph=F)
    p_batchCol <- FactoMineR:: plot.PCA(pca, axes=c(1, 2), choix="ind", 
                                        habillage=(ncol(allintable)-1),label="none",title = paste0(k," - ", i) )  #plot colored accorded to batch1 (A and B) +barycentre
    p_condCol <- FactoMineR::plot.PCA(pca, axes=c(1, 2), choix="ind", 
                                      habillage=ncol(allintable),label="none",title = paste0(k," - ", i)) 
    corcoef[[k]] <-dimdesc(pca, axes=c(1,2,3), proba = 1)   #proba no selected to have all the number
    print(p_batchCol)
    print(p_condCol)
  }
  
  dev.off()
  

  corcoef_dim1 <- list()
  corcoef_dim2 <- list()
  corcoef_dim3 <- list()
  for (j in names(exprSet_toplot)){

      corcoef_dim1 [[j]] <- corcoef[[j]]$Dim.1$quali
      corcoef_dim1 [[j]] <- round(corcoef_dim1 [[j]],digits = 3)
      rownames(corcoef_dim1 [[j]]) <- paste("dim1" ,i,j,rownames(corcoef_dim1[[j]]),  sep = "_")
      
      corcoef_dim2[[j]] <- corcoef[[j]]$Dim.2$quali
      corcoef_dim2 [[j]] <- round(corcoef_dim2 [[j]],digits = 3)
      rownames(corcoef_dim2 [[j]]) <- paste("dim2" ,i,j,rownames(corcoef_dim2[[j]]),  sep = "_")
      
      corcoef_dim3[[j]] <- corcoef[[j]]$Dim.3$quali
      corcoef_dim3 [[j]] <- round(corcoef_dim3 [[j]],digits = 3)
      rownames(corcoef_dim3 [[j]]) <- paste("dim3" ,i,j,rownames(corcoef_dim3[[j]]),  sep = "_")
    }
    
    
    for (k in names(corcoef_dim1)){
      corcoef_dim1[[k]]<- corcoef_dim1[[k]][order(rownames(corcoef_dim1[[k]])),]
      corcoef_dim2[[k]]<- corcoef_dim2[[k]][order(rownames(corcoef_dim2[[k]])),]
      corcoef_dim3[[k]]<- corcoef_dim3[[k]][order(rownames(corcoef_dim3[[k]])),]
    }
    
    corcoef_dim1<- as.data.frame(do.call("rbind", corcoef_dim1))
    corcoef_dim2<- as.data.frame(do.call("rbind", corcoef_dim2))
    corcoef_dim3<- as.data.frame(do.call("rbind", corcoef_dim3))    
    
    
    # prepare table corr
    corcoef_dims <- cbind(corcoef_dim1[,"R2"], corcoef_dim2[,"R2"],corcoef_dim3[,"R2"]) 
    rownames(corcoef_dims) <- gsub("dim1_","",rownames(corcoef_dim1))
    colnames(corcoef_dims)<- c("PC1", "PC2", "PC3")
    
    corcoef_dims_batch<- corcoef_dims[grep(paste0("_",batchCol), rownames(corcoef_dims)),] 
    rownames(corcoef_dims_batch) <-  gsub(paste0("_",batchCol),"",rownames(corcoef_dims_batch))
    corcoef_dims_cond<- corcoef_dims[grep(paste0("_",condCol), rownames(corcoef_dims)),] 
    rownames(corcoef_dims_cond) <-  gsub(paste0("_",condCol),"",rownames(corcoef_dims_cond))
    
  
    df <- cbind(as.data.frame.table(corcoef_dims_cond, responseName = "cor.cond"), cor.batch = c(corcoef_dims_batch))
    
    
    p<-ggplot(df, aes(x=Var2, y=Var1, col = cor.batch, size = cor.cond)) + 
      geom_point() + ggtitle(i) +
      scale_color_gradientn(limits=c(0,1),colors =c("#5183cf","#cf2b0a"),
                            breaks=c(0,0.5,1),labels=c("0","0.5","1") ) +           
      scale_size(name = waiver(), limits=c(0,1)) +     #mettre la size_legend à la scale désirée
      xlab("") + ylab("") + theme_minimal()
    
    pdf(paste0("PCAcorcoeff_aftCPN__",i,".pdf"),width = 5, height = 4)
    print(p)
    dev.off()

    
    corcoef_3D[[i]] <- list(PC1=corcoef_dim1,PC2=corcoef_dim2,PC3=corcoef_dim3)
  
    print ("PCA done")
  
  #==========================
  
  #   gPCA
  
  #=============================
  
  
  delta_gPCA <- do.call(rbind,lapply(exprSet_toplot, function(x) {
    tryCatch({
      out_batch<-gPCA.batchdetect(t(x),batch=as.numeric(as.factor(sampleAnnot[[1]][,batchCol])),center=FALSE,
                                  scaleY=FALSE,filt=NULL,nperm=100,seed=5)
      pvalue_batch <- out_batch$p.val
      delta_batch <- out_batch$delta
      out_sub<-gPCA.batchdetect(t(x),batch=as.numeric(as.factor(sampleAnnot[[1]][,condCol])),center=FALSE,
                                scaleY=FALSE,filt=NULL,nperm=100,seed=5)
      pvalue_sub <- out_sub$p.val
      delta_sub <- out_sub$delta
    }, error=function(e){cat("Error",conditionMessage(e), "\n")})
    return(data.frame(delta_batch = delta_batch ,pval_batch = pvalue_batch,delta_sub=delta_sub,pval_sub=pvalue_sub,stringsAsFactors = FALSE))
  }
  ))
  
  rownames(delta_gPCA) <-paste0(i,"_",rownames(delta_gPCA))
  delta_gPCA_all[[i]] <-delta_gPCA
 
  
  
  
  
  
  
   print(paste0("=====================================  ",i,  ": done"))
  
  
  
  
  #=============================
  
  # collect results at the end
  
  #=============================
  
  
} 


write.csv(goodmatch,"dend_goodmatch.csv")    # dendro stat
saveRDS(goodmatch, file="goodmatch.rds")      # dendro stat

saveRDS(corPltf_all, file="corPltf_all.rds")   # correlation stat

saveRDS(corcoef_3D, file="corcoef_3D.rds")    # PCA stat

delta_gPCA_all <-as.data.frame(do.call(rbind, delta_gPCA_all))
rownames(delta_gPCA_all) <- gsub(".*\\.", "", rownames(delta_gPCA_all))  # remove what is before the period
saveRDS(delta_gPCA_all, file= "delta_gPCA_all.rds")    # gPCA stat




# the end   => see step 4 for metric calulation scaluing and heatmap


