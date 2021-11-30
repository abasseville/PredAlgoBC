#!/usr/bin/Rscript

#=========================================================================
#=================================================================================


#         vizualisation and beginning of metric calculation


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

# should get :
# exprSet_noCPN.rds
# exprSet_TDM.rds
# exprSet_FSQN.rds
# exprSet_RBE.rds
# exprSet_Comb.rds
# exprSet_MM.rds
# exprSet_GQ.rds
# exprSet_XPN.rds

fileNames <-gsub(".rds","",fileNames)   # remove.rds to have nicer names

sampleAnnot <- readRDS("sampleAnnot_noCPN.rds")
condCol<-"subtype"   
batchCol<-"batch" 

#=====================

# function

#===================



# custom gPCA function
#=======================


# from https://rdrr.io/cran/gPCA/man/gPCA.batchdetect.html

# attention, I change batch  to be numeric : if not, it gives wrong results!!!!!

gPCA.batchdetect <-function (x, as.numeric(as.factor(batch)), filt = NULL, nperm = 1000, center = FALSE,
                             scaleY = FALSE, seed = NULL)
{
  if (!is.null(seed)) {
    set.seed(seed)
  }
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

corPltf3 <- list()     # correlation stat

goodmatch<-list()      # dendro stat

all_corcoef <-list()   # PCA stat

delta_gPCA <-list()    #gPCA stat



for (i in 1:length(fileNames)){

  exprSet_toplot <- readRDS(paste0(fileNames,".rds"))
  
  
#=======================  
  
  # correlation
  
#======================
  
  
    
    nbSpl <- ncol(exprSet_toplot[[1]])/2
    
    corPltf2 <- list()
    for (j in names(exprSet_toplot)){
      corPltf <- list()
      for (i in (1:nbSpl)){
        corPltf[[i]]  <- cor(exprSet_toplot[[j]][,i],exprSet_toplot[[j]][,nbSpl+ i])
        
      }
      
      corPltf2[[j]] <-unlist(corPltf)
      print(paste0(j,  ": done"))
    }
    
    print(paste0(k,  ": done"))
    corPltf3[[i]] <-as.data.frame(do.call(cbind, corPltf2))
    names(corPltf3[[i])<-paste0(fileNames,"_",names(corPltf3[[i]]))
    rownames(corPltf3[[i]]) <- colnames(exprSet_toplot[[1]])[(nbSpl+1):(nbSpl*2)]
    
  
  
  #=======================  
  
  # dendrogramm
  
  #====================== 
  
    
  pdf(paste0("dend_aftCPN_",fileNames,".pdf"),width = 20, height = 12)
    par(mfrow=c(3,1))
    
  nb_mismatch <- list()
  perc_goodmach  <- list()
    
  for (k in names(exprSet_toplot)){

  distance <- dist(t(exprSet_toplot[[k]]),method="euclidian")
  clusters <- hclust(distance)

  # calculate pair matching
    dend <- as.dendrogram(clusters)
    cutree_cell <-cutree_bypair(dend, talk=F)
    rownames(cutree_cell)<-cutree_cell$keyName
    cutree_cell$keyName <-gsub("_affy","",cutree_cell$keyName)
    cutree_cell$dupli <-duplicated(cutree_cell) | duplicated(cutree_cell, fromLast = TRUE)
    nb_mismatch[[k]] <- (length ((cutree_cell$dupli)[cutree_cell$dupli == FALSE]))/2
    nb_goodmatch <- (length ((cutree_cell$dupli)[cutree_cell$dupli == TRUE]))/2
    perc_goodmach[[k]] <- round(nb_goodmatch*100/(nrow(sampleAnnot[[1]]) /2),digits = 1)


  
  # plot dendro
    
  for (j in c(batchCol,condCol) ){
    col_choice<- as.factor(sampleAnnot[[1]][,j])
    colDendPlot(col_choice,sampleAnnot[[1]],clusters, 
                plotTitle=paste0(k,"-",name_toplot,"dataset by ", j),labels_cex=0.5)
  }
    
    print(paste0(fileNames,"-", k ,": OK"))
    
  } # k end 
    
  goodmatch[[i]] <- data.frame(nb_mismatch=unlist(nb_mismatch), perc_goodmach= unlist(perc_goodmach), row.names=gsub(".rds","", filenames))
    
  dev.off() 
  


#====================

# PCA

#====================
  
pdf(paste0("PCA_aftCPN_",fileNames,".pdf"),width = 12, height = 10)
  
  corcoef <- list()
  for (k in names(exprSet_toplot)){
  
  allintable <-cbind(t(exprSet_toplot[[k]]),sampleAnnot[[1]][,c(batchCol,condCol)])
  pca<- FactoMineR:: PCA(allintable, scale.unit=TRUE, ncp=5 ,
                              quali.sup= (ncol(allintable)-1):ncol(allintable)  , graph=F)
  FactoMineR:: plot.PCA(pca, axes=c(1, 2), choix="ind", habillage=(ncol(allintable)-1),label="none",title = paste0(k," - ", i) )  #plot colored accorded to batch1 (A and B) +barycentre
  FactoMineR::plot.PCA(pca, axes=c(1, 2), choix="ind", habillage=ncol(allintable),label="none",title = paste0(k," - ", i)) 
  corcoef[[k]] <-dimdesc(pca, axes=c(1,2,3), proba = 1)   #proba no selected to have all the number
  print(paste0(i,"_",k , " : PCA done"))
}
all_corcoef[[i]] <- corcoef

dev.off()


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
  
  rownames(delta_gPCA) <-paste0(fileNames,"_",rownames(delta_gPCA))
  delta_gPCA[[i]] <-delta_gPCA
  print(paste0(i,  ": done"))
  



#=============================

# collect results at the end

#=============================


} 


delta_gPCA <-as.data.frame(do.call(rbind, delta_gPCA))
saveRDS(delta_gPCA, file= "delta_gPCA.rds")

write.csv(goodmatch,"dend_goodmatch.csv")
saveRDS(goodmatch, file="goodmatch.rds")

# the end see step 4



