


#!/usr/bin/Rscript

#===============================================================================================
#=========================================================================================

#       STEP 2 : Cross-platform normalization

#=========================================================================================
#===============================================================================================


# load library 
#===========================

library(TDM)      #remotes::install_github("greenelab/TDM")
library(limma)   # RBE
library(sva)   # Combat
library(tidyr)
library(plyr)
library(ggplot2)
   

# violinplot  function
#===================

violinPlot = function (df, mainTitle,dotsize,binwidth,Xangle =0,ordered_x=unique(dataf$data)  ){
  dataf <- gather(gene_df ,key="Data", value="Val")
  dataf$Data <- factor(dataf$Data , levels=ordered_x)
  ggplot(dataf, aes(x=Data, y=Val, fill=Data)) +
    theme(axis.text.x = element_text(angle = Xangle, hjust = 1),
          axis.title.x=element_blank()) +
    ggtitle(mainTitle)+
    scale_x_discrete(limits=names(df))+    #to avoid ggplot to reorder alph automatiqualy  
    geom_violin(trim = FALSE)+
    geom_dotplot(binaxis='y', stackdir='center',dotsize=dotsize,fill = "black",binwidth = binwidth)
}

violinPlot_Ylim = function (df, mainTitle="",dotsize,binwidth,Xangle =0,ylo,yhi,ordered_x=unique(dataf$data)){
  dataf <- gather(df ,key="Data", value="Val")
  dataf$Data <- factor(dataf$Data , levels=ordered_x)
  ggplot(dataf, aes(x=Data, y=Val, fill=Data)) +
    theme(axis.text.x = element_text(angle = Xangle, hjust = 1),
          axis.title.x=element_blank())+
    ggtitle(mainTitle)+ ylim(ylo,yhi)+
    scale_x_discrete(limits=names(df))+    #to avoid ggplot to reorder alph automatiqualy  
    geom_violin(trim = FALSE) + 
    geom_dotplot(binaxis='y', stackdir='center',dotsize=dotsize,fill = "black",binwidth = binwidth)
}

# normalization custom function
#===================

# These functions were copied from  https://github.com/dy16b/Cross-Platform-Normalization 
# publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7868049/


quantileNormalizeByFeature2 <- function (matrix_to_normalize, target_distribution_matrix) 
{
  if (ncol(matrix_to_normalize) != ncol(target_distribution_matrix)) {
    cat("ERROR: Data matrices are not compatible - column lengths differ!")
  }
  else {
    data.qn <- matrix(0, nrow = nrow(matrix_to_normalize), 
                      ncol = ncol(matrix_to_normalize))
    for (i in 1:ncol(matrix_to_normalize)) {
      feature.to.normalize <- matrix_to_normalize[, i]
      target.feature.dist <- target_distribution_matrix[, 
                                                        i]
      result <- preprocessCore::normalize.quantiles.use.target(x = as.matrix(feature.to.normalize), 
                                                               target = target.feature.dist, copy = TRUE)
      data.qn[, i] <- result
    }
    rownames(data.qn) = rownames(matrix_to_normalize)
    colnames(data.qn) = colnames(matrix_to_normalize)
    return(data.qn)
  }
}

processplatforms = function(datalist, namesvec=NULL, skip.match=FALSE){
  #Convert data from various formats to the proper format for use 
  #with all the crossnorm normalization functions
  
  for(i in 1:length(datalist)){
    if(is.matrix(datalist[[i]])){
      datalist[[i]] <- as.data.frame(datalist[[i]])
    }
  }
  
  if (is.null(namesvec)){
    namesvec <- numeric(length(datalist))
    for (i in 1:length(datalist)){
      namesvec[i] <- 0
    }
  }
  
  #Put the row names in their places
  for (i in 1:length(namesvec)){
    if(namesvec[i] != 0){
      rownames(datalist[[i]]) = datalist[[i]][,namesvec[i]]
      datalist[[i]] = datalist[[i]][,-1*namesvec[i],drop=FALSE]
    }	
  }
  
  if(!skip.match){
    #Create the common genes list
    commongenes <- rownames(datalist[[1]])
    for (i in 2:length(datalist)){
      commongenes <- intersect(commongenes,rownames(datalist[[i]]))
    }
    
    
    #Put it all together
    for (i in 1:length(datalist)){
      datalist[[i]] <- datalist[[i]][commongenes,,drop=FALSE]
    }
  }
  return(datalist)
}

OLS <- function(X, Y){
  ## A fast row-wise univariate OLS regression based on linear algebra
  Xbar <- rowMeans(X); Ybar <- rowMeans(Y, na.rm = TRUE)
  X.c <- sweep(X, 1, Xbar); Y.c <- sweep(Y, 1, Ybar)
  CovXY <- rowSums(X.c * Y.c, na.rm = TRUE)
  VarX <- rowSums(X.c^2, na.rm = TRUE); VarY <- rowSums(Y.c^2, na.rm = TRUE)
  ## regression coefs.
  beta1 <- CovXY / VarX; beta0 <- Ybar - beta1 * Xbar
  betamat <- cbind("Intercept"=beta0, "Slope"=beta1)
  ## Pearson corr. coefs.
  cc <- CovXY / sqrt(VarX * VarY)
  ## predictions
  Yhat <- sweep(X, 1, beta1, FUN="*") + beta0 %o% rep(1, ncol(X))
  ## RSS vector
  RSS <- rowSums((Y - Yhat)^2)
  SS1 <- rowSums(Yhat^2)
  SST <- rowSums(Y^2)
  return(list(betamat=betamat, corr=cc, covxy=CovXY, varx=VarX, Yhat=Yhat, RSS=RSS, var.explained.prop=SS1/SST, Xbar=Xbar, Ybar=Ybar))
}


flmer <- function(Xmat, Ymat){
  Xmat <- as.matrix(Xmat); Ymat <- as.matrix(Ymat)
  ## 11/12/2018. Use new notations; p==ngenes; n==sample size
  p <- nrow(Xmat); n <- ncol(Xmat); N <- p*n
  ## calculate some useful items
  Xibar <- rowMeans(Xmat); Xbar <- mean(Xibar)
  Xmat.c <- Xmat - Xibar
  Yibar <- rowMeans(Ymat); Ybar <- mean(Yibar)
  Ymat.c <- Ymat - Yibar
  covXY <- rowSums(Xmat.c * Ymat.c)
  covXX <- rowSums(Xmat.c * Xmat.c)
  beta.yixi <- covXY / covXX
  ## The overall regression
  Xc <- Xmat - mean(Xmat); Yc <- Ymat - mean(Ymat)
  beta.yx <- sum(Xc*Yc) / sum(Xc^2)
  ## In this case, we assume  that there is no collinearity in Z
  qprime <- 2*p
  var.epsilon <- sum((Ymat.c - Xmat.c * beta.yixi)^2) / (N - qprime)
  ## Preparing for calculating S, \| \mathbf{R}_{Z|\mathbf{X}} \|^2,
  ## and \|Z' \mathbf{R}_{Z|\mathbf{X}} \|^2 terms.
  XiYibar <- covXY/n + Xibar*Yibar
  ## Xi2bar = \|Xi\|^2 / n
  Xi2bar <- covXX/n + Xibar^2
  ## Xs is Xmat standardized by the global mean/sd
  Xs <- (Xmat - Xbar) / sqrt(sum((Xmat-Xbar)^2))
  Xsbar <- rowMeans(Xs)
  ## XsX is the inner producd between Xs and X. It is denoted as
  ## \varsigma in the notes.
  XsX <- rowSums(Xs*Xmat)
  XsNorm2 <- sum(Xs^2)
  ## the normalized version of S
  S <- n^2 * (sum((Yibar - Ybar -(Xibar - Xbar)*beta.yx)^2) + sum((XiYibar -Ybar*Xibar -(Xi2bar - Xbar*Xibar)*beta.yx)^2)) / var.epsilon
  ## Calculate \| \mathbf{R}_{Z|\mathbf{X}} \|^2
  Projxz.norm2 <- n + XsNorm2 * ( n^2*sum(Xsbar^2) + sum(XsX^2)) + n*sum(Xibar^2)/p
  Rzx.norm2 <- sum(Xmat^2) + N -Projxz.norm2
  ## General terms; without n^4 yet.
  S1 <- sum(Xibar^2); S2 <- sum(Xsbar^2); S3 <- sum(XsX^2)/n^2
  ZPZ.norm2 <- p^2/N^2 + S2^2 + 2*( p*S1/N^2 + S2*S3) + S1^2/N^2 + 2*(sum(Xibar*XsX)^2)/(N*n^2) + S3^2
  ## subtractions; to multiply by *n^4
  Minus.terms <- sum((1/N + Xsbar^2)^2) + 2*sum((Xibar/N + Xsbar*XsX/n)^2) + sum((Xibar^2/N + XsX^2/n^2)^2)
  ## added terms
  XiNorm2 <- rowSums(Xmat^2)
  Add.terms <- sum((1/n -1/N -Xsbar^2)^2) + 2*sum((Xibar/n -Xibar/N -Xsbar*XsX/n)^2) + sum((XiNorm2/n^2 -Xibar^2/N -XsX^2/n^2)^2)
  ZRzx.norm2 <- n^4*(ZPZ.norm2 - Minus.terms + Add.terms)
  ## Estimate lambda based on Moment matching
  lambdahat <- max(0, (S - Rzx.norm2) / ZRzx.norm2)
  ## Inference of the fixed effects. First, compute some small
  ## "building blocks".
  Axx <- XiNorm2 - (lambdahat/(1+n*lambdahat)) * n^2 * Xibar^2
  A1x <- n*Xibar/(1+n*lambdahat)
  A1y <- n*Yibar/(1+n*lambdahat)
  Axy <- XiYibar*n - (lambdahat/(1+n*lambdahat)) * n^2 * Xibar* Yibar
  W11 <- n/(1+n*lambdahat) - lambdahat*A1x^2/(1+lambdahat*Axx)
  W1x <- A1x/(1+lambdahat*Axx)
  W1y <- A1y -lambdahat*A1x*Axy / (1+lambdahat*Axx)
  Wxx <- Axx/(1+lambdahat*Axx)
  Wxy <- Axy / (1+lambdahat*Axx)
  ## Now the actual estimation
  covBeta <- var.epsilon * solve(matrix(c(sum(W11), sum(W1x), sum(W1x), sum(Wxx)), 2))
  betahat <- drop(covBeta %*% c(sum(W1y), sum(Wxy)) / var.epsilon)
  names(betahat) <- c("Intercept", "Slope")
  dimnames(covBeta) <- list(c("Intercept", "Slope"), c("Intercept", "Slope"))
  ## Mixed effects terms via EBLUP
  gamma0hat <- lambdahat * (W1y - betahat[1]*W11 - betahat[2]*W1x)
  gamma1hat <- lambdahat * (Wxy - betahat[1]*W1x - betahat[2]*Wxx)
  ## Individual betas
  betamat <- t(rbind(gamma0hat, gamma1hat) + betahat)
  colnames(betamat) <- c("Intercept", "Slope")
  ## t-statistics for the fixed effects
  t.fixed <- betahat/sqrt(diag(covBeta))
  ## the predicted Yhat
  Yhat <- Xmat * betamat[, "Slope"] + betamat[, "Intercept"]
  return(list(betahat=betahat, betamat=betamat, Yhat=Yhat,
              lambdahat=lambdahat, var.epsilon=var.epsilon,
              covBeta=covBeta, t.fixed=t.fixed))
}


MM <- function(Xmat, Ymat){
  rr1 <- OLS(Xmat, Ymat)
  cov1 <- cov(rr1$betamat)
  s0 <- sqrt(cov1[1,1]); s1 <- sqrt(cov1[2,2]); rho <- cov1[1,2]/s0/s1
  a12 <- rho/sqrt(1-rho^2); a22 <- s1/(s0 * sqrt(1-rho^2))
  A <- matrix(c(1, a12,
                0, a22), 2, byrow=TRUE)
  ## the X transformation
  Xtilde <- a12 + a22*Xmat
  ## apply flmer() to the covariance transformed data
  rr3 <- flmer(Xtilde, Ymat)
  ## the reverse transformation. Note that each *row* of "betamat" is
  ## beta_i; so we have to transpose the reverse matrix multiplication.
  betamat <- rr3$betamat %*% t(A); colnames(betamat) <- colnames(rr3$betamat)
  ## other misc. items
  betahat <- drop(A %*% rr3$betahat); names(betahat) <- names(rr3$betahat)
  lambdahat <- rr3$lambdahat
  covGamma <-  rr3$lambdahat * rr3$var.epsilon * (A %*% t(A))
  Yhat <- rr3$Yhat
  var.epsilon <- rr3$var.epsilon
  covBeta <- A %*% rr3$covBeta %*% t(A)
  t.fixed <- betahat/sqrt(diag(covBeta))
  return(list(betahat=betahat, betamat=betamat, Yhat=Yhat,
              lambdahat=lambdahat, var.epsilon=var.epsilon,
              covGamma <- covGamma,
              covBeta=covBeta, t.fixed=t.fixed))
}


gq = function(platform1.data, platform2.data, p1.names=0, p2.names=0, skip.match=FALSE){
  #This function is basically a wrapper for normalizeGQ
  
  #Match names
  input = processplatforms(list(x=platform1.data,y=platform2.data),namesvec = c(p1.names, p2.names), skip.match=skip.match)
  
  #Prepare for normalizeGQ
  combined = cbind(input$x,input$y)
  pf = c(seq(1,1,length.out=dim(input$x)[2]),seq(2,2,length.out=dim(input$y)[2]))
  
  #Call normalizeGQ
  ngq = normalizeGQ(combined,pf)
  
  #Split the results and return
  out=split(seq(pf),pf)
  out[[1]] = ngq[,out[[1]]]
  out[[2]] = ngq[,out[[2]]]
  names(out) <- c("x","y")
  return(out)
}


normalizeGQ <- function(M, pf, ...) { 
  #This function was provided by Xiao-Qin Xia, one of the authors of webarraydb modified MRS
  # M is the data matrix
  # pf is the vector to specify the platform for each column of M.
  idx <- split(seq(pf), pf)
  if (length(pf)<=1) return(M)
  imax <- which.max(sapply(idx, length)) # for reference
  ref_med <- apply(M[, idx[[imax]]], 1, function(x) median(x, na.rm=TRUE))
  ref_med_srt <- sort(ref_med)
  idx[imax] <- NULL
  lapply(idx, function(i) {
    MTMP <- sapply(i, function(x) ref_med_srt[rank(M[,x])]); 
    M[,i] <<- MTMP - apply(MTMP, 1, median) + ref_med 
  } )
  invisible(M)
}


# zscore function
#====================

# function inspiration from from hamza
# Pool  <- sapply(noms_cohortes , function(coh) { exprSet <- get(coh);t(apply(exprSet,1 , function(x) (x-median(x))/sd(x)))})


geneTrans <- function(exprMat)  {     #exprmat = gene in row, sample in col
  results <-t(apply(exprMat,1 , function(x) (x-median(x))/sd(x)))
  return(results)
}

BatchGeneTrans <- function(dat,batch){
  if (length(dim(batch)) > 1) {
    stop("only allows one batch variable")
  }
  batch <- as.factor(batch)
  design <- model.matrix(~-1 + batch)
  n.batch <- nlevels(batch)
  batches <- list()
  for (i in 1:n.batch) {
    batches[[i]] <- which(batch == levels(batch)[i])
  }
  n.batches <- sapply(batches, length)
  if (any(n.batches == 1)) {
    stop("Note: one batch has only one sample")
  }
  n.array <- sum(n.batches)
  check <- apply(design, 2, function(x) all(x == 1))
  design <- as.matrix(design[, !check])
  
  for (i in 1:length(batches) ) {
    selrow <- batches[[i]]
    dat[, selrow] <- geneTrans(dat[, selrow] )
  }
  return(dat)
}


#========================================================================================


# load dataset


#======================================================================================

setwd("tcga")

exprSet <- readRDS("exprSet.rds")   # from STEP1
exprSet <- exprSet[-1]   # remove "count" dataset from the list (1st one)
sampleAnnot <- readRDS("sampleAnnot.rds")

sampleAnnot_All <- list()
for (i in names(exprSet) ){   
  sampleAnnot_All[[i]] <- sampleAnnot
  sampleAnnot_All[[i]]$batch <- i
}

BatchCol<-"batch" 


#==========================================================

#     Merging each RNAseq table with array

#=========================================================


exprSet_array <- exprSet[["array"]]
exprSet_Rseq <- exprSet[-length(exprSet)]  # remove the last one = "array" one
sampleAnnot_array<- sampleAnnot_All[["array"]]
sampleAnnot_Rseq <- sampleAnnot_All[-length(exprSet)]

# select randomly 4 genes and 4 samples to plot  for each cross pltform normalization
random_genes <-sample(row.names(exprSet_array),size = 4,replace =F)                                
random_sples<-sample(names(exprSet_array[,1:ncol(exprSet_array)]),size = 4,replace =F)
                     
                     
#===============================================
                     
 # no crossplatform normalization (noCPN) reference dataset
                     
#======================================================
                     
                     
CombSampleAnnot <-list()
for (i in names(sampleAnnot_Rseq) ){
sampleAnnot_array2 <-sampleAnnot_array
rownames(sampleAnnot_array2) <- paste0(rownames(sampleAnnot_array2),"_arr")
CombSampleAnnot[[i]] <- rbind(sampleAnnot_array2,sampleAnnot_Rseq[[i]])
CombSampleAnnot[[i]]$subtype <- "Lum"
CombSampleAnnot[[i]]$subtype[CombSampleAnnot[[i]]$er.status.by.ihc == "Negative" & CombSampleAnnot[[i]]$pr.status.by.ihc == "Negative" & CombSampleAnnot[[i]]$her2.status.by.ihc == "Negative"]<-"TripleNeg"
CombSampleAnnot[[i]]$subtype[CombSampleAnnot[[i]]$her2.status.by.ihc == "Positive"]<-"HER2"
 }
                     
                     
exprSet_noCPN<-list()
 for (i in names(exprSet_Rseq ) ){
print(paste0("nrow before merge:",nrow(exprSet_Rseq[[i]])))
print(paste0("ncol before merge:", ncol(exprSet_Rseq[[i]]), "  expected after: ", (ncol(exprSet_Rseq[[i]])*2)))
exprSet_noCPN[[i]] <- merge(exprSet_array, exprSet_Rseq[[i]],
                           by= "row.names")
row.names(exprSet_noCPN[[i]]) <-exprSet_noCPN[[i]][,1]; exprSet_noCPN[[i]] <-exprSet_noCPN[[i]][,-1]
names(exprSet_noCPN[[i]])<- gsub (".x$","_arr",names(exprSet_noCPN[[i]]),ignore.case = FALSE)
names(exprSet_noCPN[[i]])<- gsub (".y$","",names(exprSet_noCPN[[i]]),ignore.case = FALSE)
print(paste0("nrow after merge:",nrow(exprSet_noCPN[[i]])))
print(paste0("ncol after merge:",ncol(exprSet_noCPN[[i]])))
print("-------------------------")
}
                     
 #  violin plot
#===============
                     
                
exprSet_toplot <- exprSet_noCPN
name_toplot <-"noCPN"
n <-ncol(exprSet_toplot[[1]])  # same length for all df
                     
pdf(paste0("after_CPN_violinPlot_4genes_",name_toplot,".pdf"),height = 4, width = 6)       
                     
for (j in random_genes){
                       
gene_df <-do.call(rbind, (lapply(exprSet_toplot, function(x) x[j,((n/2)+1):n])))   #because "NOT array samples" are at the second part of the dataframe
gene_df<- as.data.frame(t(gene_df))
gene_df$array <-as.numeric(  exprSet_toplot[[1]][j,1:(n/2)]  )  # because array samples are at the first part of the dataframe
p<-violinPlot(gene_df, mainTitle=paste0(name_toplot," - gene = ",j),
                                     dotsize=0,binwidth=0.2,Xangle=45, ordered_x = names(gene_df))      
 print (p)
                       
}
dev.off()
                     
pdf(paste0("after_CPN_violinPlot_4sples_",name_toplot,".pdf"),height = 4, width = 6)   
                     
for (j in random_sples){
        
denst <-lapply( lapply( exprSet_toplot, function(y)  y=y[, ((n/2)+1):n]) , 
                                       function(x) density(x[,j])  )  
                       
                       ind_df <-do.call(rbind, (lapply(denst, function(dens) dens$y)  ) )
                       ind_df<- as.data.frame(t(ind_df))
                       densarray <- density(exprSet_array[,paste0(j)])
                       ind_df$array <- densarray$y
p<-violinPlot_Ylim(ind_df, mainTitle=paste0(name_toplot, " - pat = ",j),
                     dotsize=0,binwidth=0.2,Xangle=45, ylo = -0.1, yhi = 0.5 ,
                     ordered_x = names(ind_df))      
print (p)
}
dev.off()
                     
                     
saveRDS(CombSampleAnnot,file="sampleAnnot_noCPN.rds")
saveRDS(exprSet_noCPN,file="exprSet_noCPN.rds")
                     
                     
#========================
                     
# TDM
                     
#==================
                     
 exprSet_TDM <-list()
 for  (i in names(exprSet)[-length(exprSet)] ) {      # no transformation for the last of the list:  array
   exprSet_TDM [[i]]<-TDM::tdm_transform(ref_data = data.table(cbind(gene=rownames(exprSet[["array"]]), 
                                                                     exprSet[["array"]])),
                                         target_data = data.table(cbind(gene=rownames(exprSet[[i]]), 
                                                                        exprSet[[i]])),
                                         log_target = TRUE, tdm_factor = 1)
   
   # transform to numeric (not numric with datatable)
   exprSet_TDM[[i]] <- (as.data.frame(exprSet_TDM[[i]]) )
   row.names.tdm<-exprSet_TDM [[i]]$gene
   exprSet_TDM [[i]] <- as.data.frame(sapply(exprSet_TDM [[i]][,-1], as.numeric) )
   row.names(exprSet_TDM[[i]]) <- row.names.tdm
   
   # add array
   exprSet_TDM[[i]] <- merge(exprSet_array, exprSet_TDM[[i]],
                               by= "row.names")
   row.names(exprSet_TDM[[i]]) <-exprSet_TDM[[i]][,1]; exprSet_TDM[[i]] <-exprSet_TDM[[i]][,-1]
   names(exprSet_TDM[[i]])<- gsub (".x$","_array",names(exprSet_TDM[[i]]),ignore.case = FALSE)
   names(exprSet_TDM[[i]])<- gsub (".y$","",names(exprSet_TDM[[i]]),ignore.case = FALSE)
   
   print(paste0(i, ": done"))
 }
 
 #  violin plot
 #===============                          
 
 exprSet_toplot <- exprSet_TDM
 name_toplot <-"TDM"
 n <-ncol(exprSet_toplot[[1]])
 
 pdf(paste0("after_CPN_violinPlot_4genes_",name_toplot,".pdf"),height = 4, width = 6)       
 
 for (j in random_genes){
   
   gene_df <-do.call(rbind, (lapply(exprSet_toplot, function(x) x[j,((n/2)+1):n])))   #because "NOT array samples" are at the second part of the dataframe
   gene_df<- as.data.frame(t(gene_df))
   gene_df$array <-as.numeric(  exprSet_toplot[[1]][j,1:(n/2)]  )  # because array samples are at the first part of the dataframe
   p<-violinPlot(gene_df, mainTitle=paste0(name_toplot," - gene = ",j),
                 dotsize=0,binwidth=0.2,Xangle=45, ordered_x = names(gene_df))      
   print (p)
   
 }
 dev.off()
 
 pdf(paste0("after_CPN_violinPlot_4sples_",name_toplot,".pdf"),height = 4, width = 6)   
 
 for (j in random_sples){
   denst <-lapply( lapply( exprSet_toplot, function(y)  y=y[, ((n/2)+1):n]) , 
                   function(x) density(x[,j])  )  
   
   ind_df <-do.call(rbind, (lapply(denst, function(dens) dens$y)  ) )
   ind_df<- as.data.frame(t(ind_df))
   densarray <- density(exprSet_array[,paste0(j)])
   ind_df$array <- densarray$y
   p<-violinPlot_Ylim(ind_df, mainTitle=paste0(name_toplot, " - pat = ",j),
                      dotsize=0,binwidth=0.2,Xangle=45, ylo = -0.1, yhi = 0.5 ,
                      ordered_x = names(ind_df))      
   print (p)
 }
 dev.off()
 
 saveRDS(exprSet_TDM, file= "exprSet_TDM.rds")
 rm(exprSet_TDM)
                     
                     
                     
#===================
                     
 #     FSQN
                     
#===================
 
 
 exprSet_FSQN <-list()
 for  (i in names(exprSet)[-length(exprSet)] ) {
   exprSet_FSQN[[i]] <- quantileNormalizeByFeature2(as.matrix(exprSet[[i]]),   #test
                                                    as.matrix(exprSet[["array"]]))   #target
   
   # add array
   exprSet_FSQN[[i]] <- merge(exprSet_array, exprSet_FSQN[[i]],
                             by= "row.names")
   row.names(exprSet_FSQN[[i]]) <-exprSet_FSQN[[i]][,1]; exprSet_FSQN[[i]] <-exprSet_FSQN[[i]][,-1]
   names(exprSet_FSQN[[i]])<- gsub (".x$","_array",names(exprSet_FSQN[[i]]),ignore.case = FALSE)
   names(exprSet_FSQN[[i]])<- gsub (".y$","",names(exprSet_FSQN[[i]]),ignore.case = FALSE)
   print(paste0("exprSet ", i, ": done"))
 }
 
 #  violin plot
 #===============                          
 
 exprSet_toplot <- exprSet_FSQN
 name_toplot <-"FSQN"
 n <-ncol(exprSet_toplot[[1]])
 
 pdf(paste0("after_CPN_violinPlot_4genes_",name_toplot,".pdf"),height = 4, width = 6)       
 
 for (j in random_genes){
   
   gene_df <-do.call(rbind, (lapply(exprSet_toplot, function(x) x[j,((n/2)+1):n])))   #because "NOT array samples" are at the second part of the dataframe
   gene_df<- as.data.frame(t(gene_df))
   gene_df$array <-as.numeric(  exprSet_toplot[[1]][j,1:(n/2)]  )  # because array samples are at the first part of the dataframe
   p<-violinPlot(gene_df, mainTitle=paste0(name_toplot," - gene = ",j),
                 dotsize=0,binwidth=0.2,Xangle=45, ordered_x = names(gene_df))      
   print (p)
   
 }
 dev.off()
 
 pdf(paste0("after_CPN_violinPlot_4sples_",name_toplot,".pdf"),height = 4, width = 6)   
 
 for (j in random_sples){
   denst <-lapply( lapply( exprSet_toplot, function(y)  y=y[, ((n/2)+1):n]) , 
                   function(x) density(x[,j])  )  
   
   ind_df <-do.call(rbind, (lapply(denst, function(dens) dens$y)  ) )
   ind_df<- as.data.frame(t(ind_df))
   densarray <- density(exprSet_array[,paste0(j)])
   ind_df$array <- densarray$y
   p<-violinPlot_Ylim(ind_df, mainTitle=paste0(name_toplot, " - pat = ",j),
                      dotsize=0,binwidth=0.2,Xangle=45, ylo = -0.1, yhi = 0.5 ,
                      ordered_x = names(ind_df))      
   print (p)
 }
 dev.off()
 
 saveRDS(exprSet_FSQN, file= "exprSet_FSQN.rds")
 rm(exprSet_FSQN)
 
 
 #===================
 
 #     RBE  
 
 #===================
 
 
 exprSet_RBE <-list()
 for (i in names(exprSet_noCPN)) {
   exprSet_RBE[[i]]<-limma::removeBatchEffect(exprSet_noCPN[[i]],batch = CombSampleAnnot[[i]][,BatchCol])
   exprSet_RBE[[i]]<-exprSet_RBE[[i]]-min(exprSet_RBE[[i]])
   print(paste0(i, " done"))
 }
 
 #  violin plot
 #===============                          
 
 exprSet_toplot <- exprSet_RBE
 name_toplot <-"RBE"
 n <-ncol(exprSet_toplot[[1]])
 
 pdf(paste0("after_CPN_violinPlot_4genes_",name_toplot,".pdf"),height = 4, width = 6)       
 
 for (j in random_genes){
   
   gene_df <-do.call(rbind, (lapply(exprSet_toplot, function(x) x[j,((n/2)+1):n])))   #because "NOT array samples" are at the second part of the dataframe
   gene_df<- as.data.frame(t(gene_df))
   gene_df$array <-as.numeric(  exprSet_toplot[[1]][j,1:(n/2)]  )  # because array samples are at the first part of the dataframe
   p<-violinPlot(gene_df, mainTitle=paste0(name_toplot," - gene = ",j),
                 dotsize=0,binwidth=0.2,Xangle=45, ordered_x = names(gene_df))      
   print (p)
   
 }
 dev.off()
 
 pdf(paste0("after_CPN_violinPlot_4sples_",name_toplot,".pdf"),height = 4, width = 6)   
 
 for (j in random_sples){
   denst <-lapply( lapply( exprSet_toplot, function(y)  y=y[, ((n/2)+1):n]) , 
                   function(x) density(x[,j])  )  
   
   ind_df <-do.call(rbind, (lapply(denst, function(dens) dens$y)  ) )
   ind_df<- as.data.frame(t(ind_df))
   densarray <- density(exprSet_array[,paste0(j)])
   ind_df$array <- densarray$y
   p<-violinPlot_Ylim(ind_df, mainTitle=paste0(name_toplot, " - pat = ",j),
                      dotsize=0,binwidth=0.2,Xangle=45, ylo = -0.1, yhi = 0.5 ,
                      ordered_x = names(ind_df))      
   print (p)
 }
 dev.off()                                   
 
 saveRDS(exprSet_RBE, file= "exprSet_RBE.rds")
 rm(exprSet_RBE)
 
 
 
 #=============================
 
 # Combat   #sva dependant
 
 #=============================
 
 exprSet_Comb <-list()
 for (i in names(exprSet_noCPN)){
   exprSet_Comb[[i]]<-sva::ComBat(as.matrix(exprSet_noCPN[[i]]),batch = CombSampleAnnot[[i]][,BatchCol],
                                  mod=NULL, par.prior = TRUE, prior.plots = FALSE)
   print(paste0(i, ": done"))
 }
 
 
 #  violin plot
 #===============                          
 
 exprSet_toplot <- exprSet_Comb
 name_toplot <-"ComBat"
 n <-ncol(exprSet_toplot[[1]])
 
 pdf(paste0("after_CPN_violinPlot_4genes_",name_toplot,".pdf"),height = 4, width = 6)       
 
 for (j in random_genes){
   
   gene_df <-do.call(rbind, (lapply(exprSet_toplot, function(x) x[j,((n/2)+1):n])))   #because "NOT array samples" are at the second part of the dataframe
   gene_df<- as.data.frame(t(gene_df))
   gene_df$array <-as.numeric(  exprSet_toplot[[1]][j,1:(n/2)]  )  # because array samples are at the first part of the dataframe
   p<-violinPlot(gene_df, mainTitle=paste0(name_toplot," - gene = ",j),
                 dotsize=0,binwidth=0.2,Xangle=45, ordered_x = names(gene_df))      
   print (p)
   
 }
 dev.off()
 
 pdf(paste0("after_CPN_violinPlot_4sples_",name_toplot,".pdf"),height = 4, width = 6)   
 
 for (j in random_sples){
   denst <-lapply( lapply( exprSet_toplot, function(y)  y=y[, ((n/2)+1):n]) , 
                   function(x) density(x[,j])  )  
   
   ind_df <-do.call(rbind, (lapply(denst, function(dens) dens$y)  ) )
   ind_df<- as.data.frame(t(ind_df))
   densarray <- density(exprSet_array[,paste0(j)])
   ind_df$array <- densarray$y
   p<-violinPlot_Ylim(ind_df, mainTitle=paste0(name_toplot, " - pat = ",j),
                      dotsize=0,binwidth=0.2,Xangle=45, ylo = -0.1, yhi = 0.5 ,
                      ordered_x = names(ind_df))      
   print (p)
 }
 dev.off()                                     
 
 saveRDS(exprSet_Comb, file= "exprSet_Comb.rds")
 rm(exprSet_Comb)
 
 #=============================
 
 # MatchMixeR
 
 #=============================
 
 exprSet_MM<-list()
 
 for (i in names(exprSet_noCPN)){
   
   Babatch <- as.factor ( CombSampleAnnot[[i]][,BatchCol] )
   contrasts(Babatch) <- contr.sum(levels(Babatch))
   Babatch <- model.matrix(~Babatch)[, -1, drop = FALSE]
   arrayMat <- exprSet_noCPN[[i]][,Babatch==1]
   RseqMat <- exprSet_noCPN[[i]][,Babatch==-1]
   
   exprSet_MM[[i]]<-MM(Xmat=RseqMat,Ymat=arrayMat)
   exprSet_MM[[i]]<-cbind(arrayMat,as.data.frame(exprSet_MM[[i]]$Yhat))
   
   print(paste0(i, ": done"))
 }
 
 
 #  violin plot
 #===============                          
 
 exprSet_toplot <- exprSet_MM
 name_toplot <-"MatchMixer"
 n <-ncol(exprSet_toplot[[1]])
 
 pdf(paste0("after_CPN_violinPlot_4genes_",name_toplot,".pdf"),height = 4, width = 6)       
 
 for (j in random_genes){
   
   gene_df <-do.call(rbind, (lapply(exprSet_toplot, function(x) x[j,((n/2)+1):n])))   #because "NOT array samples" are at the second part of the dataframe
   gene_df<- as.data.frame(t(gene_df))
   gene_df$array <-as.numeric(  exprSet_toplot[[1]][j,1:(n/2)]  )  # because array samples are at the first part of the dataframe
   p<-violinPlot(gene_df, mainTitle=paste0(name_toplot," - gene = ",j),
                 dotsize=0,binwidth=0.2,Xangle=45, ordered_x = names(gene_df))      
   print (p)
   
 }
 dev.off()
 
 pdf(paste0("after_CPN_violinPlot_4sples_",name_toplot,".pdf"),height = 4, width = 6)   
 
 for (j in random_sples){
   denst <-lapply( lapply( exprSet_toplot, function(y)  y=y[, ((n/2)+1):n]) , 
                   function(x) density(x[,j])  )  
   
   ind_df <-do.call(rbind, (lapply(denst, function(dens) dens$y)  ) )
   ind_df<- as.data.frame(t(ind_df))
   densarray <- density(exprSet_array[,paste0(j)])
   ind_df$array <- densarray$y
   p<-violinPlot_Ylim(ind_df, mainTitle=paste0(name_toplot, " - pat = ",j),
                      dotsize=0,binwidth=0.2,Xangle=45, ylo = -0.1, yhi = 0.5 ,
                      ordered_x = names(ind_df))      
   print (p)
 }
 dev.off()                          
 
 saveRDS(exprSet_MM, file="exprSet_MM.rds")
 rm(exprSet_MM)
 
 
 #=============================
 
 #  GQ
 
 #=============================
 
 
 exprSet_GQ<-list()
 
 for (i in names(exprSet_noCPN)){
   
   Babatch <- as.factor ( CombSampleAnnot[[i]][,BatchCol] )
   contrasts(Babatch) <- contr.sum(levels(Babatch))
   Babatch <- model.matrix(~Babatch)[, -1, drop = FALSE]
   arrayMat <- exprSet_noCPN[[i]][,Babatch==1]
   RseqMat <- exprSet_noCPN[[i]][,Babatch==-1]
   
   exprSet_GQ[[i]]<- gq(platform1.data=arrayMat, platform2.data=RseqMat)
   exprSet_GQ[[i]]<-cbind(arrayMat,exprSet_GQ[[i]]$y)
   
   print(paste0(i, ": done"))
 }
 
 
 #  violin plot
 #===============                          
 
 exprSet_toplot <- exprSet_GQ
 name_toplot <-"GQ"
 n <-ncol(exprSet_toplot[[1]])
 
 pdf(paste0("after_CPN_violinPlot_4genes_",name_toplot,".pdf"),height = 4, width = 6)       
 
 for (j in random_genes){
   
   gene_df <-do.call(rbind, (lapply(exprSet_toplot, function(x) x[j,((n/2)+1):n])))   #because "NOT array samples" are at the second part of the dataframe
   gene_df<- as.data.frame(t(gene_df))
   gene_df$array <-as.numeric(  exprSet_toplot[[1]][j,1:(n/2)]  )  # because array samples are at the first part of the dataframe
   p<-violinPlot(gene_df, mainTitle=paste0(name_toplot," - gene = ",j),
                 dotsize=0,binwidth=0.2,Xangle=45, ordered_x = names(gene_df))      
   print (p)
   
 }
 dev.off()
 
 pdf(paste0("after_CPN_violinPlot_4sples_",name_toplot,".pdf"),height = 4, width = 6)   
 
 for (j in random_sples){
   denst <-lapply( lapply( exprSet_toplot, function(y)  y=y[, ((n/2)+1):n]) , 
                   function(x) density(x[,j])  )  
   
   ind_df <-do.call(rbind, (lapply(denst, function(dens) dens$y)  ) )
   ind_df<- as.data.frame(t(ind_df))
   densarray <- density(exprSet_array[,paste0(j)])
   ind_df$array <- densarray$y
   p<-violinPlot_Ylim(ind_df, mainTitle=paste0(name_toplot, " - pat = ",j),
                      dotsize=0,binwidth=0.2,Xangle=45, ylo = -0.1, yhi = 0.5 ,
                      ordered_x = names(ind_df))      
   print (p)
 }
 dev.off()
                     
saveRDS(exprSet_GQ, file="exprSet_GQ.rds")
rm(exprSet_GQ)
                     
                     
#==================================

                     
#  zscoring
                     
                     
#=====================================
      
                               
      
exprSet_ZSc<-list()
                     
for (i in names(exprSet_noCPN)){
exprSet_ZSc[[i]] <-BatchGeneTrans(dat=exprSet_noCPN[[i]],batch=CombSampleAnnot[[i]][,BatchCol])
                       
print(paste0(i, ": done"))
}
                     
                     
#  violin plot
#===============                          
                     
exprSet_toplot <- exprSet_ZSc
name_toplot <-"ZSc"
n <-ncol(exprSet_toplot[[1]])
                     
pdf(paste0("after_CPN_violinPlot_4genes_",name_toplot,".pdf"),height = 4, width = 6)       

for (j in random_genes){

gene_df <-do.call(rbind, (lapply(exprSet_toplot, function(x) x[j,((n/2)+1):n])))   #because "NOT array samples" are at the second part of the dataframe
gene_df<- as.data.frame(t(gene_df))
gene_df$array <-as.numeric(  exprSet_toplot[[1]][j,1:(n/2)]  )  # because array samples are at the first part of the dataframe
p<-violinPlot(gene_df, mainTitle=paste0(name_toplot," - gene = ",j),
         dotsize=0,binwidth=0.2,Xangle=45, ordered_x = names(gene_df))      
print (p)

}
dev.off()

pdf(paste0("after_CPN_violinPlot_4sples_",name_toplot,".pdf"),height = 4, width = 6)   

for (j in random_sples){
denst <-lapply( lapply( exprSet_toplot, function(y)  y=y[, ((n/2)+1):n]) , 
           function(x) density(x[,j])  )  
                       
ind_df <-do.call(rbind, (lapply(denst, function(dens) dens$y)  ) )
ind_df<- as.data.frame(t(ind_df))
densarray <- density(exprSet_array[,paste0(j)])
ind_df$array <- densarray$y
p<-violinPlot_Ylim(ind_df, mainTitle=paste0(name_toplot, " - pat = ",j),
                      dotsize=0,binwidth=0.2,Xangle=45, ylo = -0.1, yhi = 0.5 ,
                        ordered_x = names(ind_df))      
print (p)
}
dev.off()                                   
                     
saveRDS(exprSet_ZSc, file="exprSet_ZSc.rds")
                     
                     
# the end, see step 3 for data visualisation and performance calculation
                                
                                
