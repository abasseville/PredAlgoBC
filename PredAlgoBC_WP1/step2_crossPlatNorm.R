#!/usr/bin/Rscript

#===============================================================================================
#=========================================================================================

#       STEP 2 : Cross-platform normalization

#=========================================================================================
#===============================================================================================

# matchmixer function was available at https://github.com/dy16b/Cross-Platform-Normalization (publication: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7868049/)

# load dataset
#==============

setwd("tcga)
exprSet <- readRDS("exprSet.rds")
exprSet <- exprSet[-1]   # remove "count" dataset from the list (1st one)
sampleAnnot <- readRDS("sampleAnnot.rds")

sampleAnnot_All <- list()
for (i in names(exprSet) ){   
  sampleAnnot_All[[i]] <- sampleAnnot
  sampleAnnot_All[[i]]$batch <- i
}

BatchCol<-"batch" 


# load library and function
#===========================

library(preprocessCore)
library(TDM)
library(limma)   # RBE
library(sva)   # Combat
library(CONOR)   # xpn


# custom function
#===================

quantileNormalizeByFeature2 <-function (matrix_to_normalize, target_distribution_matrix) 
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
  #This function was provided by Xiao-Qin Xia, one of the authors of webarraydb.
  # modified MRS
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




#==========================================================

#     Merging each RNAseq table with affy

#=========================================================


exprSet_affy <- exprSet[["affy"]]
exprSet_Rseq <- exprSet[-(length(exprSet)]
sampleAnnot_affy<- sampleAnnot[["affy"]]
sampleAnnot_Rseq <- sampleAnnot[-(length(exprSet)]


CombSampleAnnot <-list()
for (i in names(exprSet_Rseq) ){
CombSampleAnnot[[i]] <- rbind(sampleAnnot_affy,sampleAnnot_Rseq[[i]])
CombSampleAnnot[[i]]$subtype <- "Lum"
CombSampleAnnot[[i]]$subtype[CombSampleAnnot[[i]]$er.status.by.ihc == "Negative" & CombSampleAnnot[[i]]$pr.status.by.ihc == "Negative" & CombSampleAnnot[[i]]$her2.status.by.ihc == "Negative"]<-"TripleNeg"
CombSampleAnnot[[i]]$subtype[CombSampleAnnot[[i]]$her2.status.by.ihc == "Positive"]<-"HER2"
}


CombexprSet<-list()
for (i in names(sampleAnnot_Rseq) ){
  print(paste0("nrow before merge:",nrow(sampleAnnot_Rseq[[i]])))
  print(paste0("ncol before merge:", ncol(sampleAnnot_Rseq[[i]]), " - expected after: ", (ncol(sampleAnnot_Rseq[[i]])*2)))
  CombexprSet[[i]] <- merge(exprSet_affy, sampleAnnot_Rseq[[i]],
                        by= "row.names")
  row.names(CombexprSet[[i]]) <-CombexprSet[[i]][,1]; CombexprSet[[i]] <-CombexprSet[[i]][,-1]
  print(paste0("nrow after merge:",nrow(CombexprSet[[i]])))
  print(paste0("ncol after merge:",ncol(CombexprSet[[i]])))
}

saveRDS(CombsampleAnnot,file="noCPN_sampleAnnot.rds")
saveRDS(CombexprSet,file="noCPN_exprSet.rds")


#========================

# TDM

#==================

tdm <-list()
for  (i in names(exprSet)[-(length(exprSet)] ) {      # no transformation for the last of the list:  affy
  tdm [[i]]<-TDM::tdm_transform(ref_data = data.table(cbind(gene=rownames(exprSet[["affy"]]), 
                                                            exprSet[["affy"]])),
                                target_data = data.table(cbind(gene=rownames(exprSet[[i]]), 
                                                               exprSet[[i]])),
                                log_target = TRUE, tdm_factor = 1)
                                
  # transform to numeric (not numric with datatable)
  tdm[[i]] <- as.data.frame(tdm[[i]])
  row.names.tdm<-tdm [[i]]$gene
  tdm [[i]] <- as.data.frame(sapply(tdm [[i]][,-1], as.numeric) )
  row.names(tdm[[i]]) <- row.names.tdm
  print(paste0(names(exprSet)[[i]], ": done"))
}

saveRDS(tdm, file= "exprSet_TDM.rds")
rm(tdm)

#===================

#     FSQN

#===================
 

fsqn <-list()
for  (i in names(exprSet)[-(length(exprSet)] ) {
  fsqn[[i]] <- quantileNormalizeByFeature2(as.matrix(exprSet[[i]]),   #test
                                           as.matrix(exprSet[["affy"]]))   #target
  print(paste0("exprSet ", i, ": done"))
}

saveRDS(fsqn, file= "exprSet_FSQN.rds")
rm(fsqn)

#===================

#     RBE  

#===================


exprSet_RBE <-list()
for (i in names(CombexprSet)) {
  exprSet_RBE[[i]]<-limma::removeBatchEffect(CombexprSet[[i]],batch = CombsampleAnnot[[i]][,batchCol])
  exprSet_RBE[[i]]<-exprSet_RBE[[i]]-min(exprSet_RBE[[i]])
  print(paste0(i, " done"))
}

saveRDS(exprSet_RBE, file= "exprSet_RBE.rds")
rm(exprSet_RBE)

#=============================

# Combat   #sva dependant

#=============================

exprSet_Comb <-list()
for (i in names(CombexprSet)){
  exprSet_Comb[[i]]<-sva::ComBat(as.matrix(CombexprSet[[i]]),batch = CombsampleAnnot[[i]][,batchCol],mod=NULL, par.prior = TRUE, prior.plots = FALSE)
  print(paste0(i, ": done"))
}

saveRDS(exprSet_Comb, file= "exprSet_Comb.rds")
rm(exprSet_Comb)

#=============================

# MatchMixeR

#=============================

exprSet_MM<-list()

for (i in names(CombexprSet)){
  
  Babatch <- as.factor ( CombsampleAnnot[[i]][,batchCol] )
  contrasts(Babatch) <- contr.sum(levels(Babatch))
  Babatch <- model.matrix(~Babatch)[, -1, drop = FALSE]
  affyMat <- CombexprSet[[i]][,Babatch==1]
  RseqMat <- CombexprSet[[i]][,Babatch==-1]
  
  exprSet_MM[[i]]<-MM(Xmat=RseqMat,Ymat=affyMat)
  exprSet_MM[[i]]<-cbind(affyMat,as.data.frame(exprSet_MM[[i]]$Yhat))
  
  print(paste0(i, ": done"))
}

save(exprSet_MM, file="exprSet_MM.rda")
rm(exprSet_MM)

#=============================

#  GQ

#=============================


exprSet_GQ<-list()

for (i in names(CombexprSet)){
  
  Babatch <- as.factor ( CombsampleAnnot[[i]][,batchCol] )
  contrasts(Babatch) <- contr.sum(levels(Babatch))
  Babatch <- model.matrix(~Babatch)[, -1, drop = FALSE]
  affyMat <- CombexprSet[[i]][,Babatch==1]
  RseqMat <- CombexprSet[[i]][,Babatch==-1]
  
  exprSet_GQ[[i]]<- gq(platform1.data=affyMat, platform2.data=RseqMat)
  exprSet_GQ[[i]]<-cbind(affyMat,exprSet_GQ[[i]]$y)
  
  print(paste0(i, ": done"))
}


saveRDS(exprSet_GQ, file="exprSet_GQ.rds")
rm(exprSet_GQ)

#=============================

#  XPN

#=============================

exprSet_XPN<-list()

for (i in names(CombexprSet)){
  Babatch <- as.factor ( CombsampleAnnot[[i]][,batchCol] )
  contrasts(Babatch) <- contr.sum(levels(Babatch))
  Babatch <- model.matrix(~Babatch)[, -1, drop = FALSE]
  affyMat <- CombexprSet[[i]][,Babatch==1]
  RseqMat <- CombexprSet[[i]][,Babatch==-1]
  
  exprSet_XPN[[i]]<-xpn(affyMat,RseqMat)
  exprSet_XPN[[i]] <-cbind(exprSet_XPN[[i]]$x,exprSet_XPN[[i]]$y)
  
  print(paste0(i, ": done"))
}

saveRDS(exprSet_XPN, file="exprSet_XPN.rds")


# the end



