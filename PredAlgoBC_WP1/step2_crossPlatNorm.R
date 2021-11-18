#!/usr/bin/Rscript

#===============================================================================================
#=========================================================================================

STEP 2 : Cross-platform normalization

#=========================================================================================
#===============================================================================================


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


# load library and function
#===========================

library(preprocessCore)
library(TDM)

# custom function

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


#========================

# TDM

#==================

tdm <-list()
for  (i in names(exprSet)[-(length(exprSet)] ) {      # no transformation for the last of the list:  affy
  tdm [[i]]<-TDM::tdm_transform(ref_data = data.table(cbind(gene=rownames(exprSet[["affy"]]), 
                                                            exprSet[["affy"]])),
                                target_data = data.table(cbind(gene=rownames(exprSet_nocts[[i]]), 
                                                               exprSet_nocts[[i]])),
                                log_target = TRUE, tdm_factor = 1)
                                
  # transform to numerci en num (not the case with datatable)
  tdm[[i]] <- as.data.frame(tdm[[i]])
  row.names.tdm<-tdm [[i]]$gene
  tdm [[i]] <- as.data.frame(sapply(tdm [[i]][,-1], as.numeric) )
  row.names(tdm[[i]]) <- row.names.tdm
  print(paste0(names(exprSet_nocts)[[i]], ": done"))
}

#===================

#     FSQN

#===================
 

fsqn <-list()
for  (i in names(exprSet)[-(length(exprSet)] ) {
  fsqn[[i]] <- quantileNormalizeByFeature2(as.matrix(exprSet[[i]]),   #test
                                                as.matrix(exprSet[["affy"]]))   #target
  
  print(paste0("exprSet ", i, ": done"))
}



save(fsqn, file= "fsqn.rda")


#===================

#     FSQN

#===================


