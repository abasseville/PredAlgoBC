#!/usr/bin/Rscript

setwd("tcga")


library(SummarizedExperiment)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(vsn)
library(hexbin)
library(DESeq2)
library(edgeR)

# custom functions
#=================

  test_match_order <- function(x,y) {
    if (isTRUE(all.equal(x,y))) print('Perfect match in SAME order')
    if (!isTRUE(all.equal(x,y)) && isTRUE(all.equal(sort(x),sort(y)))) print('Perfect match in WRONG order')
    if (!isTRUE(all.equal(x,y)) && !isTRUE(all.equal(sort(x),sort(y)))) {
      if( length(intersect(x,y)) ==0  )  {
        print('No matched element at all') 
        }else{
        print(paste0 ("not in vector1: ", setdiff(y,x)) )
        print(paste0 ("not in vector2: ", setdiff(x,y)) )
        } # else
    }
  }

Vst <- function(countdata){
  condition <- factor(rep("Tumour", ncol(countdata)))
  countdata <- DESeq2::DESeqDataSetFromMatrix(countdata, data.frame(condition), ~ 1, tidy = FALSE )
  Vstdata <- DESeq2::vst( countdata,blind = TRUE ,nsub=1000, fitType = "parametric")   #en fait, il existe aussi une fct vst encore plus rapide
  return(SummarizedExperiment::assay(Vstdata))
}

VstBis <- function(countdata){
  condition <- factor(rep("Tumour", ncol(countdata)))
  countdata <- DESeq2::DESeqDataSetFromMatrix(countdata, data.frame(condition), ~ 1, tidy = FALSE )
  Vstdata <- DESeq2::vst( countdata,blind = TRUE ,nsub=1000, fitType = "mean")   #en fait, il existe aussi une fct vst encore plus rapide
  return(SummarizedExperiment::assay(Vstdata))
}


violinPlot = function (df, mainTitle,dotsize,binwidth,Xangle =0){
  dataf <- gather(df ,key="Data", value="Val")
  ggplot(dataf, aes(x=Data, y=Val, fill=Data)) +
    theme(axis.text.x = element_text(angle = Xangle, hjust = 1))+
    ggtitle(mainTitle)+
    scale_x_discrete(limits=names(df))+    #to avoid ggplot to reorder alph automatiqualy  
    geom_violin(trim = FALSE)+
    geom_dotplot(binaxis='y', stackdir='center',dotsize=dotsize,fill = "black",binwidth = binwidth)
}

violinPlot_Ylim = function (df, mainTitle,dotsize,binwidth,Xangle =0,ylo,yhi){
  dataf <- gather(df ,key="Data", value="Val")
  ggplot(dataf, aes(x=Data, y=Val, fill=Data)) +
    theme(axis.text.x = element_text(angle = Xangle, hjust = 1))+
    ggtitle(mainTitle)+ ylim(ylo,yhi)+
    scale_x_discrete(limits=names(df))+    #to avoid ggplot to reorder alph automatiqualy  
    geom_violin(trim = FALSE)+
    geom_dotplot(binaxis='y', stackdir='center',dotsize=dotsize,fill = "black",binwidth = binwidth)
}




#=============================================================================================================
#=============================================================================================================

#  load TCGA data (previouly downloaded using getTCGA , see https://rdrr.io/cran/TCGA2STAT/man/getTCGA.html

#=============================================================================================================
#================================================================================================

# micorarray (Agilent)
#========================

exprSet_affy <- read.table("TCGA.microarray.txt",row.names = 1)  
exprSet_affy <- exprSet_affy[,-1]

# RNAseq
#==============

exprSet_Rsq_count <- read.csv("TCGA.tumor_FeatureCounts.1119_Breast.csv",row.names = 1,sep = ";")  
exprSet_Rsq_fpkm <- read.table("GSM1536837_TCGA.FPKM.breast_1119tum.txt",row.names = 1) 
exprSet_Rsq_tpm <- read.table("GSM1536837_TCGA.TPM.breast_1119tum.txt",row.names = 1) 
 
sampleAnnot <-  read.csv("TCGA_clinical.1119_Breast et NAT.csv",sep = ";",row.names = 1)
sampleAnnot <- sampleAnnot[-1,]
rownames(sampleAnnot) <-gsub("-",".",rownames(sampleAnnot))


# keep patients in common

Rseq_com_pat <- intersect(names(exprSet_Rsq_count), intersect(names(exprSet_Rsq_fpkm), names(exprSet_Rsq_tpm)))
com_pat <- intersect(Rseq_com_pat,names(exprSet_affy))
com_pat2 <- intersect(com_pat,rownames(sampleAnnot))  # verif it's ok wiht sampleAnnot


com_cols <- names(exprSet_affy) %in% com_pat
exprSet_affy <- exprSet_affy [,com_cols]
com_cols <- names(exprSet_Rsq_count) %in% com_pat
exprSet_Rsq_count <- exprSet_Rsq_count [,com_cols]
com_cols <- names(exprSet_Rsq_fpkm) %in% com_pat
exprSet_Rsq_fpkm <- exprSet_Rsq_fpkm [,com_cols]
com_cols <- names(exprSet_Rsq_tpm) %in% com_pat
exprSet_Rsq_tpm <- exprSet_Rsq_tpm [,com_cols]

com_rows <- rownames(sampleAnnot) %in% com_pat
sampleAnnot <- sampleAnnot [com_rows,]

test_match_order(rownames(sampleAnnot), names(exprSet_affy))


#=============================================================================================================
#=============================================================================================================

#  RNA transformation

#=============================================================================================================
#================================================================================================


#===========================================
#  count,  log2 count, FPKM and TDM
#================================

#  genes with too much 0 count values

zero_val<-(rowSums(exprSet_Rsq_count<5)>round(0.995*(ncol(exprSet_Rsq_count)), digits = 0) )
table(zero_val)
gene_out_rsq <- names(zero_val[zero_val=="TRUE"])
print(paste0("length gene_out_rsq: ",length(gene_out_rsq))


# remove genes and apply log2

lowexp_rows <- rownames(exprSet_Rsq_count) %in% gene_out_rsq
exprSet_Rsq_count_nolo <-exprSet_Rsq_count[!lowexp_rows,]

exprSet_Rsq_log <- log2(exprSet_Rsq_count+1)
exprSet_Rsq_log_nolo <- exprSet_Rsq_log[!lowexp_rows,]

lowexp_rows <- rownames(exprSet_Rsq_fpkm) %in% gene_out_rsq
exprSet_Rsq_fpkm_nolo <-exprSet_Rsq_fpkm[!lowexp_rows,]
exprSet_Rsq_fpkm_nolo <- log2(exprSet_Rsq_fpkm_nolo+1)

lowexp_rows <- rownames(exprSet_Rsq_tpm) %in% gene_out_rsq
exprSet_Rsq_tpm_nolo <-exprSet_Rsq_tpm[!lowexp_rows,]
exprSet_Rsq_tpm_nolo <- log2(exprSet_Rsq_tpm_nolo+1)


#==========================
#  VST, RLE and TMM
#=======================



exprSet_Rsq_vst_NoLo <- Vst(exprSet_Rsq_count[!lowexp_rows,])
exprSet_Rsq_vst_NoLo <-as.data.frame(exprSet_Rsq_vst_NoLo)

exprSet_Rsq_vst2_NoLo <- VstBis(exprSet_Rsq_count[!lowexp_rows,])
exprSet_Rsq_vst2_NoLo <-as.data.frame(exprSet_Rsq_vst2_NoLo)


dge_counts <- edgeR::DGEList(counts=exprSet_Rsq_count[!lowexp_rows,], group=factor(rep("Tumour", ncol(exprSet_Rsq_count[!lowexp_rows,]))))
normfact <-edgeR::calcNormFactors(dge_counts,method="RLE")
exprSet_Rsq_rle_NoLo <- edgeR::cpm(normfact, log=TRUE,prior.count=1)

normfact <-edgeR::calcNormFactors(dge_counts,method="TMM")
exprSet_Rsq_tmm_NoLo<- edgeR::cpm(normfact, log=TRUE,prior.count=1)


#==========================
# adapted affy and save all
#==========================

lowexp_rows <- rownames(exprSet_affy) %in% gene_out_rsq
exprSet_affy <-exprSet_affy[!lowexp_rows,]
save(exprSet_affy,file= "exprSet_affy.rda")


exprSet <- list(exprSet_Rsq2_count_NoLo,exprSet_Rsq2_log_NoLo,exprSet_Rsq2_vst_NoLo,exprSet_Rsq2_vst2_NoLo,
                exprSet_Rsq2_fpkm_NoLo, exprSet_Rsq2_tpm_NoLo, exprSet_Rsq2_rle_NoLo, exprSet_Rsq2_tmm_NoLo, exprSet_affy)

names(exprSet) <-c("cts", "log","vst","vst2","fpkm","tpm","rle", "tmm","affy")

saveRDS(exprSet, "exprSet.rds")
saveRDS(sampleAnnot,file= "sampleAnnot.rds")



#======================================================================================================
#======================================================================================================
#
#                          plot after RNAseq transformation
#
#======================================================================================================
#======================================================================================================


# SD plot 
#==========

SDplot2 <- list()
for ( i in 1:length(exprSet)){
  if(names(exprSet)[[i]] %in% "cts"){
    SDplot <-meanSdPlot(as.matrix(exprSet[[i]]), plot=FALSE ,bins  = 200, ylab = paste0("SD_", names(exprSet)[[i]]) ) 
    SDplot2[[i]] <-SDplot$gg + ggtitle(names(exprSet)[[i]]) + scale_y_continuous(limits = c(0,150000))
    print(paste(names(exprSet)[[i]] , "in large scale"))
  }else{
    SDplot <-meanSdPlot(as.matrix(exprSet[[i]]),plot=FALSE ,bins  = 200) 
    SDplot2[[i]] <-SDplot$gg + ggtitle(names(exprSet)[[i]])+scale_y_continuous(limits = c(0,4))
    print(paste(names(exprSet)[[i]] , "in small scale"))
  }
}


pdf("after_transf_SDPlot.pdf",height = 4, width = 4)
for ( i in 1:length(exprSet)){    
print(SDplot2[[i]])
}
dev.off()

# Violin plot fo gene density
#==============================


pdf("after_transf_violinPlot_4genes.pdf",height = 4, width = 6)   
random_genes <-sample(row.names(as.data.frame(exprSet_nocts[["affy"]][1])),size = 4,replace =F)

lgExp <-(length(exprSet_nocts))-1

for (j in random_genes){
  gene_df <-do.call(rbind, (lapply(exprSet_nocts[1:lgExp], function(x) x[j,])))
  gene_df<- as.data.frame(t(gene_df))
  gene_df$affy <-as.numeric(exprSet_nocts[["affy"]][j,])
  p<-violinPlot(gene_df, mainTitle=paste0("gene = ",j),
                dotsize=0,binwidth=0.2,Xangle=45)      #binwidth
  print (p)
}

dev.off()


# Violin plot fo sample density  
#=================================


IndDsty_med <-list()
for (i in 1:length(exprSet_nocts)){
  IndDsty <- apply(exprSet_nocts[[i]], 2, density)
  IndDsty2 <-do.call(rbind, (lapply(IndDsty, function(dens) dens$y)  ) )
  IndDsty_med[[i]] <- apply(IndDsty2, 2, median)
}

ind_df <- data.frame(matrix(unlist(IndDsty_med), nrow=length(IndDsty_med), byrow=T),row.names = names(exprSet_nocts) )
ind_df<- as.data.frame(t(ind_df))


p_all1<-violinPlot(ind_df, mainTitle="median all samples",
                   dotsize=0,binwidth=0.2,Xangle=45)      #binwidth

p_all2<-violinPlot_Ylim(ind_df, mainTitle="median all samples - adjusted y-axe",
                        dotsize=0,binwidth=0.2,Xangle=45, ylo = -0.1, yhi = 0.5) 

pdf("after_transf_violinPlot_rand4sample_med.pdf",height = 4, width = 6)

print (p_all1)
print(p_all2)

random_ind <-sample(names(exprSet_nocts[["log"]]),size = 4,replace =F)
for (j in random_ind){
  prout <- lapply(exprSet_nocts[1:lgExp], function(x) density(x[,j])  )
  ind_df <-do.call(rbind, (lapply(prout, function(dens) dens$y)  ) )
  ind_df<- as.data.frame(t(ind_df))
  densAffy <- density(exprSet_nocts[["affy"]][,paste0(j,"_affy")])
  ind_df$affy <- densAffy$y
  p<-violinPlot_Ylim(ind_df, mainTitle=paste0("pat = ",j),
                     dotsize=0,binwidth=0.2,Xangle=45, ylo = -0.1, yhi = 0.5)      
  print (p)
}


dev.off()
