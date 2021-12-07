#!/usr/bin/Rscript



dir.create("tcga",showWarnings = F)
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
library(TCGA2STAT)   # not updated in CRAN

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

#  load TCGA data (previouly downloaded using getTCGA , see https://rdrr.io/cran/TCGA2STAT/man/getTCGA.html)
# 

#=============================================================================================================
#================================================================================================


# important information: TCGA2STAT not updated on CRAN anymore, so impossible to dowload it   => used datasets available upon request
# available on GSM1536837 but very heavy file (all cancers included)

foldersource <- "path_for_folder_with_needed_files"

# sampleAnnot
#========================

sampleAnnot <-  read.csv(paste0(foldersource,"TCGA_clinical.1119_Breast et NAT.csv"),sep = ";",row.names = 1)  # from  TCGA2STAT
sampleAnnot[1:5,1:5]
sampleAnnot <- sampleAnnot[-1,] # remove duplicated first row
rownames(sampleAnnot) <-gsub("-",".",rownames(sampleAnnot))  # to be the same thean transcr.data

# micorarray (Agilent)
#========================

exprSet_array <- read.table(paste0(foldersource,"TCGA.microarray.txt"),row.names = 1)  
exprSet_array[1:5,1:5]
exprSet_array <- exprSet_array[,-1] # remove duplicated rownames column


# RNAseq  
#==============================================

exprSet_Rsq_count<- read.table(paste0(foldersource,"GSM1536837_01_27_15_TCGA_20.Illumina.tumor_Rsubread_FeatureCounts.txt"),row.names = 1)   #(download from GSM1536837 on GEO)
com_pat <- intersect(names(exprSet_Rsq_count), rownames(sampleAnnot))  # remove  other than breast cancers

exprSet_Rsq_count <- exprSet_Rsq_count [,names(exprSet_Rsq_count) %in% com_pat]  # selection
exprSet_Rsq_count<- exprSet_Rsq_count[,com_pat]   # order
saveRDS(exprSet_Rsq_count,"exprSet_Rsq_count.rds")    


exprSet_Rsq_fpkm<- read.table(paste0(foldersource,"GSM1536837_TCGA.FPKM.breast_1119tum.txt"),row.names = 1)   # from  TCGA2STAT
saveRDS(exprSet_Rsq_fpkm,"exprSet_Rsq_fpkm.rds")
exprSet_Rsq_fpkm[1:5,1:5]

exprSet_Rsq_tpm <- read.table(paste0(foldersource,"GSM1536837_TCGA.TPM.breast_1119tum3.txt"), row.names = 1) # from  TCGA2STAT
saveRDS(exprSet_Rsq_tpm,"exprSet_Rsq_tpm.rds")
exprSet_Rsq_tpm[1:5,1:5]



# keep patients in common

Rseq_com_pat <- intersect(names(exprSet_Rsq_count), intersect(names(exprSet_Rsq_fpkm), names(exprSet_Rsq_tpm)))
com_pat <- intersect(Rseq_com_pat,names(exprSet_array))
com_pat2 <- intersect(com_pat,rownames(sampleAnnot))  # verif it's ok with sampleAnnot


exprSet_array <- exprSet_array [,names(exprSet_array) %in% com_pat]  # selection
exprSet_array<- exprSet_array[,com_pat]                             # order

exprSet_Rsq_count <- exprSet_Rsq_count [,names(exprSet_Rsq_count) %in% com_pat]
exprSet_Rsq_count<- exprSet_Rsq_count[,com_pat]

exprSet_Rsq_fpkm <- exprSet_Rsq_fpkm [,names(exprSet_Rsq_fpkm) %in% com_pat]
exprSet_Rsq_fpkm<- exprSet_Rsq_fpkm[,com_pat]

exprSet_Rsq_tpm <- exprSet_Rsq_tpm [,names(exprSet_Rsq_tpm) %in% com_pat]
exprSet_Rsq_tpm<- exprSet_Rsq_tpm[,com_pat]

sampleAnnot <- sampleAnnot [rownames(sampleAnnot) %in% com_pat,]
sampleAnnot<- sampleAnnot[com_pat,]

test_match_order(rownames(sampleAnnot), names(exprSet_array))
test_match_order(rownames(sampleAnnot), names(exprSet_Rsq_count))
test_match_order(names(exprSet_Rsq_count), names(exprSet_Rsq_fpkm))
test_match_order(names(exprSet_Rsq_fpkm), names(exprSet_Rsq_tpm))



# keep gene in common

Rseq_com_gene <- intersect(rownames(exprSet_Rsq_count), intersect(rownames(exprSet_Rsq_fpkm), rownames(exprSet_Rsq_tpm)))
com_gene <- intersect(Rseq_com_gene,rownames(exprSet_array))

exprSet_array <- exprSet_array [rownames(exprSet_array) %in% com_gene,]  # selection
exprSet_array<- exprSet_array[com_gene,]                             # order

exprSet_Rsq_count <- exprSet_Rsq_count [rownames(exprSet_Rsq_count) %in% com_gene,]
exprSet_Rsq_count<- exprSet_Rsq_count[com_gene,]

exprSet_Rsq_fpkm <- exprSet_Rsq_fpkm [rownames(exprSet_Rsq_fpkm) %in% com_gene,]
exprSet_Rsq_fpkm<- exprSet_Rsq_fpkm[com_gene,]

exprSet_Rsq_tpm <- exprSet_Rsq_tpm [rownames(exprSet_Rsq_tpm) %in% com_gene,]
exprSet_Rsq_tpm<- exprSet_Rsq_tpm[com_gene,]


test_match_order(rownames(exprSet_array), rownames(exprSet_Rsq_count))
test_match_order(rownames(exprSet_Rsq_count), rownames(exprSet_Rsq_fpkm))
test_match_order(rownames(exprSet_Rsq_fpkm), rownames(exprSet_Rsq_tpm))

saveRDS(exprSet_array,file= "exprSet_array.rds")


#=============================================================================================================
#=============================================================================================================

#  RNA transformation

#=============================================================================================================
#================================================================================================


#===========================================
#  count,  log2 count, FPKM and TDM
#================================

#  genes with too much 0 count values

low_val<-(rowSums(exprSet_Rsq_count<5)>round(0.995*(ncol(exprSet_Rsq_count)), digits = 0) )
table(low_val)

gene_out_rsq <- names(low_val[low_val=="TRUE"])
print(paste0("length gene_out_rsq: ",length(gene_out_rsq)))
      
      
# remove genes and apply log2
      
lowexp_rows <- rownames(exprSet_Rsq_count) %in% gene_out_rsq
exprSet_Rsq_count <-exprSet_Rsq_count[!lowexp_rows,]
      
exprSet_Rsq_log <- log2(exprSet_Rsq_count+1)
exprSet_Rsq_log <- exprSet_Rsq_log[!lowexp_rows,]
      

exprSet_Rsq_fpkm <-exprSet_Rsq_fpkm[!lowexp_rows,]
exprSet_Rsq_fpkm <- log2(exprSet_Rsq_fpkm+1)
      
exprSet_Rsq_tpm <-exprSet_Rsq_tpm[!lowexp_rows,]
exprSet_Rsq_tpm <- log2(exprSet_Rsq_tpm+1)

exprSet_array<-exprSet_array[!lowexp_rows,]


      
#==========================
 #  VST, RLE and TMM
#=======================
      
      
exprSet_Rsq_vst <- Vst(exprSet_Rsq_count)
exprSet_Rsq_vst <-as.data.frame(exprSet_Rsq_vst)
      
exprSet_Rsq_vst2 <- VstBis(exprSet_Rsq_count)
exprSet_Rsq_vst2 <-as.data.frame(exprSet_Rsq_vst2)
      
      
dge_counts <- edgeR::DGEList(counts=exprSet_Rsq_count, group=factor(rep("Tumour", ncol(exprSet_Rsq_count))))
normfact <-edgeR::calcNormFactors(dge_counts,method="RLE")
exprSet_Rsq_rle <- edgeR::cpm(normfact, log=TRUE,prior.count=1)
      
normfact <-edgeR::calcNormFactors(dge_counts,method="TMM")
exprSet_Rsq_tmm<- edgeR::cpm(normfact, log=TRUE,prior.count=1)
      
      
#==========================
# adapted array and save all
#==========================
      
exprSet <- list(exprSet_Rsq_count,exprSet_Rsq_log,exprSet_Rsq_vst,exprSet_Rsq_vst2,
                      exprSet_Rsq_fpkm, exprSet_Rsq_tpm, exprSet_Rsq_rle, exprSet_Rsq_tmm, exprSet_array)
      
names(exprSet) <-c("cts", "log","vst","vst2","fpkm","tpm","rle", "tmm","array")
      
saveRDS(exprSet, "exprSet.rds")
saveRDS(sampleAnnot,file= "sampleAnnot.rds")
      
      
      
#======================================================================================================
#======================================================================================================
#
#                          plot after RNAseq transformation
#
#======================================================================================================
#======================================================================================================
      
      
# SD plot (with adapted scale)
#====================================
      
SDplot2 <- list()
for ( i in 1:length(exprSet)){
        if(names(exprSet)[[i]] %in% "cts"){
          SDplot <-vsn::meanSdPlot(as.matrix(exprSet[[i]]), plot=FALSE ,bins  = 200, ylab = paste0("SD_", names(exprSet)[[i]]) ) 
          SDplot2[[i]] <-SDplot$gg + ggtitle(names(exprSet)[[i]]) + scale_y_continuous(limits = c(0,150000))
          print(paste(names(exprSet)[[i]] , "in large scale"))
        }else{
          SDplot <-vsn::meanSdPlot(as.matrix(exprSet[[i]]),plot=FALSE ,bins  = 200) 
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

# remove the count matrice to have "same range" data for plot
exprSet_nocts <-exprSet[-1]        
      
pdf("after_transf_violinPlot_4genes.pdf",height = 4, width = 6)   
random_genes <-sample(row.names(as.data.frame(exprSet_nocts[["array"]][1])),size = 4,replace =F)
      
lgExp <-(length(exprSet_nocts))-1
      
      for (j in random_genes){
        gene_df <-do.call(rbind, (lapply(exprSet_nocts[1:lgExp], function(x) x[j,])))
        gene_df<- as.data.frame(t(gene_df))
        gene_df$array <-as.numeric(exprSet_nocts[["array"]][j,])
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
        densarray <- density(exprSet_nocts[["array"]][,j])
        ind_df$array <- densarray$y
        p<-violinPlot_Ylim(ind_df, mainTitle=paste0("pat = ",j),
                           dotsize=0,binwidth=0.2,Xangle=45, ylo = -0.1, yhi = 0.5)      
  print (p)
}
      
      
dev.off()
      


 # the end of RNA transfo - see step 2 for the crossplatform normalization (CPN) step
                                        
                                        
                                        
