



#!/usr/bin/Rscript

#==========================================================================
#==========================================================================



#     STEP 4:  heatmap all results   



#==========================================================================
#==========================================================================


library(ComplexHeatmap)
library(pvclust)

# custom function

autoGparFontSizeMatrix<-function(n){ # automatic caluclation of font size for heatmap ( from https://gitlab.univ-nantes.fr/E114424Z)
  n=min(n,1000)
  return(grid::gpar(fontsize=1/n*600))
}


#=================================================


# load resuslts and metric scaling


#==================================================

setwd("tcga")

goodmatch <- readRDS (file="goodmatch.rds")      # dendro stat
corPltf_all <- readRDS (  file="corPltf_all.rds")   # correlation stat
corcoef_3D<- readRDS ("corcoef_3D.rds")             # PCA correlation stat
delta_gPCA <- readRDS ("delta_gPCA_all.rds")     # gPCA  stat
                       
fileNames <- readRDS("fileNames.rds")                       


# Sample correlation scaling 
#==========================

cor_med <- list()
for (i in fileNames)  {
  
  cor_med[[i]]<-apply(corPltf_all[[i]],2,FUN=median)
  cor_med[[i]] <- 100*cor_med[[i]]
  names(cor_med[[i]])<-colnames(corPltf_all[[i]])  
 
}                     
                      
cor_med <- as.data.frame(do.call(c,cor_med)   )                   
rownames(cor_med) <- gsub(".*\\.", "", rownames(cor_med))  # remove what is before the period  
names(cor_med) <- "sample_corr"


# dendrogrram pairing scaling
#==========================

dend_pair <- list()
for (i in fileNames)  {                       
   dend_pair[[i]]<-goodmatch[[i]]$perc_goodmach
   names(dend_pair[[i]])<-rownames(goodmatch[[i]])                    
}

dend_pair <- as.data.frame(do.call(c,dend_pair)   )                   
rownames(dend_pair) <- gsub(".*\\.", "", rownames(dend_pair))  # remove what is before the period  
names(dend_pair) <- "dend_pairing"

                 
# PCA correlation scaling
#==========================

PCAcorr <- list()
for (i in fileNames)  {   
  PCAcorr[[i]] <- corcoef_3D[[i]]$PC1[grep("batch",rownames(corcoef_3D[[i]]$PC1)), "R2"]
  PCAcorr[[i]]<- 100-((PCAcorr[[i]])*100)
  names(PCAcorr[[i]]) <- rownames(corcoef_3D[[i]]$PC1)[grep("batch",rownames(corcoef_3D[[i]]$PC1))]
  names(PCAcorr[[i]]) <- gsub("dim1_|_batch","",names(PCAcorr[[i]]))
} 
 
PCAcorr <- as.data.frame(do.call(c,PCAcorr)   )                   
rownames(PCAcorr) <- gsub(".*\\.", "", rownames(PCAcorr))  # remove what is before the period  
names(PCAcorr) <- "PCAcorr"

                       
# gPCA scaling
#===============
                       
PCAdelta_batch <- data.frame(gPCA = round(100-((delta_gPCA$delta_batch)*100),digits=2),
                           row.names=rownames(delta_gPCA) )

                       
 #===========================================
 
 
 # heatmap 
 
 
 #============================================
 
 #=================
 # prepa matrix
 #=================

ht_Mat <-Reduce(merge, lapply(list(cor_med,dend_pair,PCAcorr,PCAdelta_batch), function(x) data.frame(x, rn = row.names(x))))
rownames(ht_Mat) <-gsub("exprSet_","",ht_Mat$rn) ; ht_Mat$rn <- NULL
ht_Mat<-as.data.frame(t(ht_Mat))

                                         
 #=========================
 #  clustering
 #=========================
 
 
 hclust_row<-pvclust::pvclust(ht_Mat,nboot=10,method.dist = "euclidean", #use.cor = "complete.obs",
                                 method.hclust = "ward.D2",parallel = TRUE)
 
 
 hclust_col<-pvclust::pvclust(t(ht_Mat),nboot=10,method.dist = "euclidean",
                               method.hclust = "ward.D2",parallel = TRUE)
 
 
 #===================
 #  down_annot
 #====================
 
 annot_bot<- data.frame(RNAseq_trans= gsub(".*_", "",rownames(ht_Mat)) ,
                        crossPtf_algo=gsub("_.*", "",rownames(ht_Mat)) ,
                        row.names = rownames(ht_Mat),check.rows = T)

 annot_bot<-droplevels.data.frame(annot_bot)
 
 
 colTopannot_bot<-vector("list", ncol(annot_bot))
 names(colTopannot_bot)<-colnames(annot_bot)
 for(col in colnames(annot_bot)){
   colTopannot_bot[[col]]<-rainbow(nlevels(as.factor(annot_bot[,col])))
   names(colTopannot_bot[[col]])<-levels(as.factor(annot_bot[,col]))

 }
 
 ha_do<-HeatmapAnnotation(df=as.data.frame(annot_bot),col = colTopannot_bot, na_col = "white",
                          show_annotation_name = T, annotation_name_gp = gpar(fontsize=8),
                          annotation_legend_param =list(ncol = 1),                                                                  
                          annotation_height = unit.c(unit(3, "mm"), unit(3, "mm"), 
                                                     max_text_width(rownames(ht_Mat)) + unit(5, "mm")),                
                          colname = anno_text(rownames(ht_Mat),gp = autoGparFontSizeMatrix(nrow(ht_Mat)),rot = 90, 
                                              just = "right", location = unit(1, "npc") - unit(1, "mm")))
 
 
 #===================
 #  right annot
 #====================
 
 
 annot_right<- data.frame(dataset=rep("TCGA",ncol(ht_Mat)),
                          perf_metric=colnames(ht_Mat),
                          row.names = colnames(ht_Mat),
                          check.rows = T, stringsAsFactors=F)

 annot_right<-droplevels.data.frame(annot_right)
 
 colannot_right<-vector("list", ncol(annot_right))
 names(colannot_right)<-colnames(annot_right)
 for(col in colnames(annot_right)){
   colannot_right[[col]]<-heat.colors(nlevels(as.factor(annot_right[,col])))
   names(colannot_right[[col]])<-levels(as.factor(annot_right[,col]))
 }
 
 ha_right<-rowAnnotation(df=as.data.frame(annot_right), na_col = "white",
                         col=colannot_right, 
                         show_legend = T, annotation_name_side = "top",
                         show_annotation_name = T,annotation_name_gp = gpar(fontsize=8), 
                         annotation_legend_param =list(ncol = 1)     ,                                                             
                         annotation_width = unit(1, "mm")  )  
 
 #===========================================
 #     plot performance heatmap
 #======================================
 
 
 Ht_pv<-Heatmap(t(ht_Mat),name="% algo perf", na_col = "white",
                col = circlize::colorRamp2(c(0, 50, 100), c("#4F81BD", "yellow","#C0504D")),       # "#FFEB84"
                bottom_annotation = ha_do,
                right_annotation = ha_right,
                show_column_names = F,
                heatmap_legend_param = list(ncol=1),
                row_names_gp =  autoGparFontSizeMatrix(nrow((ht_Mat))),
                column_names_gp = autoGparFontSizeMatrix(nrow((ht_Mat))), 
                cluster_rows = hclust_row$hclust,
                cluster_columns = hclust_col$hclust,
                column_title = "")
 
 
 pdf("Ht_Algo_Perf_both.pdf",width = 15, height = 6)
 print(Ht_pv)
 dev.off()  
 

 #   the end of WP1
                       
                       
