

#!/usr/bin/Rscript

#==========================================================================
#==========================================================================


 
#     STEP 4:  heatmap all results   



#==========================================================================
#==========================================================================

setwd("tcga")

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

goodmatch <- readRDS("goodmatch.rds")
load("corcoef_batch.rda")
load("corcoef_sub.rda")
corPltf3 <- readRDS ("corPltf3.rds")
delta_gPCA <- readRDS (load("delta_gPCA.rds")


row.names(goodmatch)<-gsub("_clusters","",row.names(goodmatch) )
row.names(goodmatch)<-gsub("TDM_RBE","exprSet_TDM_RBE",row.names(goodmatch) )
row.names(goodmatch)<-gsub("^FSQN","exprSet_FSQN",row.names(goodmatch) )
row.names(goodmatch)<-gsub("^RBE","exprSet_RBE",row.names(goodmatch) )
row.names(goodmatch)<-gsub("^TDM","exprSet_TDM",row.names(goodmatch) )


# correlation scaling 

cor_med<-apply(corPltf3,2,FUN=median)
names(cor_med)<-colnames(corPltf3)


# dendrogrram pairing scaling

dend_pair<-goodmatch$perc_goodmach
names(dend_pair)<-rownames(goodmatch)
dend_pair <-dend_pair[names(cor_med)]

# PC axe corrlation scaling

PCAcor_batch<- 100-((corcoef_batch$PC1)*100)
names(PCAcor_batch)<-rownames(corcoef_batch)

# gPCA scaling

PCAdelta_batch<- round(100-((delta_gPCA$delta_batch)*100),digits=2)
names(PCAdelta_batch)<- rownames(delta_gPCA)


#===========================================



# heatmap 


#============================================

#=================
# prepa matrix
#=================

cbind(names(dend_pair),names(cor_med),names(PCAcor_batch), names(PCAdelta_batch))
ht_Mat_TCGA <-data.frame(dend_pair=dend_pair,cor_ind=cor_med*100, PCA_cor=PCAcor_batch,PCA_delta=PCAdelta_batch,
                         row.names = gsub("exprSet_","",names(dend_pair)),check.rows = T)
names(ht_Mat_TCGA) <- paste0("TCGA_", names(ht_Mat_TCGA) )
ht_Mat_TCGA <-as.data.frame(t(ht_Mat_TCGA))

#=================
# prepa sample annot for heatmap
#=================

PCAcor_sub<- ((corcoef_sub$PC1)*100)
names(PCAcor_sub)<-rownames(corcoef_sub)

PCAdelta_sub<- round(((delta_gPCA2$delta_sub)*100),digits=2)
names(PCAdelta_sub)<- rownames(delta_gPCA2)

cbind(names(PCAcor_sub), names(PCAdelta_sub))  #check order

ht_sampleAnnot_TCGA <-data.frame(cancertype_cor=PCAcor_sub,cancertype_delta=PCAdelta_sub,
                                 row.names = gsub("exprSet_","",names(PCAdelta_sub)),check.rows = T)
row.names(ht_sampleAnnot_TCGA)<-gsub("fpkm","xpkm",row.names(ht_sampleAnnot_TCGA))
names(ht_sampleAnnot_TCGA) <-paste0("TCGA_", names(ht_sampleAnnot_TCGA))

ht_sampleAnnot <- ht_sampleAnnot_TCGA

#=========================
#  clustering
#=========================


hclust_sample<-pvclust::pvclust(ht_Mat,nboot=10,method.dist = "euclidean", #use.cor = "complete.obs",
                                method.hclust = "ward.D2",parallel = TRUE)


hclust_gene<-pvclust::pvclust(t(ht_Mat),nboot=10,method.dist = "euclidean",
                              method.hclust = "ward.D2",parallel = TRUE)


#===================
#  up annotation
#====================

names(ht_sampleAnnot)
annot_up<-ht_sampleAnnot


ha_up<-HeatmapAnnotation(df=as.data.frame (annot_up), annotation_height = unit(1, "mm"),show_legend = F,
                         annotation_legend_param =list(ncol = 1),
                         show_annotation_name = TRUE, annotation_name_gp = gpar(fontsize=8),
                         col= list(TCGA_cancertype_cor = circlize::colorRamp2(c(0,50,100), c("#380091","#ff9d00",  "#b8070d")),
                                   TCGA_cancertype_delta = circlize::colorRamp2(c(0,50,100), c("#380091","#ff9d00", "#b8070d")) ) )


#===================
#  annot du bas
#====================

annot_bot<- data.frame(RNAseq_trans=rownames(ht_sampleAnnot),crossPtf_algo=rownames(ht_sampleAnnot),
                       row.names = rownames(ht_sampleAnnot),check.rows = T)
for(i in c("log","vst2", "xpkm","tpm","rle","tmm")){
  annot_bot$RNAseq_trans[annot_bot$RNAseq_trans %in% grep(pattern=i,annot_bot$RNAseq_trans,value=T)]<-i
}
annot_bot$RNAseq_trans[annot_bot$RNAseq_trans %in% grep(pattern="*vst$",annot_bot$RNAseq_trans,value=T)]<-"vst"     #for the mix betzween vst and vst2


for(i in c("log","vst$", "vst2", "xpkm","tpm","rle","tmm")){
  annot_bot$crossPtf_algo <- gsub(i,"",annot_bot$crossPtf_algo,)
  annot_bot$crossPtf_algo <- gsub("_$","",annot_bot$crossPtf_algo,)
}

annot_bot$crossPtf_algo[annot_bot$crossPtf_algo %in% ""]<-"w/o norm"
annot_bot<-droplevels.data.frame(annot_bot)


colTopannot_bot<-vector("list", ncol(annot_bot))
names(colTopannot_bot)<-colnames(annot_bot)
colFun<-c(rainbow,rainbow,heat.colors,rainbow,rainbow,terrain.colors,rainbow)  #mettre autant de couleur que de doncol de ha
i<-1
for(col in colnames(annot_bot)){
  colTopannot_bot[[col]]<-colFun[[i]](nlevels(as.factor(annot_bot[,col])))
  names(colTopannot_bot[[col]])<-levels(as.factor(annot_bot[,col]))
  i<-i+1
}

ha_do<-HeatmapAnnotation(df=as.data.frame(annot_bot),col = colTopannot_bot, na_col = "white",
                         show_annotation_name = T, annotation_name_gp = gpar(fontsize=8),
                         annotation_legend_param =list(ncol = 1),                                                                  
                         annotation_height = unit.c(unit(3, "mm"), unit(3, "mm"), 
                                                    max_text_width(rownames(ht_Mat)) + unit(5, "mm")),                
                         colname = anno_text(names(ht_Mat),gp = autoGparFontSizeMatrix(nrow(ht_Mat)),rot = 90, 
                                             just = "right", location = unit(1, "npc") - unit(1, "mm")))


#===================
#  right annot
#====================


annot_right<- data.frame(dataset=rownames(ht_Mat),perf_metric=rownames(ht_Mat),
                         row.names = rownames(ht_Mat),check.rows = T, stringsAsFactors=F)
for(i in c("CCLE","TCGA")){
  annot_right$dataset[annot_right$dataset %in% grep(pattern=i,annot_right$dataset,value=T)]<-i
}

for(i in c("cor_ind","dend_pair", "PCA_cor", "PCA_delta")){
  annot_right$perf_metric[annot_right$perf_metric %in% grep(pattern=i,annot_right$perf_metric,value=T)]<-i
}

annot_right<-droplevels.data.frame(annot_right)

colannot_right<-vector("list", ncol(annot_right))
names(colannot_right)<-colnames(annot_right)
colFun<-c(topo.colors,topo.colors,heat.colors,heat.colors,rainbow,rainbow,terrain.colors,rainbow)  #mettre autant de couleur que de doncol de ha
i<-1
for(col in colnames(annot_right)){
  colannot_right[[col]]<-colFun[[i]](nlevels(as.factor(annot_right[,col])))
  names(colannot_right[[col]])<-levels(as.factor(annot_right[,col]))
  i<-i+1
}

ha_right<-rowAnnotation(df=as.data.frame(annot_right), na_col = "white",
                        col=colannot_right, 
                        show_legend = T, show_annotation_name = F,annotation_name_gp = gpar(5), 
                        annotation_legend_param =list(ncol = 1)     ,                                                             
                        annotation_width = unit(1, "mm")  )  

#===========================================
#     plot heatmap
#======================================




Ht_pv<-Heatmap(as.matrix(ht_Mat),name="% algo perf", na_col = "white",
               col = circlize::colorRamp2(c(0, 50, 100), c("#4F81BD", "yellow","#C0504D")),       # "#FFEB84"
               bottom_annotation = ha_do,
               top_annotation = ha_up,
               right_annotation = ha_right,
               show_column_names = F,
               heatmap_legend_param = list(ncol=1),
               row_names_gp =  autoGparFontSizeMatrix(nrow(ht_Mat)*10),
               column_names_gp = autoGparFontSizeMatrix(ncol(ht_Mat)*10), 
               cluster_rows = hclust_gene$hclust,
               cluster_columns = hclust_sample$hclust,
               column_title = "")


pdf("Ht_Algo_Perf_both.pdf",width = 15, height = 6)
print(Ht_pv)
dev.off()  

col_fun =  circlize::colorRamp2(c(0,50,100), c("#380091","#ff9d00",  "#b8070d"))
lgd = Legend(col_fun = col_fun, title = "cancerType_metrics", at = c(0,20,40, 60, 80, 100), 
             legend_height = unit(4, "cm"))

draw(lgd)
dev.off()
                         
                         
                         
#   the end of WP1


