
#==============================================================================================================


#                WP2 - STEP ONE  /// download, normalize and merge Affymetrix raw data  (exprSet preparation)


#=====================================================================================================================

library(R.utils)  # for gunzip
library(simpleaffy)  # for read and noramlization affy

library(biomaRt)  # gene equivalence
library(dplyr)

#========================================================

#            download and normalize data

#=======================================================

# download
#==========

dir.create("rawCEL_andNorm")
setwd("rawCEL_andNorm")

# download .CEL data from GEO in named folder:
#     Clarke_GSE42568
#     Desmedt07_GSE7390
#     Filipits_GSE26971
#     Hatzis_GSE25055
#     Kao_GSE20685
#     Li_GSE19615
#     Loi07_GSE6532
#     loi08_GSE9195
#     Minn05_GSE2603
#     Minn07_GSE5327
#     Nagalla_GSE45255
#     Pawitan_GSE1456
#     Rody_GSE31519
#     Schmidt_GSE11121
#     Sircoulomb_GSE17907
#     Symmans_GSE17705
#     Wang05_GSE2034    
#     Zhang_GSE12093
#     Zhou_GSE7378


# download .CEL data from ENA
#     Guedj_E_METAB_365

# once done, load and normalize all affy
#=========================================

celFolders <- list.dirs(recursive=F,full.names = F)

for ( i in celFolders){

setwd(i)
  
 # unzip files
files <- list.files(full.names = T,pattern = ".gz",ignore.case=TRUE)
for (j in 1:length(files)){
  gunzip(files[j], remove=FALSE)
}
files <- list.files(full.names = T,pattern = ".cel$",ignore.case=TRUE)

# read affy
celfiles <- ReadAffy(filenames = files)
hist(celfiles)
   
# normalization
celfiles.mas5 = mas5(celfiles) 
hist(exprs(celfiles.mas5))
exprSet.nologs = affy::exprs(celfiles.mas5)
  
# log transfform the data
exprSet.logs <- log(exprSet.nologs, 2)  
hist(exprSet.logs)
rm(files,celfiles,celfiles.mas5,exprSet.nologs) # remove file to be sure there will be no confusion in the next loop in case of issue       
setwd("..")
saveRDS(exprSet.logs, file=paste0(i,".rds")
      
}
        
        
#========================================================

# transform Affy probes in HUGO names and merge datasets

#=======================================================
        
# Download gene information for homo sapiens
ensembl2 <- biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")        
biomaRt::listAttributes(mart=ensembl2)  
        
# Select the column wanted
gene_equi<-biomaRt::getBM(mart=ensembl2,attributes=c("hgnc_symbol","affy_hg_u133_plus_2" ))          
head(gene_equi)        
        
for ( i  in setdiff( celFolders, c("Pawitan", "Guedj", "Loi")) ) {            # setdiff to do particular case after

NormAffy<-readRDS(file= paste0(i,".rds")
NormAffy <-as.data.frame(NormAffy)
NormAffy$ID2  <- gene_equi$hgnc_symbol[match(x=rownames(NormAffy), table=gene_equi$affy_hg_u133_plus_2)]
length(unique(NormAffy$ID2) )                                                        
NormAffy <- NormAffy[,grep(pattern = ".gz",invert = T,names(NormAffy)   )]
names(NormAffy) <-gsub(".CEL","",names(NormAffy))

NormAffy <- aggregate(NormAffy,by=list( NormAffy$ID2 ),FUN= median)
rownames(NormAffy)<-NormAffy$Group.1 ; NormAffy$Group.1 <- NULL  ;NormAffy$ID2 <-NULL
# NormAffy <- NormAffy[-1,]   #remove first line wiht aggregate of unamed genes

saveRDS(NormAffy,file= paste0(i,"_Process.rds")   
 rm( NormAffy)      
        
  }      
        
        
        
        
        
        
