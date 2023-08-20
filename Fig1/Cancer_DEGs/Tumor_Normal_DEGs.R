library(Signac)
library(Seurat)
library(presto)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(reshape)
library(reshape2)
library(EnsDb.Hsapiens.v86)
library(tidyverse)
RhpcBLASctl::blas_set_num_threads(50)
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 50 * 1024^3)

path_to_merged_rna_obj=''
path_to_RNA_cell_type_annotation=''
path_to_ATAC_catalog=''

rna=readRDS(path_to_merged_rna_obj)

meta=read.table(path_to_RNA_cell_type_annotation,sep='\t',header=T)
rownames(meta)=meta$X
meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)

orig_1=as.data.frame(rna$Piece_ID)
meta=meta[rownames(orig_1),]
rna$cell_type.detailed=meta$cell_type.detailed


samples=read.table(path_to_ATAC_catalog,sep='\t',header=T)
s_sampl=samples[,c('Piece_ID','Disease.Type','Sample.Type')]
s_sampl$Piece_ID_2=paste(s_sampl$Disease.Type,s_sampl$Piece_ID,sep='_')
met_samples=s_sampl[s_sampl$Sample.Typ=='Met',]


rna$Disease=ifelse(rna$Piece_ID %in% c("HT268B1-Th1H3", "HT029B1-S1PC", "HT035B1-S1PA",
 "HT1408-06","HT141B1-S1H1", "HT206B1-S1H4", "HT271B1-S1H3",
"HT378B1-S1H1", "HT378B1-S1H2", "HT384B1-S1H1", "HT517B1-S1H1"), "BRCA_Basal", rna$Disease)



rna$test=ifelse(rna$cell_type.detailed=='Tumor' & !(rna$Piece_ID %in% c(met_samples$Piece_ID, 'CE332E1-N1','CE336E1-S1',
'CE354E1-S1','CE357E1-S1','CE507-C1A2','PM565P1-T1N1')),
paste(rna$cell_type.detailed,rna$Disease,sep='__'),rna$cell_type.detailed)

Idents(rna)=rna$test


all_degs=NULL

cancers=c('ccRCC','PDAC','UCEC','CRC','BRCA','BRCA_Basal','GBM','MM','OV','SKCM','HNSCC','CESC')
for (can in cancers){
    if (can %in% c('ccRCC')){
    origin_c=c('Proximal Tubule')
}else if (can=='BRCA'){
    origin_c=c("Luminal mature")
}else if (can=='BRCA_Basal'){
    origin_c=c('Luminal progenitor')
}else if (can=='GBM'){
    origin_c=c('OPC')
}else if (can=='PDAC'){
    origin_c=c('Ductal-like2')
}else if (can=='UCEC'){
    origin_c=c('Secretory Endometrial epithelial cells')
}else if (can=='CRC'){
    origin_c=c('Distal Stem Cells')
}else if (can=='MM'){
    origin_c=c('B-cells')
}else if (can=='OV'){
    origin_c=c('Secretory Endometrial epithelial cells')
}else if (can=='HNSCC'){
    origin_c=c('Normal squamous cells')
}else if (can=='CESC'){
    origin_c=c('Normal squamous cells')
}else if (can=='SKCM'){
    origin_c=c('Melanocytes')
}

print (can)
    cell_t1=paste('Tumor',can,sep='__')
    for (cell_t2 in origin_c){

print (paste(cell_t1,cell_t2,sep=' '))

degs<- FindMarkers(
  object = rna,
  ident.1 = cell_t1,
  ident.2 = cell_t2,
  only.pos = FALSE,
  min.pct = 0.05,
  min.diff.pct=0,
  assay='SCT',
  logfc.threshold=0,
)

degs$Disease=can
degs$Gene=rownames(degs)
degs$cell_t1=cell_t1
degs$cell_t2=cell_t2
all_degs=rbind(all_degs,degs)

}
}

all_degs=as.data.frame(all_degs)

write.table(all_degs, paste("out/degs_Tumor_Normal.minPct0.05.logFC.0.20230318.tsv",sep=""),
sep="\t",quote=FALSE,row.names=FALSE)

