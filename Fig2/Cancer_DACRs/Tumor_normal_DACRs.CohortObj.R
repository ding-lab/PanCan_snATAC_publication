library(Signac)
library(Seurat)
library(reshape)
library(reshape2)
library(EnsDb.Hsapiens.v86)
library(plyr)
library(dplyr)
library(future)
plan("multiprocess", workers =10)
options(future.globals.maxSize = 50 * 1024^3)

all_da_peaks=NULL

path_to_ATAC_cell_type_annotation=''
path_to_ATAC_catalog=''
path_to_atac_cohort_obj=''

meta=read.table(path_to_ATAC_cell_type_annotation, sep='\t',header=T)

meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)
rownames(meta)=meta[,1]

samples=read.table(path_to_ATAC_catalog,sep='\t',header=T)
s_sampl=samples[,c('Piece_ID','Disease.Type','Sample.Type')]
s_sampl$Piece_ID_2=paste(s_sampl$Disease.Type,s_sampl$Piece_ID,sep='_')
met_samples=s_sampl[s_sampl$Sample.Typ=='Met',]


cancers_s=c('ccRCC','PDAC','UCEC','CRC','BRCA','BRCA_Basal','GBM','MM','SKCM','CESC')


cancers_s=c('CESC','PDAC')
for (can in cancers_s){

if(can %in% c('BRCA','BRCA_Basal')){
disease='BRCA'
}else{
disease=can
}

atac=readRDS(path_to_atac_cohort_obj,sep=''))

meta_s=meta[meta$Cancer==disease,]
orig_1=as.data.frame(atac$Piece_ID)
meta_s=meta_s[rownames(orig_1),]
atac$cell_type.detailed=meta_s$cell_type.detailed

atac$Disease=disease
atac$Disease=ifelse(atac$Piece_ID %in% c("BRCA_HT268B1-Th1H3", "BRCA_HT029B1-S1PC", "BRCA_HT035B1-S1PA",
 "BRCA_HT1408-06","BRCA_HT141B1-S1H1", "BRCA_HT206B1-S1H4", "BRCA_HT271B1-S1H3",
"BRCA_HT378B1-S1H1", "BRCA_HT378B1-S1H2", "BRCA_HT384B1-S1H1", "BRCA_HT517B1-S1H1"), "BRCA_Basal",
atac$Disease)


atac$test=ifelse(atac$cell_type.detailed=='Tumor' & !(atac$Piece_ID %in% c(met_samples$Piece_ID_2,'CESC_CE332E1-N1','CESC_CE336E1-S1',
'CESC_CE354E1-S1','CESC_CE357E1-S1','CESC_CE507-C1A2','PDAC_PM565P1-T1N1')),
paste(atac$cell_type.detailed,atac$Disease,sep='__'),atac$cell_type.detailed)

Idents(atac)=atac$test


ATAC=atac

DefaultAssay(ATAC) <- 'pancan'
rm(atac)

peak.data <- GetAssayData(object = ATAC, assay = 'pancan', slot = "counts")
peak.counts <- colSums(x = peak.data)
ATAC <- AddMetaData(object = ATAC, metadata = peak.counts, col.name = 'peak_RF_pancan')
DefaultAssay(ATAC)='pancan'


Idents(ATAC)=ATAC$test

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
    cell_t2=origin_c

print (paste(cell_t1,cell_t2,sep=' '))

da_peaks <- FindMarkers(
  object = ATAC,
  ident.1 = cell_t1,
  ident.2 = cell_t2,
  only.pos = FALSE,
  min.pct = 0.05,
  min.diff.pct=0,
  logfc.threshold=0,
  test.use = 'LR',
  latent.vars = 'peak_RF_pancan'
)

da_peaks$Disease=can
da_peaks$cell_t1=cell_t1
da_peaks$cell_t2=cell_t2
da_peaks$peak=rownames(da_peaks)

}

all_da_peaks=as.data.frame(all_da_peaks)

write.table(all_da_peaks,paste("out/TumorNormal.DACRs.pct.0.01.logFC.tsv",sep=""),
quote=FALSE,sep="\t",row.names=FALSE)