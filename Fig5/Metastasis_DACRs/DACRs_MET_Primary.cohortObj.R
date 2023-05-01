library(Signac)
library(Seurat)
library(org.Hs.eg.db)
library(reshape)
library(reshape2)
library(EnsDb.Hsapiens.v86)
library(plyr)
library(dplyr)

library(future)
plan("multiprocess", workers =20)
options(future.globals.maxSize = 50 * 1024^3)

path_to_ATAC_cell_type_annotation=''
path_to_ATAC_catalog=''
path_to_atac_cohort_obj=''

meta=read.table(path_to_ATAC_cell_type_annotation,sep=''),sep='\t',header=T)
rownames(meta)=meta[,1]


for (disease in c('PDAC','CRC','SKCM','UCEC')){

atac=readRDS(path_to_atac_cohort_obj')

meta_s=meta[meta$Cancer==disease,]
orig_1=as.data.frame(atac$Piece_ID)
meta_s=meta_s[rownames(orig_1),]
atac$cell_type.detailed=meta_s$cell_type.harmonized.cancer


samples=read.table(path_to_ATAC_catalog,sep='\t',header=T)
s_sampl=samples[,c('Piece_ID','Disease.Type','Sample.Type')]
s_sampl$Piece_ID_2=paste(s_sampl$Disease.Type,s_sampl$Piece_ID,sep='_')


atac$test='Other'
can=disease
atac$test=ifelse(atac$cell_type.detailed=='Tumor' & atac$Piece_ID %in% s_sampl$Piece_ID_2[s_sampl$Disease.Type==can &
s_sampl$Sample.Type=='Tumor'], paste(can,'Primary_tumor',sep='__'),atac$test)
atac$test=ifelse(atac$cell_type.detailed=='Tumor' & atac$Piece_ID %in% s_sampl$Piece_ID_2[s_sampl$Disease.Type==can &
s_sampl$Sample.Type=='Met'], paste(can,'Met_tumor',sep='__'), atac$test)

ATAC=atac
rm(atac)

###Try with FindMarkers:
###2022-12-28: multiome samples doesn't have 'passed_filters' metadata.

peak.data <- GetAssayData(object = ATAC, assay = 'pancan', slot = "counts")
peak.counts <- colSums(x = peak.data)
ATAC <- AddMetaData(object = ATAC, metadata = peak.counts, col.name = 'peak_RF_pancan_s')
DefaultAssay(ATAC)='pancan'




Idents(ATAC)=ATAC$test
all_da_peaks=NULL
cell_t1=paste(disease,'Met_tumor',sep='__')
cell_t2=paste(disease,'Primary_tumor',sep='__')
print(cell_t1)

da_peaks <- FindMarkers(
  object = ATAC,
  ident.1 = cell_t1,
  ident.2=cell_t2,
  only.pos = FALSE,
  min.pct = 0.01,
  min.diff.pct=0,
  logfc.threshold=0,
  test.use = 'LR',
  latent.vars = 'peak_RF_pancan_s'
)

da_peaks$Cell_t1=cell_t1
da_peaks$Disease=disease
da_peaks$peak=rownames(da_peaks)
all_da_peaks=rbind(all_da_peaks,da_peaks)

write.table(all_da_peaks, paste("out/da_peaks_Met_primary_CohortObj.",can,".tsv",sep=""),
sep="\t",quote=FALSE,row.names=FALSE)
}
