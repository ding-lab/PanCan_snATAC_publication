library(Signac)
library(Seurat)
library(reshape)
library(reshape2)
library(EnsDb.Hsapiens.v86)
library(future)
plan("multiprocess", workers =4)
options(future.globals.maxSize = 50 * 1024^3)


path_to_atac_obj=''
path_to_ATAC_cell_type_annotation=''

atac=readRDS(path_to_atac_obj)

meta=read.table(path_to_ATAC_cell_type_annotation),sep='\t',header=T)

meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)

orig_1=as.data.frame(atac$Piece_ID)
rownames(meta)=meta$X
meta=meta[rownames(orig_1),]
atac$cell_type.detailed=meta$cell_type.detailed



ATAC=subset(atac, cell_type.detailed=='Tumor')
rm(atac)

ATAC$Disease=ifelse(ATAC$Piece_ID %in% c("BRCA_HT268B1-Th1H3", "BRCA_HT029B1-S1PC", "BRCA_HT035B1-S1PA",
 "BRCA_HT1408-06","BRCA_HT141B1-S1H1", "BRCA_HT206B1-S1H4", "BRCA_HT271B1-S1H3",
"BRCA_HT378B1-S1H1", "BRCA_HT378B1-S1H2", "BRCA_HT384B1-S1H1", "BRCA_HT517B1-S1H1"), "BRCA_Basal",
ATAC$Disease)

###Try with FindMarkers:
peak.data <- GetAssayData(object = ATAC, assay = 'pancan', slot = "counts")
peak.counts <- colSums(x = peak.data)
ATAC <- AddMetaData(object = ATAC, metadata = peak.counts, col.name = 'peak_RF_pancan')
DefaultAssay(ATAC)='pancan'


###calculate in parallel:
Idents(ATAC)=ATAC$Disease

cell_groups=c("HNSCC","CESC","UCEC","ccRCC","OV","SKCM","CRC","PDAC","MM","GBM","BRCA_Basal","BRCA")


all_da_peaks=NULL
for (cell_type in cell_groups){
da_peaks <- FindMarkers(
  object = ATAC,
  ident.1 = cell_type,
#  ident.2='', 
  only.pos = TRUE,
  min.pct = 0.1,
  min.diff.pct=0,
  logfc.threshold=0,
  test.use = 'LR',
  latent.vars = 'peak_RF_pancan'
)

da_peaks$Disease=cell_type
da_peaks$peak=rownames(da_peaks)
all_da_peaks=rbind(all_da_peaks,da_peaks)
print(cell_type)
}


write.table(all_da_peaks, 'out/da_peaks_oneCances_vs_Others.minPct0.1.20230126.tsv',sep='\t',quote=F,row.names=F)
