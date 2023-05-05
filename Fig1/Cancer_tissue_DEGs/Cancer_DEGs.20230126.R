library(Signac)
library(Seurat)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
RhpcBLASctl::blas_set_num_threads(90)
library(future)
plan("multiprocess", workers = 20)
options(future.globals.maxSize = 50 * 1024^3)

path_to_merged_rna_obj=''
path_to_RNA_cell_type_annotation='

rna=readRDS(path_to_merged_rna_obj)

meta=read.table(path_to_RNA_cell_type_annotation,sep='\t',header=T)
rownames(meta)=meta$X
meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)

orig_1=as.data.frame(rna$Piece_ID)
meta=meta[rownames(orig_1),]
rna$cell_type.detailed=meta$cell_type.detailed

rna_all=rna
rna=subset(rna_all, cell_type.detailed=='Tumor')

rna$Piece_ID=paste(rna$Disease,rna$Piece_ID, sep='_')

rna$Disease=ifelse(rna$Piece_ID %in% c("BRCA_HT268B1-Th1H3", "BRCA_HT029B1-S1PC", "BRCA_HT035B1-S1PA",
 "BRCA_HT1408-06","BRCA_HT141B1-S1H1", "BRCA_HT206B1-S1H4", "BRCA_HT271B1-S1H3",
"BRCA_HT378B1-S1H1", "BRCA_HT378B1-S1H2", "BRCA_HT384B1-S1H1", "BRCA_HT517B1-S1H1"), "BRCA_Basal",
rna$Disease)

cancers=unique(rna$Disease)
Idents(rna)=rna$Disease
DefaultAssay(rna)='SCT'
all_degs=NULL

for (disease in cancers){
degs <- FindMarkers(
  object = rna,
  ident.1 = disease,
#  ident.2='',
  only.pos = TRUE,
  min.pct = 0.1,
  min.diff.pct=0,
  assay='SCT',
  logfc.threshold=0
)

degs$Disease=disease
degs$Gene=rownames(degs)
all_degs=rbind(all_degs,degs)
print(disease)
}
write.table(all_degs, paste("out/degs_oneCances_vs_Others.minPct0.1.20230126.tsv",sep=""),
sep="\t",quote=FALSE,row.names=FALSE)