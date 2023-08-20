library(Signac)
library(Seurat)
library(Matrix)
library(SCENIC)
library(SCopeLoomR)


path_to_merged_rna_obj=''
path_to_RNA_cell_type_annotation=''
wdir=''

data_dir='/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0'

rna=readRDS(path_to_merged_rna_obj)

meta=read.table(path_to_RNA_cell_type_annotation,sep='\t',header=T)


rownames(meta)=meta$X
meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)

orig_1=as.data.frame(rna$Piece_ID)
meta=meta[rownames(orig_1),]
rna$cell_type.harmonized.cancer=meta$cell_type.harmonized.cancer
rna$cell_type.detailed=meta$cell_type.detailed


cell_types=unique(rna$cell_type.detailed)
sel_cell_types=cell_types[!(cell_types %in% c('','Fibroblasts','Macrophages','T-cells','Endothelial',
'DC','Plasma','Unknown','Pericytes','Microglia','NK','Mast','Tregs','Alveolar',
'Hepatocytes','Cholangiocytes','Neurons','Skeletal Muscle','Erythrocytes','Neurons','Skeletal Muscle','Low quality','vSMCs',
'Adipocytes','Immune','Neutrophils','Monocytes','Pre-B-cells','Other_doublets','Smooth muscle'))]

rna_all=rna
rna=subset(rna_all, cell_type.detailed %in% sel_cell_types)


###Downsample cells:

Idents(rna)=paste(rna$Piece_ID,rna$cell_type.detailed,sep='__')
rna_downs=subset(x = rna, downsample = 200)

assay=GetAssay(rna_downs,assay='SCT')
data=assay@data
genes_list=read.csv(hgnc-gene_list.txt',sep='\t')
genes_prot_c=genes_list[genes_list$Locus.type=='Gene with protein product',]

data=data[rownames(data) %in% genes_prot_c$Symbol,]


genes=rownames(data)
genes_r=genes[grepl('^MT-|^RPL|^RPS', genes)]
data=data[!(rownames(data) %in% genes_r),]

cellInfo_1 <- data.frame(Piece_ID=rna_downs$Piece_ID, Cell_type=rna_downs$cell_type.detailed)
exprMatrix=data

dir.create(paste(wdir,"/data",sep=''))
loom <- build_loom(paste(wdir,"/data/Pancan_obj.loom",sep=''), dgem=exprMatrix)
loom <- add_cell_annotation(loom, cellInfo_1)
close_loom(loom)

cellInfo_1$Barcode=row.names(cellInfo_1)
write.table(cellInfo_1,paste(wdir,'/cellInfo_1_PancanObj_200cellsPerPieceID_SelCellTypes.tsv',sep=''),
sep='\t',quote=F,row.names=F)

