# Alla Karpova create pancan combo object with RNA and ATAC

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =4)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))



######################
### FUNCTIONS #####

normalize_multiome <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj <- SCTransform(obj, 
                       vars.to.regress = c("nCount_RNA","percent.mito"),
                       return.only.var.genes = T, verbose = FALSE, conserve.memory = TRUE) %>% 
    RunPCA(npcs = 50, verbose = FALSE) %>% 
    RunUMAP(dims = 1:50, reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
  
  DefaultAssay(obj) <- "pancan"
  
  obj <- obj %>% 
    RunTFIDF() %>%
    FindTopFeatures( min.cutoff = 20) %>%
    RunSVD(
      reduction.key = 'LSI_',
      reduction.name = 'lsi',
      irlba.work = 400) %>%
    RunUMAP(dims = 2:50,reduction = 'lsi', reduction.name = "atac.umap", reduction.key = "atacUMAP_")
  
  obj <- FindMultiModalNeighbors(obj, 
                                   reduction.list = list("pca", "lsi"), 
                                   dims.list = list(1:50, 2:50))
  obj <- RunUMAP(obj, nn.name = "weighted.nn", 
                   reduction.name = "wnn.umap", 
                   reduction.key = "wnnUMAP_")
  obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3, verbose = T)
  
  return(obj)
}


###options###
######################
option_list = list(
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("--metadata.rna"),
              type="character",
              default='/diskmnt/Projects/snATAC_primary/04_celltyped_rds/cell_type_snRNA_merged/v7.0_data_freeze/All_snRNA_samples_metadata_data_freeze_v7.0.tsv', 
              help="path to rna metadata"),
  make_option(c("--metadata.atac"),
              type="character",
              default='/diskmnt/Projects/snATAC_primary/05_merged_rds/update_cell_types/v7.0_data_freeze_ATAC/All_225_samples_metadata_data_freeze_v7.0.tsv', 
              help="path to atac metadata"),
  make_option(c("--cancer"),
              type="character",
              default=NULL, 
              help="BRCA, CRC, CESC, PDAC, OV, UCEC")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
out_path <- opt$output
meta.rna.path <- opt$metadata.rna
meta.atac.path <- opt$metadata.atac
cancer = opt$cancer

dir.create(out_path, showWarnings = F)
setwd(out_path)

dir.create('out')
dir.create('plots')

rna.path <- paste0('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/snRNA/Merged_objects/', cancer, '/', cancer, '_merge_obj_no_doublets_v7.rds')
atac.path <- paste0('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/snATAC/Merged_objects_PanCan_peaks/RDS.withUMAPs.50PCs.noDoublets/', cancer, '_snATAC_Merged.PancanSet.noDoublets.20230114.rds')

annotation <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')
colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/Colors_panatac_v1.0.rds')

rna.meta <- fread(meta.rna.path) %>% data.frame(row.names = 1)
atac.meta <- fread(meta.atac.path) %>% data.frame(row.names = 1)


rna.meta.renamed <- rna.meta %>% 
  dplyr::select(-orig.ident, -Piece_ID, -Cancer) %>% 
  rename(cell_type.harmonized.cancer.rna = cell_type.harmonized.cancer, 
        seurat_clusters.rna = seurat_clusters) 
print(head(rna.meta.renamed))

atac.meta.renamed <- atac.meta %>% 
  dplyr::select( -UMAP_1, -UMAP_2, -Sample_type, -Chemistry) %>% 
  rename(
         seurat_clusters.atac = seurat_clusters) 
print(head(atac.meta.renamed))

combo.meta <- merge(atac.meta.renamed,rna.meta.renamed, by = 0) %>% 
  data.frame(row.names = 1) %>%
  filter(data.type == '10x_SC_Multi_ATAC_SEQ')

print(head(combo.meta))
  print(cancer)
  
if (file.exists(paste0('out/', cancer, '_all_cells_multiome_obj.20230120.rds'))) {
  r.obj <- readRDS(paste0('out/', cancer, '_all_cells_multiome_obj.20230120.rds'))
} else {
  

  cat('opening objects \n')
  a.obj <- readRDS(atac.path)
  r.obj <- readRDS(rna.path)
  cat('done \n')
  
  Annotation(a.obj) <-annotation 
 
  rna.cells <- rownames(rna.meta %>% filter(data.type.atac == '10x_SC_Multi_ATAC_SEQ' & Cancer == cancer))
  atac.cells <- rownames(atac.meta %>% filter(data.type == '10x_SC_Multi_ATAC_SEQ' & Cancer == cancer))
  combo.cells <- rownames(combo.meta %>% filter(data.type.atac == '10x_SC_Multi_ATAC_SEQ' & Cancer == cancer))
  
  a.obj <- subset(a.obj, cells = combo.cells)
  r.obj <- subset(r.obj, cells = combo.cells)
  
  a.obj <- AddMetaData(a.obj, combo.meta)
  r.obj <- AddMetaData(r.obj, combo.meta)
  
  r.obj[['pancan']] <- a.obj[['pancan']]
  
  cat('normalizing all \n')
  r.obj <- normalize_multiome(r.obj)
  cat('done \n')
  
  saveRDS(r.obj, paste0("out/",cancer,"_all_cells_multiome_obj.", format(Sys.Date(), format="%Y%m%d"),".rds"))
}  
  p1 <- DimPlot(r.obj, reduction = "rna.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
  p2 <- DimPlot(r.obj, reduction = "atac.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
  p3 <- DimPlot(r.obj, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  
  pdf(paste0("plots/",cancer,"_all_cell_types_seurat_clusters.pdf"),height=8,width=24)
  p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  dev.off()
  
  DimPlot(r.obj, reduction = "wnn.umap", group.by = "cell_type.harmonized.cancer.rna", label = TRUE, label.size = 2.5, repel = TRUE, cols = colors$cell_type)
  ggsave(paste0("plots/",cancer,"_all_cell_types_cell_type.harmonized.cancer.rna.pdf"), height=8,width=11)
  
  DimPlot(r.obj, reduction = "wnn.umap", group.by = "cell_type.harmonized.cancer.atac", label = TRUE, label.size = 2.5, repel = TRUE, cols = colors$cell_type)
  ggsave(paste0("plots/",cancer,"_all_cell_types_cell_type.harmonized.cancer.atac.pdf"), height=8,width=11)
  
  DimPlot(r.obj, reduction = "wnn.umap", group.by = "Piece_ID", label = TRUE, label.size = 2.5, repel = TRUE)
  ggsave(paste0("plots/",cancer,"_all_cell_types_Piece_ID.pdf"), height=8,width=11)
  
  #### subset tumor and don't include Normal samples ###########
  cat('subsetting tumor \n')
  r.obj <- subset(r.obj, cell_type.harmonized.cancer.atac == "Tumor" & cell_type.harmonized.cancer.rna == "Tumor" & Sample_type != 'Normal')
  cat('normalizing tumor \n')
  r.obj <- normalize_multiome(r.obj)
  cat('done \n')
  saveRDS(r.obj, paste0("out/",cancer,"_atac.tumor_cells_multiome_obj.", format(Sys.Date(), format="%Y%m%d"),".rds"))
  
  p1 <- DimPlot(r.obj, reduction = "rna.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
  p2 <- DimPlot(r.obj, reduction = "atac.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
  p3 <- DimPlot(r.obj, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  
  pdf(paste0("plots/",cancer,"_atac.tumor_cell_types_seurat_clusters.pdf"),height=8,width=24)
  p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  dev.off()
  
  DimPlot(r.obj, reduction = "wnn.umap", group.by = "cell_type.harmonized.cancer.rna", label = TRUE, label.size = 2.5, repel = TRUE, cols = colors$cell_type)
  ggsave(paste0("plots/",cancer,"_atac.tumor_cell_types_cell_type.harmonized.cancer.rna.pdf"), height=8,width=11)
  
  DimPlot(r.obj, reduction = "wnn.umap", group.by = "cell_type.harmonized.cancer.atac", label = TRUE, label.size = 2.5, repel = TRUE, cols = colors$cell_type)
  ggsave(paste0("plots/",cancer,"_atac.tumor_cell_types_cell_type.harmonized.cancer.atac.pdf"), height=8,width=11)
  
  DimPlot(r.obj, reduction = "wnn.umap", group.by = "Piece_ID", label = TRUE, label.size = 2.5, repel = TRUE)
  ggsave(paste0("plots/",cancer,"_atac.tumor_cell_types_Piece_ID.pdf"), height=8,width=11)
  


















