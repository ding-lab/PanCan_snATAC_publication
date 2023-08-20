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
  make_option(c("-i", "--input"),
              type="character",
              default=NULL, 
              help="Input multiome object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("--metadata"),
              type="character",
              default='/diskmnt/Projects/snATAC_primary/04_celltyped_rds/cell_type_snRNA_merged/v7.0_data_freeze/All_snRNA_samples_metadata_data_freeze_v7.0.tsv', 
              help="path to metadata"),
  make_option(c("--cancer"),
              type="character",
              default=NULL, 
              help="cancer type")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input_path <- opt$input
out_path <- opt$output
meta.rna.path <- opt$metadata
cancer <- opt$cancer

dir.create(out_path, showWarnings = F)
setwd(out_path)

dir.create('out')
dir.create('plots')

colors <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v7.0/Colors_panatac_v1.0.rds')

rna.meta <- fread(meta.rna.path) %>% data.frame(row.names = 1)
print(cancer)

cat('opening objects \n')
r.obj <- readRDS(input_path)
cat('done \n')

normal.ct <- c('Normal epithelial cells', 'Ciliated Endometrial epithelial cells' ,'Secretory Endometrial epithelial cells',
               'Acinar', 'Epithelial','Goblet',
               'Luminal mature', 'Luminal progenitor', 'Basal progenitor')
#### subset tumor and don't include Normal samples ###########
cat('subsetting normal \n')
r.obj <- subset(r.obj, cell_type.harmonized.cancer.atac %in% normal.ct & cell_type.harmonized.cancer.rna %in% normal.ct)
cat('normalizing normal \n')
r.obj <- normalize_multiome(r.obj)
cat('done \n')
saveRDS(r.obj, paste0("out/",cancer,"_atac.normal_cells_multiome_obj.", format(Sys.Date(), format="%Y%m%d"),".rds"))

p1 <- DimPlot(r.obj, reduction = "rna.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(r.obj, reduction = "atac.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(r.obj, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

pdf(paste0("plots/",cancer,"_atac.tumor_cell_types_seurat_clusters.pdf"),height=8,width=24)
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

DimPlot(r.obj, reduction = "wnn.umap", group.by = "cell_type.harmonized.cancer.rna", label = TRUE, label.size = 2.5, repel = TRUE, cols = colors$cell_type)
ggsave(paste0("plots/",cancer,"_atac.normal_cell_types_cell_type.harmonized.cancer.rna.pdf"), height=4,width=9)

DimPlot(r.obj, reduction = "wnn.umap", group.by = "cell_type.harmonized.cancer.atac", label = TRUE, label.size = 2.5, repel = TRUE, cols = colors$cell_type)
ggsave(paste0("plots/",cancer,"_atac.normal_cell_types_cell_type.harmonized.cancer.atac.pdf"), height=4,width=9)

DimPlot(r.obj, reduction = "wnn.umap", group.by = "Piece_ID", label = TRUE, label.size = 2.5, repel = TRUE)
ggsave(paste0("plots/",cancer,"_atac.normal_cell_types_Piece_ID.pdf"), height=4,width=9)



















