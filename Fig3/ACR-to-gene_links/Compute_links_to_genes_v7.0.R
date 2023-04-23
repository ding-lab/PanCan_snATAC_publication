# Alla Karpova create pancan combo object with RNA and ATAC

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(future))
plan("multicore", workers =10)
options(future.globals.maxSize = 50 * 1024^3)
suppressMessages(library(optparse))
 
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))


###options###
######################
option_list = list(
  make_option(c("-i", "--input"),
              type="character",
              default=NULL, 
              help="input object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-r", "--doremoval"),
              type="logical",
              default=TRUE, 
              help="Remove multi case clusters? TRUE if you use tumor cells, FALSE if you use normal cells",
              metavar="logical")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input_path <- opt$input
out_path <- opt$output

dir.create(out_path, showWarnings = F)
setwd(out_path)

dir.create('out')
dir.create('plots')

out.obj <- str_replace(input_path, pattern = 'rds', replacement = 'linked.rds')
out.obj <- str_split(out.obj, '[/]')[[1]][str_count(out.obj, '/')+1]

cat('opening object \n')
obj <- readRDS(input_path)
DefaultAssay(obj) <- "pancan"
cat('done \n')

if (opt$doremoval) {
  obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 1, verbose = T, resolution=2)
  
  DimPlot(obj, group.by = 'seurat_clusters', label=TRUE,reduction = "wnn.umap")
  ggsave(paste0('plots/Dimplot_seurat_clusters_res2_', out.obj, '.pdf'), width = 10, height = 10)
  
  DimPlot(obj, group.by = 'Case_ID', label=TRUE,reduction = "wnn.umap")
  ggsave(paste0('plots/Dimplot_Case_ID_', out.obj, '.pdf'), width = 10, height = 10)
  
  
  tb <- table(obj$seurat_clusters, obj$Case_ID) %>% 
    as.data.frame.matrix() 
  tb <- tb/rowSums(tb)
  print(tb)
  
  shared.clusters <- rownames(tb)[rowSums(tb < 0.7) > (ncol(tb)-1)]
  print(shared.clusters)
  
  tohighlight <- rownames(filter(obj@meta.data, seurat_clusters %in% shared.clusters))
  DimPlot(obj, cells.highlight = tohighlight, reduction = "wnn.umap")
  ggsave(paste0('plots/Dimplot_doublet_clusters_', out.obj, '.pdf'), width = 10, height = 10)
  
  print(dim(obj))
  message('subsetting bad clusters out')
  obj <- subset(obj, seurat_clusters %in% shared.clusters, invert=TRUE)
  print(dim(obj))
}

min.cells.num <- 0.00005*ncol(obj)
print(min.cells.num)

# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
obj <- LinkPeaks(
  object = obj,
  peak.assay = "pancan",
  distance = 5e+05,
  n_sample = 1000,
  #min.cells = min.cells.num,
  expression.assay = "SCT"
)


print(out.obj)
saveRDS(obj, paste0('out/', out.obj))
