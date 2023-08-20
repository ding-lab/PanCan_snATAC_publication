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

out.obj <- str_replace(input_path, pattern = 'rds', replacement = 'diffuse.links.only.rds')
out.obj <- str_split(out.obj, '[/]')[[1]][str_count(out.obj, '/')+1]

cat('opening object \n')
obj <- readRDS(input_path)
DefaultAssay(obj) <- "pancan"
cat('done \n')


annot <- readRDS('/diskmnt/Projects/snATAC_primary/PanCan_ATAC_data_freeze/v2.0/snATAC/merged_no_recalling_upd/Annotations.EnsDb.Hsapiens.v86.rds')
Annotation(obj) <- annot


chr.size=fread('/diskmnt/Projects/snATAC_primary/02_atac_rds_signca/v3.0/hg38.chrom.sizes.txt', header = FALSE,col.names = c('seqnames','chr_length') ) %>% 
  data.frame()
chr.size <- chr.size %>% filter(seqnames %in% paste0('chr', 1:22) | seqnames == 'chrX')
chr.size.vector <- as.numeric(chr.size$chr_length)
names(chr.size.vector) <- chr.size$seqnames
bins <- GenomicRanges::tileGenome(seqlengths = chr.size.vector, tilewidth = 100000, cut.last.tile.in.chrom = TRUE)

frag <- Fragments(obj@assays$pancan)
#this will remove fragment objects with just 1 or 0 cells because they fail FeatureMatrix function
frag.filtered <- frag[do.call( 'c', lapply(frag, function(x) length(x@cells) > 1))]

if (opt$doremoval) {
  #this part finds and removes shitty clusters
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
min.cells.num <- 0.05*ncol(obj)
print(min.cells.num)


# remove pancan assay
DefaultAssay(obj)<-'SCT'
obj[['pancan']] <- NULL

plan("multicore", workers = 30)
options(future.globals.maxSize = 300 * 1024^3) # for 500 Gb RAM

print(length(bins))
pro_n <-  round(length(bins)/30)

exlude.samples <- obj@meta.data %>% group_by(Sample) %>% tally() %>% filter(n<=1) %>% pull(Sample)
#now exclude these samples from the object
obj <- subset(obj, subset = Sample %in% exlude.samples, invert=TRUE)

matrix.counts <- FeatureMatrix(
  fragments = frag.filtered,
  features = bins,
  process_n = pro_n,
  sep = c("-","-"),
  cells = colnames(obj)
)

obj[['X100kb']] <- CreateChromatinAssay(counts = matrix.counts,
                                                 annotation = annot,
                                                 #genome = 'hg38',
                                                 fragments = frag.filtered, 
                                                  min.features = -1)

DefaultAssay(obj)<-'X100kb'


# first compute the GC content for each peak
obj <- RegionStats(obj, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
obj <- LinkPeaks(
  object = obj,
  peak.assay = "X100kb",
  distance = 6e+05,
  score_cutoff = 0,
  n_sample = 1000,
  expression.assay = "SCT"
)

toreturn <- Links(obj)

print(out.obj)
saveRDS(toreturn, paste0('out/', out.obj))
