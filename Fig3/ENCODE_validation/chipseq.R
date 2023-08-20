library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(GenomicRanges)
library(ggplot2)
library(biomaRt)
library(gridExtra)

# Find the corresponding chip-seq file
metadata <- read.csv("tableS5b.txt", sep="\t", skip=1)
peak_dir <- "./files/"

peak_files <- list.files(peak_dir, pattern="\\.bed\\.gz$", full.names=TRUE) # list all the files in the dir
peak_ids <- gsub("^.*/(\\w+)\\.bed\\.gz$", "\\1", peak_files) # get ids

metadata$Files_found <- sapply(strsplit(metadata$Files, ","), function(x) {
  ids <- gsub("/files/(\\w+)/.*", "\\1", x)
  check <- intersect(ids, peak_ids)
  y <- paste0(peak_dir, check, ".bed.gz")                                                                                                   })

# Create a list of peak data frames for each TF
# Initialize a list to store the peaks for each TF
peaks_by_tf <- list()

# Loop over the files and extract the peaks for each TF
for (i in 1:length(metadata$Files_found)) {
  filename <- metadata$Files_found[i]
  peaks <- readPeakFile(filename)
  tf <- metadata$Target.gene.symbol[i]
  print(tf)
  if (is.null(peaks_by_tf[[tf]])) {
    peaks_by_tf[[tf]] <- peaks
  } else {
    peaks_by_tf[[tf]] <- c(peaks_by_tf[[tf]], peaks)
  }
}

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annotated_peaks_by_tf <- list()
for (tf in names(peaks_by_tf)) {
  gr <- makeGRangesFromDataFrame(peaks_by_tf[[tf]])
  annotated_peaks_by_tf[[tf]] <- gr
}

tf_metadata2 <- aggregate(cbind(Accession, Biosample.term.name, Files_found) ~ Target.gene.symbol, data=metadata, FUN=function(x) paste(unique(x), collapse=", "))

# Read in list of target genes
scenic <- read.table("tableS4a.txt", sep="\t", header=T)
scenic$TF <- gsub("\\(\\+\\)", "", scenic$TF)
tf_metadata2$Gene <- lapply(tf_metadata2$Target.gene.symbol, function(x)
                                scenic[(scenic$TF==x),]$Gene_target)

# Combine with cancer_type info
tfs <- read.csv("tableS4f.txt", sep="\t", header=T)
tfs$ATAC_TF <- gsub("\\(\\+\\)", "", tfs$Regulon)
tf_metadata <- merge(tf_metadata2, tfs, by.x = "Target.gene.symbol", by.y="ATAC_TF", all.x=TRUE)

atac_peaks_by_ct <- list()
# Loop over the files and extract the peaks for each TF
ct <- c('ccRCC', 'BRCA', 'PDAC', 'OV', 'CRC', 'CESC', 'GBM', 'MM','UCEC', 'HNSCC', 'SKCM')
for (i in ct) {
  filename <- paste0("atac.",i, ".bed.gz" )
  peaks <- readPeakFile(filename)
  print(i)
  atac_peaks_by_ct[[i]] <- peaks
}
annotated_atac_peaks_by_ct <- list()
for (ct in names(atac_peaks_by_ct)) {
  gr <- makeGRangesFromDataFrame(atac_peaks_by_ct[[ct]])
  annotated_atac_peaks_by_ct[[ct]] <- gr
}

tf_metadata$number_of_Target_Genes_chip <- "NA"
tf_metadata$number_of_Target_Genes_atac <- "NA"
tf_metadata$number_of_Target_Genes_window_name <- "NA"

# Convert HGNC symbol to Entrez IDs.
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
transcripts <- transcriptsBy(txdb, "gene")
attributes2 = c("entrezgene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "hgnc_symbol")
all_genes <- names(transcripts)
all_genes_coord <- getBM(filters="entrezgene_id", attributes=attributes2, values=all_genes, mart=ensembl)
all_genes_on_chr_coord <- all_genes_coord[all_genes_coord$chromosome_name %in% as.character(c(1:22, "X", "Y")),]
all_genes_on_chr <- as.character(all_genes_on_chr_coord$entrezgene_id)
# fIterate over the TFs and plot the average ChIP-seq signal around their target genes
for (tf in names(annotated_peaks_by_tf)) {
  ct <- tf_metadata$Cancer[tf_metadata$Target.gene.symbol==tf]
  print(paste0("Making plot for ", tf, " in ", ct))
  # Read in the peak files for the current TF
  peaks <- annotated_peaks_by_tf[[tf]]
  atac_peaks <- annotated_atac_peaks_by_ct[[ct]]
  # Obtain the genomic coordinates for the TSS of each target gene
  target_transcripts_coord <- all_genes_on_chr_coord[all_genes_on_chr_coord$hgnc_symbol %in% tf_metadata$Gene[tf_metadata$Target.gene.symbol==tf][[1]], ]
  if (nrow(target_transcripts_coord) >10) {
    target_genes <- GRanges(seqnames = target_transcripts_coord$chromosome_name,
                            ranges = IRanges(start = target_transcripts_coord$start_position, end = target_transcripts_coord$end_position),
                            strand = target_transcripts_coord$strand,
                            names = target_transcripts_coord$entrezgene_id)
    seqlevelsStyle(peaks) <- "Ensembl"
    seqlevelsStyle(atac_peaks) <- "Ensembl"
    seqlevelsStyle(target_genes) <- "Ensembl"
    target_region <- makeBioRegionFromGranges(target_genes, type="start_site", by="gene", upstream=5000, downstream=5000)
    seqlevelsStyle(target_region) <- "Ensembl"
    tagMatrix <- getTagMatrix(peaks, windows=target_region)
    atac_tagMatrix <- getTagMatrix(atac_peaks, windows=target_region)
    if (is.null(ncol(tagMatrix))) {
      print(paste0(tf, " don't have good chip-seq signals from ChIP-seq data"))
      next
    } else if (nrow(tagMatrix) == 0) {
      print(paste0(tf, " has an empty tag matrix"))
      next
    } else if (tf=="NRF1") {
      cut_peak <- readPeakFile("cut.NRF1_U251.bed.gz")
      seqlevelsStyle(cut_peak) <- "Ensembl"
      cut_tagMatrix <- getTagMatrix(cut_peak, windows=target_region)
      tagMatrixList <- list(tagMatrix, atac_tagMatrix, cut_tagMatrix)
      names(tagMatrixList) <- c("ChIP-seq (ENCODE)", paste0("snATAC-seq (", ct, ")"), "Cut & Run (U251)")
      pdf(paste0(tf, "_5000bp_avg_line_plot.pdf"), height=3.5, width=4.5)
      avgplot <- plotAvgProf(tagMatrixList, xlim=c(-5000, 5000),xlab=paste0(tf, " Regulon"), ylab = "Read Count Frequency", facet="row")
      grid.arrange(avgplot,  ncol=1)
      dev.off()
      tf_metadata$number_of_Target_Genes_chip[tf_metadata$Target.gene.symbol==tf] <- nrow(tagMatrix)
      tf_metadata$number_of_Target_Genes_atac[tf_metadata$Target.gene.symbol==tf] <- nrow(atac_tagMatrix)
      tf_metadata$number_of_Target_Genes_window_name[tf_metadata$Target.gene.symbol==tf] <- paste0(target_transcripts_coord$hgnc_symbol, collapse=",") 
    } else {
      tagMatrixList <- list(tagMatrix, atac_tagMatrix)
      names(tagMatrixList) <- c("ChIP-seq (ENCODE)", paste0("snATAC-seq (", ct, ")"))
      pdf(paste0(tf, "_5000bp_avg_line_plot.pdf"), height=3.5, width=4.5)
      avgplot <- plotAvgProf(tagMatrixList, xlim=c(-5000, 5000),xlab=paste0(tf, " Regulon"), ylab = "Read Count Frequency", facet="row")
      grid.arrange(avgplot,  ncol=1)
      dev.off()
      tf_metadata$number_of_Target_Genes_chip[tf_metadata$Target.gene.symbol==tf] <- nrow(tagMatrix)
      tf_metadata$number_of_Target_Genes_atac[tf_metadata$Target.gene.symbol==tf] <- nrow(atac_tagMatrix)
      tf_metadata$number_of_Target_Genes_window_name[tf_metadata$Target.gene.symbol==tf] <- paste0(target_transcripts_coord$hgnc_symbol, collapse=",")
     }
   }
}
tf_metadata$Gene <- sapply(tf_metadata$Gene, paste, collapse = ";")
write.table(tf_metadata, "tableS5c.txt", sep="\t", quote=F)
