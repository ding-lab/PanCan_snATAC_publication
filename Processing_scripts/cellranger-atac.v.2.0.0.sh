#Basic commands to run cellranger-atac for single-nuclei ATAC-sequencing Data

#References: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/what-is-cell-ranger-atac

#Reference Genome:https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/downloads/latest?

#Sample_ID
sample=''

#Path to reference
reference='refdata-cellranger-arc-GRCh38-2020-A-2.0.0'

#Path to FASTQ files
fastq_dir='FASTQ/'

cellranger-atac-2.0.0/cellranger-atac count --id $sample --fastqs $fastq_dir --reference $reference --localcores=12 --localmem=100
