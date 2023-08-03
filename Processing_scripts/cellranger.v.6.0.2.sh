#Basic commands to run cellranger for single-cell or single-nuclei RNA-sequencing Data

#References: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger

#Reference Genome: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

#Sample_ID
sample=''

#Path to reference
reference='refdata-gex-GRCh38-2020-A'

#Path to FASTQ files
fastq_dir='FASTQ/'


cellranger-6.0.2/cellranger count --id $sample --fastqs $fastq_dir --include-introns --localmem=120 --localcores=16 --transcriptome $reference
