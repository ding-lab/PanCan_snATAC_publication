#Basic commands to run cellranger for single-nuclei Multiome-sequencing Data

#References: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/what-is-cell-ranger-arc

#Reference Genome:https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/downloads/latest?

#Sample_ID
sample=''

#Path to reference
reference='refdata-cellranger-arc-GRCh38-2020-A-2.0.0'

#Path to libraries
libraries='library.csv'

cellranger-arc-2.0.0/cellranger-arc count --id $sample --reference  $reference --libraries $libraries --disable-ui --localcores 12 --localmem 100
