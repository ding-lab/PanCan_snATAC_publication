library(Signac)
library(Seurat)
library(Matrix)
library(SCENIC)
library(SCopeLoomR)
library(GSEABase)
library(AUCell)

###Import from the loom file

disease='Pancan'
wdir=''
path_to_ATAC_catalog=''
loom=open_loom(paste(wdir,'/',disease,'_obj/iter_1/',disease,'_pyscenic_output.loom',sep=''),mode='r')



#Only keep regulons and their targets if they were found in at least 8 iterations.

tfs=read.table('../TF_frequency.10Iter.tsv',sep='\t',header=T)
tfs_s=tfs[tfs$Freq>=8,]
targ=read.table('../Target_frequency_byTF.10Iter.tsv',sep='\t',header=T)
targ_s=targ[targ$Freq>=8,]

regulons_new=vector(mode='list',length=nrow(tfs_s))
for (i in 1: length(tfs_s$all_tfs)){
    tf=tfs_s$all_tfs[i]
    genes=targ_s$Gene[targ_s$TF==tf]
    regulons_new[[i]]=genes
    names(regulons_new)[[i]]=tf
}
exprMatrix=get_dgem(loom)



#Calculate AUC scores for the new set of regulons

cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=F, nCores=40)

cells_AUC=AUCell_calcAUC(
  regulons_new,
  cells_rankings,
  nCores = 40,
  normAUC = TRUE,
  aucMaxRank = ceiling(0.05 * nrow(cells_rankings)),
  verbose = TRUE
)
mat=getAUC(cells_AUC)
write.table(mat, "out/UpdatedByFreq.0.8_regulons_CellsAUC.20230124.tsv",sep='\t',quote=F,row.names=T)

saveRDS(regulons_new, "out/UpdatedByFreq.0.8_regulons.20220124.rds")



#Annotate regulons_new:

all_st=NULL
for (i in 1:length(regulons_new)){
    reg=names(regulons_new)[[i]]
    genes_n=length(regulons_new[[i]])
    st=cbind(reg,genes_n)
    all_st=rbind(all_st,st)
}
all_st=as.data.frame(all_st)
colnames(all_st)=c('Regulon','Genes_N')
write.table(all_st,'out/Regulons_new_annot.20230124.tsv',sep='\t',quote=F,row.names=F)
