library(Signac)
library(Seurat)
library(Matrix)
library(AUCell)
library(SCENIC)
library(SCopeLoomR)
library(GSEABase)


disease='Pancan'
wdir=''

regulons=vector(length=10,mode = "list")

#Collect regulons data from all 10 iterations

for (iter in 1:10){
loom=open_loom(paste(wdir,'/',disease,'_obj/iter_',iter,'/',disease,'_pyscenic_output.loom',sep=''),mode='r')

# Read information from loom file:
regulons_incidMat <- get_regulons(loom,column.attr.name = "Regulons")

regulons[[iter]] <- regulonsToGeneLists(regulons_incidMat)
print(iter)
}

all_tfs=NULL
for (iter in 1:10){
    all_tfs=c(all_tfs,names(regulons[[iter]]))
}
stat=table(all_tfs)
stat=as.data.frame(stat)

saveRDS(regulons, "TF_taargets_across10Iterations.rds")
write.table(stat, 'TF_frequency.10Iter.tsv',sep='\t',quote=F,row.names=F)

#Only keep TFs that were present in at least 8 iterations

sel_tfs=stat$all_tfs[stat$Freq>=8]

all_st=NULL
for (tf in sel_tfs){
    all_gn=NULL
    for (iter in 1:10){
    	if(tf %in% names(regulons[[iter]])){
	        g_n=regulons[[iter]][names(regulons[[iter]])==tf][[1]]
		all_gn=c(all_gn,g_n)
		}
	}
	st=as.data.frame(table(all_gn))
	st$TF=tf
	all_st=rbind(all_st,st)
	print(tf)
}
all_st=as.data.frame(all_st)
colnames(all_st)[1]=c('Gene')

write.table(all_st, "Target_frequency_byTF.10Iter.tsv",sep='\t',row.names=F,quote=F)