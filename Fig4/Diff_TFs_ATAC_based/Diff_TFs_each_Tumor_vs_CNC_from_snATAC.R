library(Signac)
library(Seurat)
library(reshape)
library(reshape2)
library(plyr)
library(dplyr)


data_dir=''
path_to_atac_pancan_obj=''
path_to_ATAC_cell_type_annotation=''
path_to_ATAC_catalog=''
path_to_TFs_annot=''


#Read snATAC pan-cancer object

atac=readRDS(path_to_atac_pancan_obj)



#Read snATAC cell type annotation. Add annotation to snATAC object.

meta=read.table(path_to_ATAC_cell_type_annotation, sep='\t',header=T)

meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)

orig_1=as.data.frame(atac$Piece_ID)
rownames(meta)=meta$X
meta=meta[rownames(orig_1),]
atac$cell_type.detailed=meta$cell_type.detailed



#Read ATAC catalog

samples=read.table(path_to_ATAC_catalog,sep='\t',header=T)
s_sampl=samples[,c('Piece_ID','Disease.Type','Sample.Type')]
s_sampl$Piece_ID_2=paste(s_sampl$Disease.Type,s_sampl$Piece_ID,sep='_')
met_samples=s_sampl[s_sampl$Sample.Typ=='Met',]



#Annotate BRCA basal subtype samples

atac$Disease=ifelse(atac$Piece_ID %in% c("BRCA_HT268B1-Th1H3", "BRCA_HT029B1-S1PC", "BRCA_HT035B1-S1PA",
 "BRCA_HT1408-06","BRCA_HT141B1-S1H1", "BRCA_HT206B1-S1H4", "BRCA_HT271B1-S1H3",
"BRCA_HT378B1-S1H1", "BRCA_HT378B1-S1H2", "BRCA_HT384B1-S1H1", "BRCA_HT517B1-S1H1"), "BRCA_Basal",
atac$Disease)



#Set idents for cells; they will be used in the differential analysis

atac$test=ifelse(atac$cell_type.detailed=='Tumor' & !(atac$Piece_ID %in% c(met_samples$Piece_ID_2,'CESC_CE332E1-N1','CESC_CE336E1-S1',
'CESC_CE354E1-S1','CESC_CE357E1-S1','CESC_CE507-C1A2','PDAC_PM565P1-T1N1')),
paste(atac$cell_type.detailed,atac$Disease,sep='__'),atac$cell_type.detailed)


ATAC=atac



#Extract table with chromvar TF scores

cell_types=unique(ATAC$test)
DefaultAssay(ATAC) <- 'chromvar'

chromv= GetAssayData(object = ATAC)



#Read annotation for TFs, and add TF names

jaspar=read.table(path_to_TFs_annot, sep='\t', header=TRUE)

mtx0=chromv
res=merge(mtx0,jaspar,by=0,all.x=TRUE)
rownames(res)=res$motif.name
res=res[,-1]
res=res[,1:(ncol(res)-2)]



#Extract annotation of cells, that will be used in the differential analysis

ann_col0 = data.frame('cell_types' = ATAC$test)
ann_col0$cel_2=ann_col0$cell_types
ann_col0=ann_col0[order(ann_col0$cell_types),]
ann_col1=data.frame("cell_types"=ann_col0$cell_types)
rownames(ann_col1)=rownames(ann_col0)



#Set the same order in the TF score matrix

res=res[,rownames(ann_col0)]



#Run differential analysis: tumor cells vs its CNC (which can be found in Fig. 2a)

final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL


cancers=c('ccRCC','PDAC','UCEC','CRC','BRCA','BRCA_Basal','GBM','MM','OV','SKCM','HNSCC','CESC')
for (can in cancers){
    if (can %in% c('ccRCC')){
    origin_c=c('Proximal Tubule')
}else if (can=='BRCA'){
    origin_c=c("Luminal mature")
}else if (can=='BRCA_Basal'){
    origin_c=c('Luminal progenitor')
}else if (can=='GBM'){
    origin_c=c('OPC')
}else if (can=='PDAC'){
    origin_c=c('Ductal-like2')
}else if (can=='UCEC'){
    origin_c=c('Secretory Endometrial epithelial cells')
}else if (can=='CRC'){
    origin_c=c('Distal Stem Cells')
}else if (can=='MM'){
    origin_c=c('B-cells')
}else if (can=='OV'){
    origin_c=c('Secretory Endometrial epithelial cells')
}else if (can=='HNSCC'){
    origin_c=c('Normal squamous cells')
}else if (can=='CESC'){
    origin_c=c('Normal squamous cells')
}else if (can=='SKCM'){
    origin_c=c('Melanocytes')
}  

print (can)
    cell_t1=paste('Tumor',can,sep='__')  
    for (cell_t2 in origin_c){
    if(length(rownames(ann_col1)[ann_col1$cell_types==cell_t1])>=50 & 
length(rownames(ann_col1)[ann_col1$cell_types==cell_t2])>=50)
    res_1=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t1]]
    res_2=res[,colnames(res) %in% rownames(ann_col1)[ann_col1$cell_types==cell_t2]]
    all_wilcoxon_stat=NULL
        for (motif in 1:nrow(res)){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[motif,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[motif,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[motif,]))),
as.numeric(as.character(unlist(res_2[motif,]))))
           stat=cbind(cell_t1,rownames(res)[motif],mean_score1,mean_score2,w_test$p.value)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]

colnames(all_wilcoxon_stat)[1:2]=c('cell_t1','TF')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
all_wilcoxon_stat$FDR=p.adjust(all_wilcoxon_stat$pvalue, method="fdr")
all_wilcoxon_stat$cell_t2=cell_t2
all_wilcoxon_stat$Disease=can

final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
print(cell_t2)
}
print(can)
}


final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
colnames(final_wilcoxon_stat)[1:2]=c('cell_t1','TF_Name')

write.table(final_wilcoxon_stat, "out/TF_score_difference.Tumor_vs_CNC_ATAC_based.tsv", quote=FALSE,sep="\t",row.names=FALSE)