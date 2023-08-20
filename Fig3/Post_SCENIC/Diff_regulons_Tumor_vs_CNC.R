library(Matrix)
library(data.table)
library(tidyverse)

path_to_ATAC_catalog=''
path_to_RNA_cell_type_annotation=''



#Read regulons scores

res=fread('UpdatedByFreq.0.8_regulons_CellsAUC.20230124.tsv',sep=''))
res=as.data.frame(res)
rownames(res)=res[,1]
res=res[,-1]



#List of all CNCs

sel_cell_o=c('Secretory Endometrial epithelial cells','Proximal Tubule','Ductal-like2','Distal Stem Cells',
'OPC','Luminal progenitor','Luminal mature','B-cells','Normal squamous cells','Melanocytes')



#Read catalog

cat=read_delim(path_to_ATAC_catalog,delim='\t')
cat=as.data.frame(cat)
colnames(cat)=gsub(' ','_',colnames(cat))
cat$Disease_Type=ifelse(cat$Disease_Type=='PKD','ccRCC',cat$Disease_Type)
 
cat$Piece_ID=paste(cat$Disease_Type,cat$Piece_ID,sep='_')
cat=cat %>% dplyr::select ('Piece_ID','Sample_Type')
colnames(cat)[2]=c('Sample_type')



#Read RNA cell type annotation

meta=read.table(path_to_RNA_cell_type_annotation,sep='\t',header=T)
meta$cell_type.detailed=ifelse(!is.na(meta$cell_type.normal),meta$cell_type.normal,
meta$cell_type.harmonized.cancer)
meta_s=meta[,c('X','cell_type.detailed')]
colnames(meta_s)[1]='Barcode'



#Read loom object cell type annotation

info=read.table('../../Pancan_obj/cellInfo_1_PancanObj_200cellsPerPieceID_SelCellTypes.20230121.tsv',sep='\t',
header=T)
rownames(info)=info$Barcode
info=info[colnames(res),]
info$Disease=gsub('(.*)_(.*)_(.*)','\\1',info$Barcode)
info$Disease=gsub('(.*)_(.*)','\\1',info$Disease)
info$Piece_ID_1=paste(info$Disease,info$Piece_ID,sep='_')



#Identify BRCA basal samples

info$Disease=ifelse(info$Piece_ID_1 %in% c("BRCA_HT268B1-Th1H3", "BRCA_HT029B1-S1PC", "BRCA_HT035B1-S1PA",
 "BRCA_HT1408-06","BRCA_HT141B1-S1H1", "BRCA_HT206B1-S1H4", "BRCA_HT271B1-S1H3",
"BRCA_HT378B1-S1H1", "BRCA_HT378B1-S1H2", "BRCA_HT384B1-S1H1", "BRCA_HT517B1-S1H1"),'BRCA_Basal', info$Disease)

info_2=merge(info,meta_s,all.x=T)
info_2$cell_type.detailed=ifelse(info_2$cell_type.detailed %in% c('Proximal Tubule convoluted',
'Proximal Tubule'),'Proximal Tubule',info_2$cell_type.detailed)



# Subset tumor cells from primary tumor samples and from CNC. Also exclude CEAD samples

info_s=info_2[(info_2$cell_type.detailed %in% c(sel_cell_o)) | (info_2$Piece_ID_1 %in% cat$Piece_ID[cat$Sample_type=='Tumor'] & 
!(info_2$Piece_ID_1 %in% c('CESC_CE332E1-N1','CESC_CE336E1-S1','CESC_CE354E1-S1','CESC_CE357E1-S1','CESC_CE507-C1A2','PDAC_PM565P1-T1N1')) & 
info_2$cell_type.detailed=='Tumor'),]

rownames(info_s)=info_s$Barcode

res_s=res[,info_s$Barcode]



#Set idents for cells that will be used in the differential analysis

info_s$test=ifelse(info_s$cell_type.detailed=='Tumor',paste('Tumor',info_s$Disease,sep='_'),
info_s$cell_type.detailed)



#Run differential analysis

final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL
cancers=c('ccRCC','PDAC','UCEC','CRC','BRCA','BRCA_Basal','GBM','OV','MM','HNSCC','CESC','SKCM')
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
}else if (can=='OV'){
    origin_c=c('Secretory Endometrial epithelial cells')
}else if (can=='MM'){
    origin_c=c('B-cells')
}else if (can=='SKCM'){
    origin_c=c('Melanocytes')
}else if (can %in% c('CESC','HNSCC')){
    origin_c=c('Normal squamous cells')
}

print (can)
cell_t1=paste('Tumor',can,sep='_')
cell_t2=origin_c
print(paste(cell_t1,cell_t2, sep=' '))

    res_1=res_s[,colnames(res_s) %in% rownames(info_s)[info_s$test==cell_t1]]
    res_2=res_s[,colnames(res_s) %in% rownames(info_s)[info_s$test==cell_t2]]
    all_wilcoxon_stat=NULL
    for (index in 1:nrow(res_s)){
           mean_score1=mean(as.numeric(as.character(unlist(res_1[index,]))),na.rm=TRUE)
           mean_score2=mean(as.numeric(as.character(unlist(res_2[index,]))),na.rm=TRUE)
           w_test=wilcox.test(as.numeric(as.character(unlist(res_1[index,]))),
as.numeric(as.character(unlist(res_2[index,]))))
           stat=cbind(cell_t1,rownames(res_s)[index],mean_score1,mean_score2,w_test$p.value)
           all_wilcoxon_stat=rbind(all_wilcoxon_stat,stat)
        }
        all_wilcoxon_stat=as.data.frame(all_wilcoxon_stat)
        all_wilcoxon_stat$V5=as.numeric(as.character(unlist(all_wilcoxon_stat$V5)))
        all_wilcoxon_stat=all_wilcoxon_stat[order(all_wilcoxon_stat$V5),]
        all_wilcoxon_stat$mean_score1=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score1)))
        all_wilcoxon_stat$mean_score2=as.numeric(as.character(unlist(all_wilcoxon_stat$mean_score2)))
        all_wilcoxon_stat$diff=all_wilcoxon_stat$mean_score1-all_wilcoxon_stat$mean_score2
        all_wilcoxon_stat=all_wilcoxon_stat[order(-all_wilcoxon_stat$diff),]

colnames(all_wilcoxon_stat)[1:2]=c('cell_t1','Regulon')
colnames(all_wilcoxon_stat)[5]=c('pvalue')
all_wilcoxon_stat$FDR=p.adjust(all_wilcoxon_stat$pvalue, method="fdr")
all_wilcoxon_stat$cell_t2=cell_t2
all_wilcoxon_stat$Disease=can
final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
final_wilcoxon_stat$cell_t2='Other'
write.table(final_wilcoxon_stat,'out/Regulons_AUC_difference.tumor_vs_CNC.tsv',sep='\t',
quote=F,row.names=F)



################################################
####Also select 200 cells per group and save:###
################################################

#These are used for plotting (Fig. 4b)
info_s$Disease_1=ifelse(info_s$Disease %in% c('BRCA','BRCA_Basal') & info_s$cell_type.detailed!='Tumor','BRCA',info_s$Disease)
info_s$ID=paste(info_s$Disease_1,info_s$cell_type.detailed,sep='__')

n_cells=200
all_cells_s=NULL
for (id in unique(info_s$ID)){
    info_1=info_s[info_s$ID==id,]
    cells_s=sample(rownames(info_1),min(n_cells,nrow(info_1)))
    all_cells_s=c(all_cells_s,cells_s)
}

res_s=res_s[,all_cells_s]
info_s_sel=info_s[info_s$Barcode %in% colnames(res_s),]
write.table(res_s, 'out/AUC_Pancan_200_RandomTumorCellsPerCancer_withNormalCells.tsv',sep='\t',quote=F)
write.table(info_s_sel, 'out/Annot_cells_AUC_Pancan_200_RandomTumorCellsPerCancer_withNormalCells.',sep='\t',quote=F)