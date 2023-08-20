library(SCENIC)
library(SCopeLoomR)
library(data.table)


#Calcultate differential regulons for each cancer vs all others, using tumor cells only

disease='Pancan'
wdir=''
path_to_ATAC_catalog=''



#Read regulons scores

mat=fread('out/UpdatedByFreq.0.8_regulons_CellsAUC.20230124.tsv')
mat=as.data.frame(mat)
rownames(mat)=mat[,1]
mat=mat[,-1]



#Read catalog

cat=read_delim(path_to_ATAC_catalog,delim='\t')
cat=as.data.frame(cat)
colnames(cat)=gsub(' ','_',colnames(cat))
cat$Disease_Type=ifelse(cat$Disease_Type=='PKD','ccRCC',cat$Disease_Type)

cat$Piece_ID=paste(cat$Disease_Type,cat$Piece_ID,sep='_')
cat=cat %>% dplyr::select ('Piece_ID','Sample_Type')
colnames(cat)[2]=c('Sample_type')

res=mat



#Read loom object cell type annotation

info=read.table('../../../Pancan_obj/cellInfo_1_PancanObj_200cellsPerPieceID_SelCellTypes.tsv',sep='\t',
header=T)
rownames(info)=info$Barcode
info=info[colnames(res),]
info$Disease=gsub('(.*)_(.*)_(.*)','\\1',info$Barcode)
info$Disease=gsub('(.*)_(.*)','\\1',info$Disease)
info$Piece_ID_1=paste(info$Disease,info$Piece_ID,sep='_')

info$Disease=ifelse(info$Piece_ID_1 %in% c("BRCA_HT268B1-Th1H3", "BRCA_HT029B1-S1PC", "BRCA_HT035B1-S1PA",
 "BRCA_HT1408-06","BRCA_HT141B1-S1H1", "BRCA_HT206B1-S1H4", "BRCA_HT271B1-S1H3",
"BRCA_HT378B1-S1H1", "BRCA_HT378B1-S1H2", "BRCA_HT384B1-S1H1", "BRCA_HT517B1-S1H1"),'BRCA_Basal', info$Disease)



# Subset tumor cells

info_s=info[info$Cell_type=='Tumor' & info$Piece_ID_1 %in% cat$Piece_ID[cat$Sample_type=='Tumor'] &
!(info$Piece_ID_1 %in% c('CESC_CE332E1-N1','CESC_CE336E1-S1','CESC_CE354E1-S1','CESC_CE357E1-S1','CESC_CE507-C1A2','PDAC_PM565P1-T1N1')),]
res_s=res[,info_s$Barcode]

cell_types=unique(info_s$Disease)



#Run differential analysis

final_wilcoxon_stat=NULL
final_wilcoxon_stat_all=NULL
for (cell_t1 in cell_types){
    print(cell_t1)
    res_1=res_s[,colnames(res_s) %in% rownames(info_s)[info_s$Disease==cell_t1]]
    res_2=res_s[,colnames(res_s) %in% rownames(info_s)[info_s$Disease!=cell_t1]]
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
final_wilcoxon_stat=rbind(final_wilcoxon_stat,all_wilcoxon_stat)
}

final_wilcoxon_stat=as.data.frame(final_wilcoxon_stat)
final_wilcoxon_stat$cell_t2='Other'
write.table(final_wilcoxon_stat,'out/Regulons_AUC_difference.eachCancer_vs_Others.tsv',sep='\t',quote=F,row.names=F)
