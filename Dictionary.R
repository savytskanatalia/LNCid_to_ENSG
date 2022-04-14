library(data.table)
library(stringr)
all_genes_ann<-fread( 'gencode_v29_lncpedia_GenesOnly_ENSG.gtf' , select=c(9), sep="\t")
all_genes_ann <- unique(all_genes_ann$V9[all_genes_ann$V9 %like% "ENSG"])
all_genes_ann<-str_split_fixed(all_genes_ann, "; ",n=12)
all_genes_ann<-as.data.frame(all_genes_ann)
all_genes_ann$V1<-str_remove(str_remove(all_genes_ann$V1, "gene_id \""),"\"")
all_genes_ann$V3<-str_remove(str_remove(all_genes_ann$V3, "gene_alias_[0-9] \""),"\"")
all_genes_ann$V4<-str_remove(str_remove(all_genes_ann$V4, "gene_alias_[0-9] \""),"\"")
all_genes_ann$V5<-str_remove(str_remove(all_genes_ann$V5, "gene_alias_[0-9] \""),"\"")
all_genes_ann$V6<-str_remove(str_remove(all_genes_ann$V6, "gene_alias_[0-9] \""),"\"")
all_genes_ann$V7<-str_remove(str_remove(all_genes_ann$V7, "gene_alias_[0-9] \""),"\"")
all_genes_ann$V8<-str_remove(str_remove(all_genes_ann$V8, "gene_alias_[0-9] \""),"\"")
all_genes_ann$V9<-str_remove(str_remove(all_genes_ann$V9, "gene_alias_[0-9] \""),"\"")

all_genes_ann$V10<-str_remove(str_remove(all_genes_ann$V10, "gene_alias_[0-9] \""),"\"")
all_genes_ann$V11<-str_remove(str_remove(all_genes_ann$V11, "gene_alias_[0-9] \""),"\"")
all_genes_ann$V12<-str_remove(str_remove(all_genes_ann$V12, "gene_alias_[0-9] \""),"\"")
all_genes_ann<-all_genes_ann[,-2]
# For each column except first replace value with "-" each time there is no ENSG match, then merge all columns except for 1
for (j in 1:nrow(all_genes_ann)){
  for (i in 2:ncol(all_genes_ann)){
    if ((all_genes_ann[j,i] %like% "ENSG")==FALSE){
      all_genes_ann[j,i]<-"-"
    } else {
      if ((all_genes_ann[j,i-1] %like% "ENSG")==TRUE){
        all_genes_ann[j,i]<-"-"
      } else {next}
      
    }
}
}

for (i in 1:nrow(all_genes_ann)){
  if ((all_genes_ann[i,2]%like% "ENSG")==TRUE){
    next
  } else {
    all_genes_ann[i,2]<-paste0(all_genes_ann[i,3],paste0(all_genes_ann[i,4],paste0(all_genes_ann[i,5],paste0(all_genes_ann[i,6],paste0(all_genes_ann[i,7],paste0(all_genes_ann[i,8],paste0(all_genes_ann[i,9],paste0(all_genes_ann[i,10],paste0(all_genes_ann[i,11],all_genes_ann[i,12])))))))))
  }
}
all_genes_ann$V3<-str_remove(all_genes_ann$V3,"\\-\\-\\-")
all_genes_ann$V3<-str_remove(all_genes_ann$V3,"\\-\\-\\-")
all_genes_ann$V3<-str_remove(all_genes_ann$V3,"\\-\\-")
all_genes_ann$V3<-str_remove(all_genes_ann$V3,"\\-")
# gsub("\\..*","",all_genes)
all_genes_ann[all_genes_ann$V3 %like% "transcript_alias","V3"]<-str_remove(str_remove(gsub("\\|.*","",all_genes_ann[all_genes_ann$V3 %like% "transcript_alias","V3"]),"transcript_alias_[0-9] "),"\"")
all_genes_ann[all_genes_ann$V3 %like% "ENST","V3"]<-gsub("\\-ENST.*","",all_genes_ann[all_genes_ann$V3 %like% "ENST","V3"])
all_genes_ann[all_genes_ann$V3 %like% "ENST","V3"]<-gsub("transcript_alias.*","",all_genes_ann[all_genes_ann$V3 %like% "ENST","V3"])
all_genes_ann[all_genes_ann$V3 %like% "ENST","V3"]<-gsub("ENST.*; ","",all_genes_ann[all_genes_ann$V3 %like% "ENST","V3"])
all_genes_ann[all_genes_ann$V3 %like% "ENSG.*ENSG","V3"]<-gsub("ENSG.*ENSG","ENSG",all_genes_ann[all_genes_ann$V3 %like% "ENSG.*ENSG","V3"])
all_genes_ann$V3<-gsub("\\..*","",all_genes_ann$V3)
# DLEU1 = ENSG00000176124
all_genes_ann[all_genes_ann$V1=="DLEU1",2]<-"ENSG00000176124"
all_genes_ann<-all_genes_ann[,1:2]
# Now dictionary for elements where V1 is not ENSG ID
colnames(all_genes_ann)<-c("ID","ENSEMBL")
# Add also IDs missing
all_genes_ann<-all_genes_ann[!(all_genes_ann$ID %like% "ENSG"),]
add_df<-data.frame(ID=c("HCG23","LINC00678","lnc-PIWIL4-1","lnc-SNRPN-1","P3H2-AS1","LNCPRESS2","LINC01950","SNHG4","lnc-STK32A-3","PVT1","lnc-PRAG1-3","lnc-DIAPH2-16",
                        "OLMALINC","lnc-SLC38A2-1","lnc-PCED1B-5","GAU1","lnc-PCDH8-10","G2E3-AS1","lnc-FAM103A1-2","lnc-SHISA9-3","lnc-CDH5-3","lnc-DHX38-3","lnc-VSTM2B-5",
                        "lnc-A1BG-1","lnc-MYO7B-2","lnc-RBM11-2","LINC01424","LINC01487","LINC02211","CKMT2-AS1","CASC15","lnc-KLHL32-4","lnc-ARL4A-54","lnc-KCNB2-12","FIRRE"),
                   ENSEMBL=c("ENSG00000228962","ENSG00000254934","ENSG00000255666","ENSG00000279192","ENSG00000225764","ENSG00000249152","ENSG00000251027","ENSG00000281398","ENSG00000280780",
                             "ENSG00000249859","ENSG00000253893"," ENSG00000281566","ENSG00000235823","ENSG00000257496","ENSG00000272369","ENSG00000255474","ENSG00000274090",
                             "ENSG00000257636","ENSG00000260608","ENSG00000262801","ENSG00000260364","ENSG00000260664","ENSG00000267498","ENSG00000268516","ENSG00000231731",
                             "ENSG00000224905","ENSG00000236519","ENSG00000241336","ENSG00000245662","ENSG00000247572","ENSG00000272168","ENSG00000271860","ENSG00000229618",
                             "ENSG00000254277","ENSG00000213468"))

all_genes_ann<-rbind(all_genes_ann,add_df)

fwrite(all_genes_ann, 'LNCID_to_ENSG_dictionary.txt',sep="\t",quote = F)
