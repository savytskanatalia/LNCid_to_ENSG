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
fwrite(all_genes_ann, 'LNCID_to_ENSG_dictionary.txt',sep="\t",quote = F)

