setwd("/ifs/scratch/c2b2/ac_lab/npg2108/HMMicro/")

library(DeepBlueR)
library(GenomicRanges)

cat("loading data...\n")
load("/ifs/scratch/c2b2/ac_lab/npg2108/HMMicro/RDAs/DeepBlueR_chr22_score_matrix.RDa")


lapply(score_matrix_chr22,colnames)


####
cat("replacing NA for 0...\n")
sub_score_matrix_chr22<-NULL
for(i in 1:length(score_matrix_chr22)){
  sm<-score_matrix_chr22[[i]]
  sub_score_matrix_chr22[[i]]<-sm
  nonvalinds<-which(is.na(sm[,4]))
  sub_score_matrix_chr22[[i]][nonvalinds,4]<-0
}

####
cat("converting data table to data frame class...\n")
sub_df_score_matrix_chr22<-lapply(sub_score_matrix_chr22,as.data.frame)

####
cat("make original dfs into signals col List...\n")
cols_chr22<-lapply(sub_df_score_matrix_chr22,function(x){x[,4]})

####
cat("combining list of columns into final matrix...\n")
mat_chr22<-do.call(cbind,cols_chr22)
dim(mat_chr22)

####
#noting ChIP-Seq TFs in which the data represent and adding as mat column names.
#Had to look up on encodeproject.org
#should use later on.
cat("extra processing...\n") 
encode<-c("Kap1","Pol2ra","TCFL2","Control.1","Control.2","ZNF263",
          "H3K4me3.1","H3K4me3.2","Control.3","CTCF.1","ELK4","CTCF.2",
          "Control.4","DNase")

####
#changing column names, 
#merging replicates and removing controls
newmat_chr22<-mat_chr22
colnames(newmat_chr22)<-encode
newmat_chr22<-newmat_chr22[,-c(4,5,9,13)]
H3K4me3_chr22<-apply(newmat_chr22[,c(5,6)],1,mean)
CTCF_chr22<-apply(newmat_chr22[,c(7,9)],1,mean)
newmat_chr22<-cbind(newmat_chr22[,c(1:4,8,12)],CTCF_chr22,H3K4me3_chr22)

####
cat("output final matrix...\n")
write.table(newmat_chr22,file="/ifs/scratch/c2b2/ac_lab/npg2108/HMMicro/data/final_matrix_chr22.txt",sep="\t",row.names=F,col.names=T)

####
cat("saving image...\n")
save.image("/ifs/scratch/c2b2/ac_lab/npg2108/HMMicro/data/DeepBlueR_chr22.RDa")
