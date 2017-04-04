####################
#Getting HEK293 epigenomic data from
#ENCODE via the DeepBlueR R package.
#
#DeepBlueR allows downloading of curated
#datasets from 4 epigenomics consortia
#
#https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkw211
#
#Getting HEK293 epigenetic signals in
#genomic intervals across the whole genome
#
#Output is a matrix with rows as 1000 bp-wide genomic windows with
#signal for every epigenetic mark in HEK293 cells and
#columns for each epigenetic mark assayed in HEK293 cells
####################

#####Install/load libraries####
source("https://bioconductor.org/biocLite.R")
biocLite("DeepBlueR")
library(DeepBlueR)
library(GenomicRanges)
biocLite("rtracklayer")
library(rtracklayer)
library(gplots)
library(RColorBrewer)
library(matrixStats)
library(stringr)
#####Making and outputting signal matrix:-run interactively#####
####
setwd("~/GitHub/HMMicro/")
####
#list experiments via DeepBlueR server
#projects<-deepblue_list_projects()[["name"]]
experiments = deepblue_list_experiments(type="signal", genome = "hg19",
                                        biosource="HEK293",
                                        project="ENCODE")

####
#get experiments and put into proper format
exps <- list(nrow(experiments))
exp_columns<-NULL
for(i in 1:nrow(experiments)){
  exp_columns[[i]] <- "VALUE"
}
names(exp_columns) <- deepblue_extract_names(experiments)
#filter for repeat file names
exp_columns <- deepblue_select_column(experiments, "VALUE")
#####chr1#####
####
#making genomic windows
tiling_regions_chr1 = deepblue_tiling_regions(size=1, genome="hg19",chromosome="chr1")

####
#requesting signals-in-windows job in DeepBlueR server
request_id_chr1<-NULL
for(i in 1:length(exp_columns)){
  #submitting requesting signals in windows job
  request_id_chr1[[i]] <- deepblue_score_matrix(
    experiments_columns = exp_columns[i],
    aggregation_function = "mean",
    aggregation_regions_id = tiling_regions_chr1)
}

####
#checking example job
deepblue_info(request_id_chr1[[6]])$state

####
#dowmloading score matrix once jobs are done
#as indicated above
score_matrix_chr1<-NULL
for(i in 1:length(request_id_chr1)){
  if(deepblue_info(request_id_chr1[[i]])$state=="done"){
    #downloading matrix 
    score_matrix_chr1[[i]] <- deepblue_download_request_data(request_id = request_id_chr1[[i]])
  }
}

####
#removing NA regions
sub_score_matrix_chr1<-NULL
for(i in 1:length(score_matrix_chr1)){
  sm<-score_matrix_chr1[[i]]
  valinds<-which(!is.na(sm[,4]))
  sub_score_matrix_chr1[[i]]<-sm[valinds,]
}

####
#converting data table to data frame class
sub_df_score_matrix_chr1<-lapply(sub_score_matrix_chr1,as.data.frame)

####
#output beds to directory
for(i in 1:length(sub_df_score_matrix_chr1)){
  write.table(sub_df_score_matrix_chr1[[i]],
              paste0("~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/bed_chr1/",
                     colnames(sub_df_score_matrix_chr1[[i]])[4],".bed"),
              quote=F,col.names=F,row.names=F,sep="\t")
}

####
#intersect with bedops-use in terminal
#download bedops from https://github.com/bedops/bedops
#paths conformed for use on local
#bedops --intersect ~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/bed_chr1/* > ~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/intersected_regions_chr1.bed

####
#upload from bedops and store as GRanges
#paths conformed to use on local
intersected_regions_chr1<-read.table("~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/intersected_regions_chr1.bed",header=F,sep="\t")
colnames(intersected_regions_chr1)<-c("chrom","start","end")
GRanges_intersected_regions_chr1<-makeGRangesFromDataFrame(intersected_regions_chr1)
#reading in bed files if needed
dir="~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/bed_chr1"
files<-list.files(path = dir,pattern=".bed$")
sub_df_score_matrix_chr1<-NULL
for(i in 1:length(files)){
  f<-paste0(dir,"/",files[i])
  sub_df_score_matrix_chr1[[i]]<-read.table(f,header=F)
}
for(i in 1:length(sub_df_score_matrix_chr1)){
  colnames(sub_df_score_matrix_chr1[[i]])<-c("chr","start","end","signal")
}

####
#make original dfs into GRanges List
ListOfGRanges_chr1<-lapply(sub_df_score_matrix_chr1,function(x){
  makeGRangesFromDataFrame(x,keep.extra.columns = T)})

####
#getting conserved columns across all GRanges
cols_chr1<-NULL
for(i in 1:length(ListOfGRanges_chr1)){
  target<-ListOfGRanges_chr1[[i]]
  conserved<-GRanges_intersected_regions_chr1
  meta<-mcols(target[which(conserved %over% target),])
  cols_chr1[[i]]<-as.data.frame(meta)
}

####
#combining list of columns into final matrix
mat_chr1<-do.call(cbind,cols_chr1)

####
#noting ChIP-Seq TFs in which the data represent and adding as mat column names.
#Had to look up on encodeproject.org
#should use later on. 
colnames(mat_chr1)
encode<-c("Kap1","Pol2ra","TCFL2","Control.1","Control.2","ZNF263",
          "H3K4me3.1","H3K4me3.2","Control.3","CTCF.1","ELK4","CTCF.2",
          "Control.4")

####
#changing column names, 
#merging replicates and removing controls
newmat_chr1<-mat_chr1
colnames(newmat_chr1)<-encode
newmat_chr1<-newmat_chr1[,-c(4,5,9,13)]
H3K4me3_chr1<-apply(newmat_chr1[,c(5,6)],1,mean)
CTCF_chr1<-apply(newmat_chr1[,c(7,9)],1,mean)
newmat_chr1<-cbind(newmat_chr1[,c(1:4,8)],CTCF_chr1,H3K4me3_chr1)

####
#output final matrix
write.table(newmat_chr1,file="~/GitHub/HMMicro/data/final_matrix_chr1.txt",sep="\t",row.names=F,col.names=T)

####
#saving environment
save.image("~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/DeepBlueR_chr1.RDa")
#####chr22#####
####
#making genomic windows
tiling_regions_chr22 = deepblue_tiling_regions(size=1, genome="hg19",chromosome="chr22")

####
#requesting signals-in-windows job in DeepBlueR server
request_id_chr22<-NULL
for(i in 1:length(exp_columns)){
  #submitting requesting signals in windows job
  request_id_chr22[[i]] <- deepblue_score_matrix(
    experiments_columns = exp_columns[i],
    aggregation_function = "mean",
    aggregation_regions_id = tiling_regions_chr22)
}

####
#checking example job
deepblue_info(request_id_chr22[[9]])$state

####
#dowmloading score matrix once jobs are done
#as indicated above
score_matrix_chr22<-NULL
for(i in 1:length(request_id)){
  if(deepblue_info(request_id_chr22[[i]])$state=="done"){
    #downloading matrix 
    score_matrix_chr22[[i]] <- deepblue_download_request_data(request_id = request_id_chr22[[i]])
  }
}

####
#removing NA regions
sub_score_matrix_chr22<-NULL
for(i in 1:length(score_matrix_chr22)){
  sm<-score_matrix_chr22[[i]]
  valinds<-which(!is.na(sm[,4]))
  sub_score_matrix_chr22[[i]]<-sm[valinds,]
}

####
#converting data table to data frame class
sub_df_score_matrix_chr22<-lapply(sub_score_matrix_chr22,as.data.frame)

####
#output beds to directory
for(i in 1:length(sub_df_score_matrix_chr22)){
  write.table(sub_df_score_matrix_chr22[[i]],
              paste0("~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/bed_chr22/",
                     colnames(sub_df_score_matrix_chr22[[i]])[4],".bed"),
              quote=F,col.names=F,row.names=F,sep="\t")
}

####
#intersect with bedops-use in terminal
#download bedops from https://github.com/bedops/bedops
#paths conformed for use on local
#bedops --intersect ~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/bed_chr22/* > ~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/intersected_regions_chr22.bed

####
#upload from bedops and store as GRanges
#paths conformed to use on local
intersected_regions_chr22<-read.table("~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/intersected_regions_chr22.bed",header=F,sep="\t")
colnames(intersected_regions_chr22)<-c("chrom","start","end")
GRanges_intersected_regions_chr22<-makeGRangesFromDataFrame(intersected_regions_chr22)
#reading in bed files if needed
dir="~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/bed_chr22"
files<-list.files(path = dir,pattern=".bed$")
sub_df_score_matrix_chr22<-NULL
for(i in 1:length(files)){
  f<-paste0(dir,"/",files[i])
  sub_df_score_matrix_chr22[[i]]<-read.table(f,header=F)
}
for(i in 1:length(sub_df_score_matrix_chr22)){
  colnames(sub_df_score_matrix_chr22[[i]])<-c("chr","start","end","signal")
}

####
#make original dfs into GRanges List
ListOfGRanges_chr22<-lapply(sub_df_score_matrix_chr22,function(x){
  makeGRangesFromDataFrame(x,keep.extra.columns = T)})

####
#getting conserved columns across all GRanges
cols_chr22<-NULL
for(i in 1:length(ListOfGRanges_chr22)){
  target<-ListOfGRanges_chr22[[i]]
  conserved<-GRanges_intersected_regions_chr22
  meta<-mcols(target[which(conserved %over% target),])
  cols_chr22[[i]]<-as.data.frame(meta)
}

####
#combining list of columns into final matrix
mat_chr22<-do.call(cbind,cols_chr22)

####
#noting ChIP-Seq TFs in which the data represent and adding as mat column names.
#Had to look up on encodeproject.org
#should use later on. 
colnames(mat_chr22)
encode<-c("Kap1","Pol2ra","TCFL2","Control.1","Control.2","ZNF263",
          "H3K4me3.1","H3K4me3.2","Control.3","CTCF.1","ELK4","CTCF.2",
          "Control.4")

####
#changing column names, 
#merging replicates and removing controls
newmat_chr22<-mat_chr22
colnames(newmat_chr22)<-encode
newmat_chr22<-newmat_chr22[,-c(4,5,9,13)]
H3K4me3_chr22<-apply(newmat_chr22[,c(5,6)],1,mean)
CTCF_chr22<-apply(newmat_chr22[,c(7,9)],1,mean)
newmat_chr22<-cbind(newmat_chr22[,c(1:4,8)],CTCF_chr22,H3K4me3_chr22)

####
#output final matrix
write.table(newmat_chr22,file="~/GitHub/HMMicro/data/final_matrix_chr22.txt",sep="\t",row.names=F,col.names=T)

####
#saving environment
save.image("~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/DeepBlueR_chr22.RDa")
#####windows#####
####
#making genomic windows
tiling_regions = deepblue_tiling_regions(size=1000, genome="hg19")

####
#requesting signals-in-windows job in DeepBlueR server
request_id<-NULL
for(i in 1:length(exp_columns)){
  #submitting requesting signals in windows job
  request_id[[i]] <- deepblue_score_matrix(
    experiments_columns = exp_columns[i],
    aggregation_function = "mean",
    aggregation_regions_id = tiling_regions)
}

####
#checking example job
deepblue_info(request_id[[4]])$state

####
#dowmloading score matrix once jobs are done
#as indicated above
score_matrix<-NULL
for(i in 1:length(request_id)){
  if(deepblue_info(request_id[[i]])$state=="done"){
    #downloading matrix 
    score_matrix[[i]] <- deepblue_download_request_data(request_id = request_id[[i]])
  }
}

####
#removing NA regions
sub_score_matrix<-NULL
for(i in 1:length(score_matrix)){
  sm<-score_matrix[[i]]
  valinds<-which(!is.na(sm[,4]))
  sub_score_matrix[[i]]<-sm[valinds,]
}

####
#converting data table to data frame class
sub_df_score_matrix<-lapply(sub_score_matrix,as.data.frame)

####
#output beds to directory
for(i in 1:length(sub_df_score_matrix)){
  write.table(sub_df_score_matrix[[i]],
              paste0("~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/bed/",
                     colnames(sub_df_score_matrix[[i]])[4],".bed"),
              quote=F,col.names=F,row.names=F,sep="\t")
}

####
#intersect with bedops-use in terminal
#download bedops from https://github.com/bedops/bedops
#paths conformed for use on local
#bedops --intersect ~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/bed/* > ~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/intersected_regions.bed

####
#upload from bedops and store as GRanges
#paths conformed to use on local
intersected_regions<-read.table("~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/intersected_regions.bed",header=F,sep="\t")
colnames(intersected_regions)<-c("chrom","start","end")
GRanges_intersected_regions<-makeGRangesFromDataFrame(intersected_regions)

####
#make original dfs into GRanges List
ListOfGRanges<-lapply(sub_df_score_matrix,function(x){
  makeGRangesFromDataFrame(x,keep.extra.columns = T)})

####
#getting conserved columns across all GRanges
cols<-NULL
for(i in 1:length(ListOfGRanges)){
  target<-ListOfGRanges[[i]]
  conserved<-GRanges_intersected_regions
  meta<-mcols(target[which(conserved %over% target),])
  cols[[i]]<-as.data.frame(meta)
}

####
#combining list of columns into final matrix
mat<-do.call(cbind,cols)

####
#noting ChIP-Seq TFs in which the data represent and adding as mat column names.
#Had to look up on encodeproject.org
#should use later on. 
colnames(mat)
encode<-c("Kap1","Pol2ra","TCFL2","Control.1","Control.2","ZNF263",
          "H3K4me3.1","H3K4me3.2","Control.3","CTCF.1","ELK4","CTCF.2",
          "Control.4")

####
#changing column names, 
#merging replicates and removing controls
newmat<-mat
colnames(newmat)<-encode
newmat<-newmat[,-c(4,5,9,13)]
H3K4me3<-apply(newmat[,c(5,6)],1,mean)
CTCF<-apply(newmat[,c(7,9)],1,mean)
newmat<-cbind(newmat[,c(1:4,8)],CTCF,H3K4me3)

####
#output final matrix
write.table(newmat,file="~/GitHub/HMMicro/data/final_matrix.txt",sep="\t",row.names=F,col.names=T)

####
#saving environment
save.image("~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/DeepBlueR.RDa")

#####saving requests env#####
save.image("~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/DeepBlueR_requests.RDa")