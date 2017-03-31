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
#output final matrix
write.table(mat,file="~/GitHub/HMMicro/data/final_matrix.txt",sep="\t",row.names=F,col.names=T)
