####################
#Getting HEK293 epigenomic data from
#ENCODE via the DeepBlueR R package.
#
#DeepBlueR allows downloading of curated
#datasets from 4 epigenomics consortia
#
#https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkw211
#
#submitting job requests for 
#epigenetic signal at nucleotide 
#resolution along chr22 in HEK293
#cells in ENCODE
#
#also downloading jobs from DeepBlueServer
#takes a while and chews up memory/cpu...
#
#run interactively e.g. RStudio
#
#used local for downloading (this script) and c2b2 
#(.../HMMicro/Scripts/R/Rjob_chr22.R) for processing
####################

#####Install/load libraries####
source("https://bioconductor.org/biocLite.R")
biocLite("DeepBlueR")
library(DeepBlueR)
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
for(i in 1:length(request_id_chr22)){
  if(deepblue_info(request_id_chr22[[i]])$state=="done"){
    #downloading matrix 
    score_matrix_chr22[[i]] <- deepblue_download_request_data(request_id = request_id_chr22[[i]])
  }
}

#####saving requests env#####
#outputted to a tmp location and kept there because too big for github repo
save.image("~/Documents/Columbia/Courses/COMPUTATIONAL_GENOMICS/HMMicro_notgit/data/DeepBlueR_chr22_score_matrix.RDa")
