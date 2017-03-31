# HMMicro
Spring 2017 Computational Genomics course project

## Overview

We are developing a Hidden Markov Model (HMM) for predicting miRNA binding sites. We hypothesize that the epigenetic arrangement of different epigenetic sequencing data types gives (redundant) regulatory information for miRNA binding in a given cell type.

Here, we have experimentally validated miRNA binding sites (list reference) and we hypothesize that the signature of epigenetic marks along the genome contain a hidden 'regulatory logic' that is present in guiding miRNA binding along the transcriptome. 

We will provide a HMM that will predict binding and non-binding states, trained on different epigenomic data types from HEK293 cells assayed by ENCODE. 

**************************************************
### HEK293 ENCODE data downloading and processing of epigenetic matrix

### Dimesion reduction of epigenetic matrix

### HMM training

### HMM testing

### HMM evaluation
**************************************************
## Data and Scripts

*intersected_regions.bed* contains the genomic windows considered. They are the rows in the *final_matrix.txt*. The columns in *final_matrix.txt* are the epigenetic marks in HEK293 cells from ENCODE. Both files were generated the the *DeepBlueR.R* script. *dimReduce.ipynb* reduces the matrix in *final_matrix.txt* to it's principal components, outputting a square matrix into *matrix_reduce.txt* and, along with our training dataset, will be used for training a HMM. 
