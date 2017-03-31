# HMMicro
Spring 2017 Computational Genomics course project

## Overview

We are developing a Hidden Markov Model (HMM) for predicting miRNA binding sites. We hypothesize that the epigenetic arrangement of different epigenetic sequencing data types gives (redundant) regulatory information for miRNA binding in a given cell type.

Here, we have experimentally validated miRNA binding sites (list reference) and we hypothesize that the signature of epigenetic marks along the genome contain a hidden 'regulatory logic' that is present in guiding miRNA binding along the transcriptome. 

We will provide a HMM that will predict binding and non-binding states, trained on different epigenomic data types from HEK293 cells assayed by ENCODE. 

**************************************************
### HEK293 ENCODE data downloading and processing of epigenetic matrix

We downloaded epigenetic data derived from HEK293 cells assayed by ENCODE via the DeepBlue server accessed through the DeepBlueR package [1]. We retrieved signal from bigWig files for various ChIP-Seq experiments, surveying 1000bp genomic windows genome-wideyeilding signatures for each epigemetic mark. 

See *DeepBlueR.R* script for detailed procedure.

### Dimesion reduction of epigenetic matrix

The dimensions of the epigenetic mark signatures in HEK293 cells were reduced to its principal components.

See *dimReduce.ipynb* script for detailed procedure.

### HMM training

### HMM testing

### HMM evaluation
**************************************************
## Data and Scripts

*intersected_regions.bed* contains the genomic windows considered. They are the rows in the *final_matrix.txt*. The columns in *final_matrix.txt* are the epigenetic marks (signal from BigWig files) in HEK293 cells from ENCODE. These are Kap1 ChIP-Seq, Pol2ra ChIPSeq, TCFL2 ChIP-Seq, Control rep 1 ChIP-Seq, Control rep 2 ChIP-Seq, ZNF263 ChIP-Seq, H3K4me3 rep 1 ChIP-Seq, H3K4me3 rep 2 ChIP-Seq, Control rep 3 ChIP-Seq, CTCF rep 1 ChIP-Seq, ELK4 ChIP-Seq, CTCF rep 2 ChIP-Seq, and Control rep 4 ChIP-Seq. The procedure for generating both files can be found in the *DeepBlueR.R* script. *dimReduce.ipynb* reduces the matrix in *final_matrix.txt* to it's principal components, outputting a square matrix into *matrix_reduce.txt*. Each row is a coefficient in the linear combination for a given principal component, and each column is a principal component, where n in the n x n matrix corresponds to the number of epigenetic mark signatures. Along with our training dataset, will be used for training a HMM. 

## References

1. Elbrecht et al. (2017). DeepBlueR: large-scale epigenomic analysis in R. *Bioinformatics* 1-2. 
