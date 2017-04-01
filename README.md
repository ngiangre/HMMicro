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

### Evaluation data

Training data was taken from (Hafner, et al.)[http://dx.doi.org/10.1016/j.cell.2010.03.009]. Specifically, these are crosslink-centered regions (CCRs) from AGO-PAR-CLIP sequencing libraries, i.e. purportedly miRNA sites, from HEK293 cells. The files are located at *CLIP_data/HEK293*. Steps/files are as follows:

1.  Hafner.combined_AGO.raw.txt, Hafner.combined_TNRC6.raw.txt
    -   Table S4 from the Hafner paper. Fasta-list of the CCRs from the combined AGO-PAR-CLIP libraries and the combined TNRC6-PAR-CLIP libraries used for the analysis
2.  process.raw.py
    -   Python script for converting raw fasta-list of CCRs to a more readable tsv format with header
    -   Outputs:
        *   Hafner.combined_AGO.proc.tsv
3.  create_beds.sh
    -   Create bed files to be lifted over from hg18 to hg19
    -   Outputs:
        *   Hafner.combined_AGO.proc.bed
        *   Hafner.combined_AGO.proc.start.end.bed
        *   Hafner.combined_AGO.proc.start.t2c_pos.bed
4.  Liftover was performed using the UCSC Genome Browser Utilities (liftOver tool)[https://genome.ucsc.edu/cgi-bin/hgLiftOver].
    -   Outputs:
        *   Hafner.combined_AGO.proc.start.end.hg19.bed
        *   Hafner.combined_AGO.proc.start.t2c_pos.hg19.bed
5.  merge_hg19_beds.sh
    -   Merge the lifted-over beds into a single bed file
    -   Outputs:
        Hafner.combined_AGO.proc.hg19.bed
### HMM training

### HMM testing

### HMM evaluation

**************************************************
## Data and Scripts

*intersected_regions.bed* contains the genomic windows considered. They are the rows in the *final_matrix.txt*. The columns in *final_matrix.txt* are the epigenetic marks (signal from BigWig files) in HEK293 cells from ENCODE. These are Kap1 ChIP-Seq, Pol2ra ChIPSeq, TCFL2 ChIP-Seq, ZNF263 ChIP-Seq, ELK4 ChIP-Seq, CTCF ChIP-Seq and H3K4me3 ChIP-Seq. The procedure for generating both files can be found in the *DeepBlueR.R* script. *dimReduce.ipynb* reduces the matrix in *final_matrix.txt* to its principal components, outputting a square matrix into *matrix_reduce.txt*. Each row is a coefficient in the linear combination for a given principal component, and each column is a epigenetic mark's (named in order above) principal component, where n in the n x n matrix corresponds to the number of epigenetic mark signatures. Along with our training dataset, the components will be used for training a HMM. 

## References

1. Elbrecht et al. (2017). DeepBlueR: large-scale epigenomic analysis in R. *Bioinformatics* 1-2. 
