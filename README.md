# HMMicro
Spring 2017 Computational Genomics course project
**************************************************
## Overview

We developed a Hidden Markov Model (HMM) for predicting miRNA binding sites. Underlying the HMM is the hypothess that the arrangement of different epigenetic sequencing data types gives (redundant) regulatory information in a given cell type, suitable for predicting other epigenetic patterns such as miRNA binding.

Here, we have experimentally validated miRNA binding sites [1] and we hypothesize that the signature of epigenetic marks, provided from ENCODE HEK293 cells [2], along the genome contain a hidden 'regulatory logic' that is present in guiding miRNA binding along the transcriptome. 

We provide an HMM that predicts binding and non-binding states, trained on the different epigenomic data types from HEK293 cells assayed by ENCODE. 

For reproducing without consulting the explanations given below, pretty much follow the ordering of script execution: *DeepBlueR.R* --> *Rjob_chr22.R* --> *dimReduce.ipynb*
**************************************************
### HEK293 ENCODE data downloading and processing of epigenetic matrix

We downloaded epigenetic data derived from HEK293 cells assayed by ENCODE via the DeepBlue server accessed through the DeepBlueR R package [3]. We retrieved signal from bigWig files from seven different ChIP-Seq experiments, with nucleotide resolution along chromosome 22, yeilding signatures for each epigemetic mark. 

See *DeepBlueR.R* and *Rjob_chr22.R* script for detailed procedures.

### Dimension reduction of epigenetic matrix

The dimensions of the epigenetic mark signatures in HEK293 cells were reduced to its principal components to provide a prior emissions matrix for training the HMM.

See *dimReduce.ipynb* script for detailed procedure.

### Evaluation data

Our evaluation set was the experimentally validated miRNA binding sites [2], which specifically included crosslink-centered regions (CCRs) from AGO-PAR-CLIP sequencing libraries, from HEK293 cells. The files are located at *CLIP_data/HEK293*. Steps/files are as follows:

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
4.  Liftover was performed using the UCSC Genome Browser Utilities [liftOver tool](https://genome.ucsc.edu/cgi-bin/hgLiftOver).
    -   Outputs:
        *   Hafner.combined_AGO.proc.start.end.hg19.bed
        *   Hafner.combined_AGO.proc.start.t2c_pos.hg19.bed
5.  merge_hg19_beds.sh
    -   Merge the lifted-over beds into a single bed file
    -   Outputs:
        Hafner.combined_AGO.proc.hg19.bed
### HMM training

We developed and trained an HMM using the *Pomegranate* Python package [4]. The principal components in *reduced_matrix.txt*  are used as the prior emission matrix. Uniformly disyributed random numbers populated a prior transition matrix to model probably transitions between binding and non-binding states. These prior emissions and transmissions probabilities were used for hard coding and then training the HMM. The detailed procedures are laid out in *train_hmm.ipynb*.

### HMM testing and evaluation

### Results and Future Directions

**************************************************
## Data and Scripts

*final_matrix.txt* (not in the repo because the size is too large) contains normalized epigenetic signal at nucleotide resolution across chromosome 22 for the seven different epigenetic marks from ENCODE HEK293 cells. The columns in *final_matrix.txt* are the epigenetic marks (signal from BigWig files) in HEK293 cells from ENCODE. These are Kap1 ChIP-Seq, Pol2ra ChIPSeq, TCFL2 ChIP-Seq, ZNF263 ChIP-Seq, ELK4 ChIP-Seq, CTCF ChIP-Seq and H3K4me3 ChIP-Seq. 

The *DeepBlueR.R* and *Rjob_chr22.R* scripts detail the procedure for downloading the data and generating the matrix file from *final_matrix.txt*. 

*dimReduce.ipynb* reduces the matrix in *final_matrix.txt* to its principal components, outputting a square matrix into *matrix_reduce.txt*. Each row is a coefficient in the linear combination for a given principal component, and each column is a principal component from *final_matrix.txt*, where n in the n x n matrix corresponds to the number of epigenetic mark signatures, herein seven. 

*findingStatesDistribution.ipynb* details the procedure for finding an approximate distribution of the evaluation dataset's miRNA start sites along chromosome 22. Unfortunately, we didn't have time to strategize how to implement this knowledge. This would be a future direction. 

*train_hmm_chr22.py* details the procedure for developing and training an HMM using the *Pomegranate* Python package [4]. The principal components in *reduced_matrix.txt*  are used as the prior emission matrix. Uniformly disyributed random numbers populated a prior transition matrix to model probably transitions between binding and non-binding states. These prior emissions and transmissions probabilities were used for hard coding and then training the HMM. 
**************************************************
## System Requirements

**R: version 3.3.1; *DeepBlueR* package.**

**Python: 2.7+; *sklearn*, *pandas*, *numpy*, *matplotlib* and *pomegranate* packages**

**HPC/Laptop requirements: system that can handle >80GB of memory and processing-for processing of epigenetic data detailed in *Rjob_chr22.R* script. Local laptop with 8GB of RAM and >20GB of memory on hard disk will suffice for submitting and downloading from DeepBlue server, though one generally can't do anything else...so run overnight with laptop never sleeping.**
**************************************************
## References

1. Lipchina, I. et al. (2011). Genome-wide identification of microRNA targets in human ES cells reveals a role for miR-302 in modulating BMP response. *Genes and Development.* 25: 2173–2186. 

2. The ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. *Nature*. 489(7414): 57-74. 

3. Elbrecht et al. (2017). DeepBlueR: large-scale epigenomic analysis in R. *Bioinformatics* 1-2. 

4. Jacob Schreiber (2016). pomegranate. *GitHub*. https://github.com/jmschrei/pomegranate. 
