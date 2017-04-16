# tRNAs.sql was originally generated by the autoSql program, which also 
# generated tRNAs.c and tRNAs.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#transfer RNA genes
CREATE TABLE tRNAs (
    chrom varchar(255) not null,	# chromosome
    chromStart int unsigned not null,	# Start position in chromosome
    chromEnd int unsigned not null,	# End position in chromosome
    name varchar(255) not null,	# transfer RNA gene name
    score int unsigned not null,	# Score from 900-1000.  1000 is best
    strand char(1) not null,	# Value should be + or -
    aa varchar(255) not null,	# Amino acid for the tRNA
    ac varchar(255) not null,	# Anticodon for the tRNA
    intron varchar(255) not null,	# Coordinates for intron
    trnaScore float not null,	# tRNAScanSE score
    genomeUrl varchar(255) not null,	# GtRNAdb genome summary URL
    trnaUrl varchar(255) not null,	# GtRNAdb tRNA alignment URL
              #Indices
    PRIMARY KEY(chrom, chromStart, chromEnd)
);
