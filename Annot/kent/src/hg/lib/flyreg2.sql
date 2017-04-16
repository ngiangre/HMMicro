# flyreg2.sql was originally generated by the autoSql program, which also 
# generated flyreg2.c and flyreg2.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#Flyreg data (version 2) from Casey Bergman
CREATE TABLE flyreg2 (
    bin smallint unsigned not null,   # bin for speed
    chrom varchar(255) not null,	# Chromosome
    chromStart int unsigned not null,	# Start position in chromosome
    chromEnd int unsigned not null,	# End position in chromosome
    name varchar(255) not null,	# Factor
    target varchar(255) not null,	# Target
    pmid int unsigned not null,	# PubMed ID
    fpid int unsigned not null,	# Footprint ID -- stable ID across versions
              #Indices
    INDEX (chrom(12), bin),
    INDEX (chromStart, bin)
);
