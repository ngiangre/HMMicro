# affyOffset.sql was originally generated by the autoSql program, which also 
# generated affyOffset.c and affyOffset.h.  This creates the database representation of
# an object which can be loaded and saved from RAM in a fairly 
# automatic way.

#File format giving offset of Affymetrix probe sets into contigs.
CREATE TABLE affyOffset (
    piece varchar(255) not null,	# Name of 'piece' of genome, something like ctg21fin1piece100
    tStart int unsigned not null,	# Start of 'piece' in contig.
              #Indices
    PRIMARY KEY(piece)
);
