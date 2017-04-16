/* rmskOut.h was originally generated by the autoSql program, which also 
 * generated rmskOut.c and rmskOut.sql.  This header links the database and the RAM 
 * representation of objects. */

/* Copyright (C) 2003 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#ifndef RMSKOUT_H
#define RMSKOUT_H

#ifndef LINEFILE_H
#include "linefile.h"
#endif

struct rmskOut
/* RepeatMasker .out record */
    {
    struct rmskOut *next;  /* Next in singly linked list. */
    unsigned swScore;	/* Smith Waterman alignment score */
    unsigned milliDiv;	/* Base mismatches in parts per thousand */
    unsigned milliDel;	/* Bases deleted in parts per thousand */
    unsigned milliIns;	/* Bases inserted in parts per thousand */
    char *genoName;	/* Genomic sequence name */
    unsigned genoStart;	/* Start in genomic sequence */
    unsigned genoEnd;	/* End in genomic sequence */
    int genoLeft;	/* Size left in genomic sequence */
    char strand[2];	/* Relative orientation + or - */
    char *repName;	/* Name of repeat */
    char *repClass;	/* Class of repeat */
    char *repFamily;	/* Family of repeat */
    int repStart;	/* Start in repeat sequence */
    int repEnd;		/* End in repeat sequence */
    int repLeft;	/* Size left in repeat sequence */
    char id[2];		/* '*' or ' '.  I don't know what this means */
    };

void rmskOutStaticLoad(char **row, struct rmskOut *ret);
/* Load a row from rmskOut table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct rmskOut *rmskOutLoad(char **row);
/* Load a rmskOut from row fetched with select * from rmskOut
 * from database.  Dispose of this with rmskOutFree(). */

struct rmskOut *rmskOutCommaIn(char **pS, struct rmskOut *ret);
/* Create a rmskOut out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new rmskOut */

void rmskOutFree(struct rmskOut **pEl);
/* Free a single dynamically allocated rmskOut such as created
 * with rmskOutLoad(). */

void rmskOutFreeList(struct rmskOut **pList);
/* Free a list of dynamically allocated rmskOut's */

void rmskOutOutput(struct rmskOut *el, FILE *f, char sep, char lastSep);
/* Print out rmskOut.  Separate fields with sep. Follow last field with lastSep. */

#define rmskOutTabOut(el,f) rmskOutOutput(el,f,'\t','\n');
/* Print out rmskOut as a line in a tab-separated file. */

#define rmskOutCommaOut(el,f) rmskOutOutput(el,f,',',',');
/* Print out rmskOut as a comma separated list including final comma. */

/* ------------ End of AutoSQL generated code. ------------------ */

struct rmskOut *rmskOutRead(char *fileName);
/* Read all records in .out file and return as list. */

void rmskOutOpenVerify(char *fileName, struct lineFile **retFile, boolean *retEmpty);
/* Open repeat masker .out file and verify that it is good.
 * Set retEmpty if it has header characteristic of an empty file. */

struct rmskOut *rmskOutReadNext(struct lineFile *lf);
/* Read next record from repeat masker file.  Return NULL at EOF. */

void rmskOutWriteHead(FILE *f);
/* Write out rmsk header lines. */

void rmskOutWriteOneOut(struct rmskOut *rmsk, FILE *f);
/* Write one rmsk in .out format to file. */

void rmskOutWriteAllOut(char *fileName, struct rmskOut *rmskList);
/* Write .out format file containing all in rmskList. */

struct binKeeper *readRepeats(char *chrom, char *rmskFileName, struct hash *tSizeHash);
/* read all repeats for a chromosome of size size, returns results in binKeeper structure for fast query*/

struct hash *readRepeatsAll(char *sizeFileName, char *rmskDir);
/* read all repeats for a all chromosomes getting sizes from sizeFileNmae , returns results in hash of binKeeper structure for fast query*/
#endif /* RMSKOUT_H */


