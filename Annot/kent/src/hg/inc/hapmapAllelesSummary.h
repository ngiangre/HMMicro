/* hapmapAllelesSummary.h was originally generated by the autoSql program, which also 
 * generated hapmapAllelesSummary.c and hapmapAllelesSummary.sql.  This header links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2007 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#ifndef HAPMAPALLELESSUMMARY_H
#define HAPMAPALLELESSUMMARY_H

#define HAPMAPALLELESSUMMARY_NUM_COLS 27

struct hapmapAllelesSummary
/* HapMap allele summaries for filtering */
    {
    struct hapmapAllelesSummary *next;  /* Next in singly linked list. */
    char *chrom;	/* Chromosome */
    unsigned chromStart;	/* Start position in chrom (0 based) */
    unsigned chromEnd;	/* End position in chrom (1 based) */
    char *name;	/* Reference SNP identifier from dbSnp */
    unsigned score;	/* Use for heterozygosity */
    char strand[2];	/* Which genomic strand contains the observed alleles */
    char *observed;	/* Observed string from genotype file */
    char allele1[2];	/* This allele has been observed */
    char *allele2;	/* This allele may not have been observed */
    unsigned popCount;	/* How many populations have data */
    char *isMixed;	/* Are there different major alleles? */
    char *majorAlleleCEU;	/* major allele for the CEU population */
    unsigned majorAlleleCountCEU;	/* major allele count for the CEU population */
    unsigned totalAlleleCountCEU;	/* total allele count for the CEU population */
    char *majorAlleleCHB;	/* major allele for the CHB population */
    unsigned majorAlleleCountCHB;	/* major allele count for the CHB population */
    unsigned totalAlleleCountCHB;	/* total allele count for the CHB population */
    char *majorAlleleJPT;	/* major allele for the JPT population */
    unsigned majorAlleleCountJPT;	/* major allele count for the JPT population */
    unsigned totalAlleleCountJPT;	/* total allele count for the JPT population */
    char *majorAlleleYRI;	/* major allele for the YRI population */
    unsigned majorAlleleCountYRI;	/* major allele count for the YRI population */
    unsigned totalAlleleCountYRI;	/* total allele count for the YRI population */
    char *chimpAllele;	/* chimp allele */
    unsigned chimpAlleleQuality;	/* quality score (0-100) */
    char *macaqueAllele;	/* macaque allele */
    unsigned macaqueAlleleQuality;	/* quality score (0-100) */
    };

void hapmapAllelesSummaryStaticLoad(char **row, struct hapmapAllelesSummary *ret);
/* Load a row from hapmapAllelesSummary table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct hapmapAllelesSummary *hapmapAllelesSummaryLoad(char **row);
/* Load a hapmapAllelesSummary from row fetched with select * from hapmapAllelesSummary
 * from database.  Dispose of this with hapmapAllelesSummaryFree(). */

struct hapmapAllelesSummary *hapmapAllelesSummaryLoadAll(char *fileName);
/* Load all hapmapAllelesSummary from whitespace-separated file.
 * Dispose of this with hapmapAllelesSummaryFreeList(). */

struct hapmapAllelesSummary *hapmapAllelesSummaryLoadAllByChar(char *fileName, char chopper);
/* Load all hapmapAllelesSummary from chopper separated file.
 * Dispose of this with hapmapAllelesSummaryFreeList(). */

#define hapmapAllelesSummaryLoadAllByTab(a) hapmapAllelesSummaryLoadAllByChar(a, '\t');
/* Load all hapmapAllelesSummary from tab separated file.
 * Dispose of this with hapmapAllelesSummaryFreeList(). */

struct hapmapAllelesSummary *hapmapAllelesSummaryCommaIn(char **pS, struct hapmapAllelesSummary *ret);
/* Create a hapmapAllelesSummary out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new hapmapAllelesSummary */

void hapmapAllelesSummaryFree(struct hapmapAllelesSummary **pEl);
/* Free a single dynamically allocated hapmapAllelesSummary such as created
 * with hapmapAllelesSummaryLoad(). */

void hapmapAllelesSummaryFreeList(struct hapmapAllelesSummary **pList);
/* Free a list of dynamically allocated hapmapAllelesSummary's */

void hapmapAllelesSummaryOutput(struct hapmapAllelesSummary *el, FILE *f, char sep, char lastSep);
/* Print out hapmapAllelesSummary.  Separate fields with sep. Follow last field with lastSep. */

#define hapmapAllelesSummaryTabOut(el,f) hapmapAllelesSummaryOutput(el,f,'\t','\n');
/* Print out hapmapAllelesSummary as a line in a tab-separated file. */

#define hapmapAllelesSummaryCommaOut(el,f) hapmapAllelesSummaryOutput(el,f,',',',');
/* Print out hapmapAllelesSummary as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* HAPMAPALLELESSUMMARY_H */

