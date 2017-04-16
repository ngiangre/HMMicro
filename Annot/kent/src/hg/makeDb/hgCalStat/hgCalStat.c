/* hgCalStat - calculate genome assembly statistics */

/* Copyright (C) 2013 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */
#include "common.h"
#include "hCommon.h"
#include "hdb.h"

void usage()
/* Explain usage and exit. */
{
errAbort(
  "hgCalStat - calculate genome assembly statistics\n"
  "usage:\n"
  "   hgCalStat statFile db outFile\n"
  "      statFile the output file generated by stats.pl\n"
  "      db is the genome db\n"
  "      outFile is the output file\n"
  "example:\n"
  "   hgCalStat stats.pl.out hg18 hg18.out\n");
}

int main(int argc, char *argv[])
{
struct sqlConnection *conn2, *conn3;
struct sqlResult *sr2;
char **row2;
char query2[256];
char cond_str[500];

char *statFileName;
boolean chromDone = FALSE;
boolean chromRandomDone = FALSE;
boolean hasRandom;
char *chromRandom;
int chromRandomSize;
int orderedCnt = 0;
int randomCnt  = 0;
int n50Total   = 0;
int n50RandomTotal = 0;

int  halfTotalSize;
int  halfTotalRandomSize;

char *outfileName;
FILE *inf, *outf;

char line[1000];
char line2[1000];
char overallLine[1000];

char *chrom;
int chromSize;
char *ctgSize;

char *chp, *chp1;

char *genomeDBname;
char *contig;
unsigned long sumL;

int sum, sumRandom;
char *ctgCnt, *ctgRandomCnt;

if (argc != 4) usage();

statFileName = argv[1];
genomeDBname = argv[2];
outfileName  = argv[3];

hSetDb(genomeDBname);

outf = mustOpen(outfileName, "w");
    
conn2= hAllocConn();
conn3= hAllocConn();
    
inf   = mustOpen(statFileName, "r");

/* skip the labels line */
fgets(line, 1000, inf);

while (fgets(line, 1000, inf) != NULL)
    {
    *(line+strlen(line)-1) = '\0';
    strcpy(line2, line);

    chp = strstr(line2, "\t");
    *chp = '\0';

    /* special processing needed for overall line */
    if (sameWord(line2, "overallChrom")) 
    	{
	strcpy(overallLine, line);
	break;
	}
    chp1 = chp + 1;
    
    chrom = strdup(line2);
    sprintf(line, "%s_random", chrom);
    chromRandom = strdup(line);
    
    sqlSafefFrag(cond_str, sizeof cond_str, "chrom='%s'", chrom);

    /* reformat and output the first 9 columns */
    chp = strstr(chrom, "chr")+strlen("chr");
    fprintf(outf,"%5s     ", chp);
    chp = strstr(chp1, "\t");
    while (chp != NULL)
    	{
	*chp = '\0';
	fprintf(outf,"%11s", chp1);
	chp1 = chp+1;
        chp = strstr(chp1, "\t");
	}
    fprintf(outf,"%19s", chp1);

    chromSize = atoi(sqlGetField(genomeDBname, "chromInfo", "size", cond_str));

    sqlSafefFrag(cond_str, sizeof cond_str, "chrom='%s'", chromRandom);
    if (sqlGetField(genomeDBname, "chromInfo", "size", cond_str) != NULL)
    	{
	hasRandom = TRUE;
        chromRandomSize = atoi(sqlGetField(genomeDBname, "chromInfo", "size", cond_str));
	}
    else
    	{
	hasRandom = FALSE;
        chromRandomSize = 0;
	}
    
    sqlSafefFrag(cond_str, sizeof cond_str, "chrom='%s'", chrom);
    ctgCnt = sqlGetField(genomeDBname, "ctgPos", "count(*)", cond_str);

    fprintf(outf,"%12s", ctgCnt);
    orderedCnt = orderedCnt + atoi(ctgCnt);

    /* calculate N50 for ordered contigs */
    sqlSafef(query2, sizeof query2, "select contig, size from %s.ctgPos where chrom='%s' order by size desc", genomeDBname, chrom);
    	
    sr2 = sqlMustGetResult(conn2, query2);
    row2 = sqlNextRow(sr2);
    sum = 0;
    chromDone = FALSE;
    while ((row2 != NULL) && !chromDone)
	{
 	contig  = row2[0];
	ctgSize = row2[1];
	
	sum = sum + atoi(ctgSize);

	if (sum >= (chromSize/2))
	    {
	    fprintf(outf,"%10s", ctgSize);
	    chromDone = TRUE;
	    }

	row2 = sqlNextRow(sr2);
	}
	
    if (!chromDone) fprintf(outf,"%10s", "?");
    sqlFreeResult(&sr2);
    
    /* calculate N50 for random contigs */
    if (hasRandom)
	{
        sqlSafefFrag(cond_str, sizeof cond_str, "chrom='%s'", chromRandom);
        ctgRandomCnt = sqlGetField(genomeDBname, "ctgPos", "count(*)", cond_str);
	
    	sqlSafef(query2, sizeof query2,
		"select contig, size from %s.ctgPos where chrom='%s' order by size desc", 
		genomeDBname, chromRandom);
    	
    	sr2 = sqlMustGetResult(conn2, query2);
    	row2 = sqlNextRow(sr2);
    	sumRandom = 0;
    	chromRandomDone = FALSE;
    	while ((row2 != NULL) && !chromRandomDone)
	    {
 	    contig  = row2[0];
	    ctgSize = row2[1];
	
	    sumRandom = sumRandom + atoi(ctgSize);

	    if (sumRandom >= (chromRandomSize/2))
	    	{
	    	fprintf(outf,"%12s%9s", ctgRandomCnt, ctgSize);
	        randomCnt = randomCnt + atoi(ctgRandomCnt);
	    	chromRandomDone = TRUE;
		}

	    row2 = sqlNextRow(sr2);
	    }
        sqlFreeResult(&sr2);
	}
    else
	{
	fprintf(outf,"%12s%9s", "-", "-");
	}
    fprintf(outf,"\n"); 
    }

/* first calculate total chromosome size */
sqlSafef(query2, sizeof query2,
    "select chrom, size from %s.ctgPos where chrom not like '%crandom%c' and chrom not like '%chap%c' order by size desc", 
    genomeDBname, '%', '%', '%', '%');
    	
sr2 = sqlMustGetResult(conn2, query2);
row2 = sqlNextRow(sr2);
sumL = 0;
while ((row2 != NULL) && !chromDone)
    {
    chrom   = row2[0];
    ctgSize = row2[1];
	
    sumL = sumL + atoi(ctgSize);
	
    row2 = sqlNextRow(sr2);
    }
sqlFreeResult(&sr2);
halfTotalSize = (int)(sumL/2L);

/* calculate N50 totoal */
sqlSafef(query2, sizeof query2,
    "select chrom, size from %s.ctgPos where chrom not like '%crandom%c' and chrom not like '%chap%c' order by size desc", 
    genomeDBname, '%', '%', '%', '%');
    	
sr2 = sqlMustGetResult(conn2, query2);
row2 = sqlNextRow(sr2);
sumL = 0L;
while ((row2 != NULL) && !chromDone)
    {
    chrom   = row2[0];
    ctgSize = row2[1];
	
    sumL = sumL + atoi(ctgSize);
	
    if ((int)sumL >= halfTotalSize)
    	{
	n50Total = atoi(ctgSize);
	break;
	}

    row2 = sqlNextRow(sr2);
    }
	
sqlFreeResult(&sr2);

/* calculate total random chromosome size */
sqlSafef(query2, sizeof query2,
    "select chrom, size from %s.ctgPos where chrom like '%crandom%c' order by size desc", 
    genomeDBname, '%', '%');
    	
sr2 = sqlMustGetResult(conn2, query2);
row2 = sqlNextRow(sr2);
sumL = 0;
while ((row2 != NULL) && !chromDone)
    {
    chrom   = row2[0];
    ctgSize = row2[1];
	
    sumL = sumL + atoi(ctgSize);
	
    row2 = sqlNextRow(sr2);
    }
sqlFreeResult(&sr2);
halfTotalRandomSize = (int)(sumL/2L);

/* calculate N50 of randome total */
sqlSafef(query2, sizeof query2,
    "select chrom, size from %s.ctgPos where chrom like '%crandom%c' order by size desc", 
    genomeDBname, '%', '%');
    	
sr2 = sqlMustGetResult(conn2, query2);
row2 = sqlNextRow(sr2);
sumL = 0L;
while ((row2 != NULL) && !chromDone)
    {
    chrom   = row2[0];
    ctgSize = row2[1];
	
    sumL = sumL + atoi(ctgSize);

    if ((int)sumL >= halfTotalRandomSize)
    	{
	n50RandomTotal = atoi(ctgSize);
	break;
	}

    row2 = sqlNextRow(sr2);
    }
	
sqlFreeResult(&sr2);

fprintf(outf,"%s\n", " ----------------------------------------------------------------------------------------------------------------------------------------------------");

/* process overall line */
fprintf(outf," Overall\n");
fprintf(outf,"   Chrom  ");

chp1 = strstr(overallLine, "\t");
chp1++;
chp = strstr(chp1, "\t");
while (chp != NULL)
    {
    *chp = '\0';
    fprintf(outf,"%11s", chp1);
    chp1 = chp+1;
    chp = strstr(chp1, "\t");
    }
fprintf(outf,"%19s", chp1);

/* print out overall N50s */
fprintf(outf,"%12d%10d%12d%9d\n", orderedCnt, n50Total, randomCnt, n50RandomTotal);
fprintf(outf,"\n");	

hFreeConn(&conn2);
	
fclose(outf);
return(0);
}

