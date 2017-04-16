/* pfamDat.h was originally generated by the autoSql program, which also 
 * generated pfamDat.c and pfamDat.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef PFAMDAT_H
#define PFAMDAT_H

#define PFAMHIT_NUM_COLS 5

struct pfamHit
/* Data about a pfam global hit. */
    {
    struct pfamHit *next;  /* Next in singly linked list. */
    char *model;	/* Model hit with sequence. */
    char *descript;	/* Description of Model hit with sequence. */
    double score;	/* Bit score of hit. */
    double eval;	/* E-values of hit. */
    unsigned numTimesHit;	/* Number of times hit was found in sequence. */
    };

void pfamHitStaticLoad(char **row, struct pfamHit *ret);
/* Load a row from pfamHit table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct pfamHit *pfamHitLoad(char **row);
/* Load a pfamHit from row fetched with select * from pfamHit
 * from database.  Dispose of this with pfamHitFree(). */

struct pfamHit *pfamHitLoadAll(char *fileName);
/* Load all pfamHit from whitespace-separated file.
 * Dispose of this with pfamHitFreeList(). */

struct pfamHit *pfamHitLoadAllByChar(char *fileName, char chopper);
/* Load all pfamHit from chopper separated file.
 * Dispose of this with pfamHitFreeList(). */

#define pfamHitLoadAllByTab(a) pfamHitLoadAllByChar(a, '\t');
/* Load all pfamHit from tab separated file.
 * Dispose of this with pfamHitFreeList(). */

struct pfamHit *pfamHitCommaIn(char **pS, struct pfamHit *ret);
/* Create a pfamHit out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new pfamHit */

void pfamHitFree(struct pfamHit **pEl);
/* Free a single dynamically allocated pfamHit such as created
 * with pfamHitLoad(). */

void pfamHitFreeList(struct pfamHit **pList);
/* Free a list of dynamically allocated pfamHit's */

void pfamHitOutput(struct pfamHit *el, FILE *f, char sep, char lastSep);
/* Print out pfamHit.  Separate fields with sep. Follow last field with lastSep. */

#define pfamHitTabOut(el,f) pfamHitOutput(el,f,'\t','\n');
/* Print out pfamHit as a line in a tab-separated file. */

#define pfamHitCommaOut(el,f) pfamHitOutput(el,f,',',',');
/* Print out pfamHit as a comma separated list including final comma. */

#define PFAMDHIT_NUM_COLS 12

struct pfamDHit
/* Data about a pfam domain hit. */
    {
    struct pfamDHit *next;  /* Next in singly linked list. */
    char *model;	/* Model of domain. */
    int domain;	/* Domain number of hit. */
    int numDomain;	/* Number of domains for a hit. */
    unsigned seqStart;	/* Start in sequence of hit. */
    unsigned seqEnd;	/* End in sequence of hit. */
    char *seqRep;	/* String representation of where hit is located in seq, '[.','..','.]','[]' */
    unsigned hmmStart;	/* Start in profile hmm of hit. */
    unsigned hmmEnd;	/* End in profile hmm of hit. */
    char *hmmRep;	/* String representation of where hit is located in seq, '[.','..','.]','[]' */
    double dScore;	/* Score for domain hit. */
    double dEval;	/* Evalue for domain. */
    char *alignment;	/* Text based alignment. */
    };

void pfamDHitStaticLoad(char **row, struct pfamDHit *ret);
/* Load a row from pfamDHit table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct pfamDHit *pfamDHitLoad(char **row);
/* Load a pfamDHit from row fetched with select * from pfamDHit
 * from database.  Dispose of this with pfamDHitFree(). */

struct pfamDHit *pfamDHitLoadAll(char *fileName);
/* Load all pfamDHit from whitespace-separated file.
 * Dispose of this with pfamDHitFreeList(). */

struct pfamDHit *pfamDHitLoadAllByChar(char *fileName, char chopper);
/* Load all pfamDHit from chopper separated file.
 * Dispose of this with pfamDHitFreeList(). */

#define pfamDHitLoadAllByTab(a) pfamDHitLoadAllByChar(a, '\t');
/* Load all pfamDHit from tab separated file.
 * Dispose of this with pfamDHitFreeList(). */

struct pfamDHit *pfamDHitCommaIn(char **pS, struct pfamDHit *ret);
/* Create a pfamDHit out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new pfamDHit */

void pfamDHitFree(struct pfamDHit **pEl);
/* Free a single dynamically allocated pfamDHit such as created
 * with pfamDHitLoad(). */

void pfamDHitFreeList(struct pfamDHit **pList);
/* Free a list of dynamically allocated pfamDHit's */

void pfamDHitOutput(struct pfamDHit *el, FILE *f, char sep, char lastSep);
/* Print out pfamDHit.  Separate fields with sep. Follow last field with lastSep. */

#define pfamDHitTabOut(el,f) pfamDHitOutput(el,f,'\t','\n');
/* Print out pfamDHit as a line in a tab-separated file. */

#define pfamDHitCommaOut(el,f) pfamDHitOutput(el,f,',',',');
/* Print out pfamDHit as a comma separated list including final comma. */

#define PFAMDAT_NUM_COLS 4

struct pfamDat
/* Structure to hold results of one hmmpfam run. Distributed by Sean Eddy. See http://www.genetics.wustl.edu/eddy/software/ */
    {
    struct pfamDat *next;  /* Next in singly linked list. */
	char *header;      /* Usually NULL, just used for testing in diffs. */
    char *seqName;	/* Sequence name. */
    char *sequence;	/* Sequence run against library. */
    struct pfamHit *pfamHitList;	/* Global hits. */
    struct pfamDHit *pfamDHitList;	/* Domain hits. */
    };

struct pfamDat *pfamDatLoad(char **row);
/* Load a pfamDat from row fetched with select * from pfamDat
 * from database.  Dispose of this with pfamDatFree(). */

struct pfamDat *pfamDatLoadAll(char *fileName);
/* Load all pfamDat from whitespace-separated file.
 * Dispose of this with pfamDatFreeList(). */

struct pfamDat *pfamDatLoadAllByChar(char *fileName, char chopper);
/* Load all pfamDat from chopper separated file.
 * Dispose of this with pfamDatFreeList(). */

#define pfamDatLoadAllByTab(a) pfamDatLoadAllByChar(a, '\t');
/* Load all pfamDat from tab separated file.
 * Dispose of this with pfamDatFreeList(). */

struct pfamDat *pfamDatCommaIn(char **pS, struct pfamDat *ret);
/* Create a pfamDat out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new pfamDat */

void pfamDatFree(struct pfamDat **pEl);
/* Free a single dynamically allocated pfamDat such as created
 * with pfamDatLoad(). */

void pfamDatFreeList(struct pfamDat **pList);
/* Free a list of dynamically allocated pfamDat's */

void pfamDatOutput(struct pfamDat *el, FILE *f, char sep, char lastSep);
/* Print out pfamDat.  Separate fields with sep. Follow last field with lastSep. */

#define pfamDatTabOut(el,f) pfamDatOutput(el,f,'\t','\n');
/* Print out pfamDat as a line in a tab-separated file. */

#define pfamDatCommaOut(el,f) pfamDatOutput(el,f,',',',');
/* Print out pfamDat as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

struct pfamDat *pfamDatFromPfamFile(char *fileName);
/* Parse a hmmpfam generated file and return a pfamDat structure 
   containing hits. */

void pfamDatWritePfamFile(struct pfamDat *pfamDat, char *fileName);
/* Write out a pfam file, mimicking pfam file format enough to diff. */

#endif /* PFAMDAT_H */

