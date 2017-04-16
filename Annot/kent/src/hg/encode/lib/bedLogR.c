/* bedLogR.c was originally generated by the autoSql program, which also 
 * generated bedLogR.h and bedLogR.sql.  This module links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2011 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "bedLogR.h"


void bedLogRStaticLoad(char **row, struct bedLogR *ret)
/* Load a row from bedLogR table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->chrom = row[0];
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = row[3];
ret->score = sqlUnsigned(row[4]);
safecpy(ret->strand, sizeof(ret->strand), row[5]);
ret->thickStart = sqlUnsigned(row[6]);
ret->thickEnd = sqlUnsigned(row[7]);
ret->reserved = sqlUnsigned(row[8]);
ret->logR = sqlFloat(row[9]);
}

struct bedLogR *bedLogRLoad(char **row)
/* Load a bedLogR from row fetched with select * from bedLogR
 * from database.  Dispose of this with bedLogRFree(). */
{
struct bedLogR *ret;

AllocVar(ret);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->score = sqlUnsigned(row[4]);
safecpy(ret->strand, sizeof(ret->strand), row[5]);
ret->thickStart = sqlUnsigned(row[6]);
ret->thickEnd = sqlUnsigned(row[7]);
ret->reserved = sqlUnsigned(row[8]);
ret->logR = sqlFloat(row[9]);
return ret;
}

struct bedLogR *bedLogRLoadAll(char *fileName) 
/* Load all bedLogR from a whitespace-separated file.
 * Dispose of this with bedLogRFreeList(). */
{
struct bedLogR *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[10];

while (lineFileRow(lf, row))
    {
    el = bedLogRLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct bedLogR *bedLogRLoadAllByChar(char *fileName, char chopper) 
/* Load all bedLogR from a chopper separated file.
 * Dispose of this with bedLogRFreeList(). */
{
struct bedLogR *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[10];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = bedLogRLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct bedLogR *bedLogRCommaIn(char **pS, struct bedLogR *ret)
/* Create a bedLogR out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new bedLogR */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->chrom = sqlStringComma(&s);
ret->chromStart = sqlUnsignedComma(&s);
ret->chromEnd = sqlUnsignedComma(&s);
ret->name = sqlStringComma(&s);
ret->score = sqlUnsignedComma(&s);
sqlFixedStringComma(&s, ret->strand, sizeof(ret->strand));
ret->thickStart = sqlUnsignedComma(&s);
ret->thickEnd = sqlUnsignedComma(&s);
ret->reserved = sqlUnsignedComma(&s);
ret->logR = sqlFloatComma(&s);
*pS = s;
return ret;
}

void bedLogRFree(struct bedLogR **pEl)
/* Free a single dynamically allocated bedLogR such as created
 * with bedLogRLoad(). */
{
struct bedLogR *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freez(pEl);
}

void bedLogRFreeList(struct bedLogR **pList)
/* Free a list of dynamically allocated bedLogR's */
{
struct bedLogR *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    bedLogRFree(&el);
    }
*pList = NULL;
}

void bedLogROutput(struct bedLogR *el, FILE *f, char sep, char lastSep) 
/* Print out bedLogR.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->chrom);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->chromStart);
fputc(sep,f);
fprintf(f, "%u", el->chromEnd);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->name);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->score);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->strand);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%u", el->thickStart);
fputc(sep,f);
fprintf(f, "%u", el->thickEnd);
fputc(sep,f);
fprintf(f, "%u", el->reserved);
fputc(sep,f);
fprintf(f, "%g", el->logR);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

