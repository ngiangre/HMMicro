/* knownToSuper.c was originally generated by the autoSql program, which also 
 * generated knownToSuper.h and knownToSuper.sql.  This module links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2014 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "knownToSuper.h"


void knownToSuperStaticLoad(char **row, struct knownToSuper *ret)
/* Load a row from knownToSuper table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->gene = row[0];
ret->superfamily = sqlSigned(row[1]);
ret->start = sqlSigned(row[2]);
ret->end = sqlSigned(row[3]);
ret->eVal = atof(row[4]);
}

struct knownToSuper *knownToSuperLoad(char **row)
/* Load a knownToSuper from row fetched with select * from knownToSuper
 * from database.  Dispose of this with knownToSuperFree(). */
{
struct knownToSuper *ret;

AllocVar(ret);
ret->gene = cloneString(row[0]);
ret->superfamily = sqlSigned(row[1]);
ret->start = sqlSigned(row[2]);
ret->end = sqlSigned(row[3]);
ret->eVal = atof(row[4]);
return ret;
}

struct knownToSuper *knownToSuperLoadAll(char *fileName) 
/* Load all knownToSuper from a whitespace-separated file.
 * Dispose of this with knownToSuperFreeList(). */
{
struct knownToSuper *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[5];

while (lineFileRow(lf, row))
    {
    el = knownToSuperLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct knownToSuper *knownToSuperLoadAllByChar(char *fileName, char chopper) 
/* Load all knownToSuper from a chopper separated file.
 * Dispose of this with knownToSuperFreeList(). */
{
struct knownToSuper *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[5];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = knownToSuperLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct knownToSuper *knownToSuperCommaIn(char **pS, struct knownToSuper *ret)
/* Create a knownToSuper out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new knownToSuper */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->gene = sqlStringComma(&s);
ret->superfamily = sqlSignedComma(&s);
ret->start = sqlSignedComma(&s);
ret->end = sqlSignedComma(&s);
ret->eVal = sqlFloatComma(&s);
*pS = s;
return ret;
}

void knownToSuperFree(struct knownToSuper **pEl)
/* Free a single dynamically allocated knownToSuper such as created
 * with knownToSuperLoad(). */
{
struct knownToSuper *el;

if ((el = *pEl) == NULL) return;
freeMem(el->gene);
freez(pEl);
}

void knownToSuperFreeList(struct knownToSuper **pList)
/* Free a list of dynamically allocated knownToSuper's */
{
struct knownToSuper *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    knownToSuperFree(&el);
    }
*pList = NULL;
}

void knownToSuperOutput(struct knownToSuper *el, FILE *f, char sep, char lastSep) 
/* Print out knownToSuper.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->gene);
if (sep == ',') fputc('"',f);
fputc(sep,f);
fprintf(f, "%d", el->superfamily);
fputc(sep,f);
fprintf(f, "%d", el->start);
fputc(sep,f);
fprintf(f, "%d", el->end);
fputc(sep,f);
fprintf(f, "%f", el->eVal);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

