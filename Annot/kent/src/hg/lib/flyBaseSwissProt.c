/* flyBaseSwissProt.c was originally generated by the autoSql program, which also 
 * generated flyBaseSwissProt.h and flyBaseSwissProt.sql.  This module links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2014 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "flyBaseSwissProt.h"


void flyBaseSwissProtStaticLoad(char **row, struct flyBaseSwissProt *ret)
/* Load a row from flyBaseSwissProt table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->flyBaseId = row[0];
ret->swissProtId = row[1];
ret->spGeneName = row[2];
ret->spSymbol = row[3];
}

struct flyBaseSwissProt *flyBaseSwissProtLoad(char **row)
/* Load a flyBaseSwissProt from row fetched with select * from flyBaseSwissProt
 * from database.  Dispose of this with flyBaseSwissProtFree(). */
{
struct flyBaseSwissProt *ret;

AllocVar(ret);
ret->flyBaseId = cloneString(row[0]);
ret->swissProtId = cloneString(row[1]);
ret->spGeneName = cloneString(row[2]);
ret->spSymbol = cloneString(row[3]);
return ret;
}

struct flyBaseSwissProt *flyBaseSwissProtLoadAll(char *fileName) 
/* Load all flyBaseSwissProt from a whitespace-separated file.
 * Dispose of this with flyBaseSwissProtFreeList(). */
{
struct flyBaseSwissProt *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[4];

while (lineFileRow(lf, row))
    {
    el = flyBaseSwissProtLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct flyBaseSwissProt *flyBaseSwissProtLoadAllByChar(char *fileName, char chopper) 
/* Load all flyBaseSwissProt from a chopper separated file.
 * Dispose of this with flyBaseSwissProtFreeList(). */
{
struct flyBaseSwissProt *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[4];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = flyBaseSwissProtLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct flyBaseSwissProt *flyBaseSwissProtCommaIn(char **pS, struct flyBaseSwissProt *ret)
/* Create a flyBaseSwissProt out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new flyBaseSwissProt */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->flyBaseId = sqlStringComma(&s);
ret->swissProtId = sqlStringComma(&s);
ret->spGeneName = sqlStringComma(&s);
ret->spSymbol = sqlStringComma(&s);
*pS = s;
return ret;
}

void flyBaseSwissProtFree(struct flyBaseSwissProt **pEl)
/* Free a single dynamically allocated flyBaseSwissProt such as created
 * with flyBaseSwissProtLoad(). */
{
struct flyBaseSwissProt *el;

if ((el = *pEl) == NULL) return;
freeMem(el->flyBaseId);
freeMem(el->swissProtId);
freeMem(el->spGeneName);
freeMem(el->spSymbol);
freez(pEl);
}

void flyBaseSwissProtFreeList(struct flyBaseSwissProt **pList)
/* Free a list of dynamically allocated flyBaseSwissProt's */
{
struct flyBaseSwissProt *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    flyBaseSwissProtFree(&el);
    }
*pList = NULL;
}

void flyBaseSwissProtOutput(struct flyBaseSwissProt *el, FILE *f, char sep, char lastSep) 
/* Print out flyBaseSwissProt.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->flyBaseId);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->swissProtId);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->spGeneName);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->spSymbol);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

