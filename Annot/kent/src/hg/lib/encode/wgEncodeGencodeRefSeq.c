/* wgEncodeGencodeRefSeq.c was originally generated by the autoSql program, which also 
 * generated wgEncodeGencodeRefSeq.h and wgEncodeGencodeRefSeq.sql.  This module links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2011 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "encode/wgEncodeGencodeRefSeq.h"


void wgEncodeGencodeRefSeqStaticLoad(char **row, struct wgEncodeGencodeRefSeq *ret)
/* Load a row from wgEncodeGencodeRefSeq table into ret.  The contents of ret will
 * be replaced at the next call to this function. */
{

ret->transcriptId = row[0];
ret->rnaAcc = row[1];
ret->pepAcc = row[2];
}

struct wgEncodeGencodeRefSeq *wgEncodeGencodeRefSeqLoad(char **row)
/* Load a wgEncodeGencodeRefSeq from row fetched with select * from wgEncodeGencodeRefSeq
 * from database.  Dispose of this with wgEncodeGencodeRefSeqFree(). */
{
struct wgEncodeGencodeRefSeq *ret;

AllocVar(ret);
ret->transcriptId = cloneString(row[0]);
ret->rnaAcc = cloneString(row[1]);
ret->pepAcc = cloneString(row[2]);
return ret;
}

struct wgEncodeGencodeRefSeq *wgEncodeGencodeRefSeqLoadAll(char *fileName) 
/* Load all wgEncodeGencodeRefSeq from a whitespace-separated file.
 * Dispose of this with wgEncodeGencodeRefSeqFreeList(). */
{
struct wgEncodeGencodeRefSeq *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[3];

while (lineFileRow(lf, row))
    {
    el = wgEncodeGencodeRefSeqLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct wgEncodeGencodeRefSeq *wgEncodeGencodeRefSeqLoadAllByChar(char *fileName, char chopper) 
/* Load all wgEncodeGencodeRefSeq from a chopper separated file.
 * Dispose of this with wgEncodeGencodeRefSeqFreeList(). */
{
struct wgEncodeGencodeRefSeq *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[3];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = wgEncodeGencodeRefSeqLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct wgEncodeGencodeRefSeq *wgEncodeGencodeRefSeqCommaIn(char **pS, struct wgEncodeGencodeRefSeq *ret)
/* Create a wgEncodeGencodeRefSeq out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new wgEncodeGencodeRefSeq */
{
char *s = *pS;

if (ret == NULL)
    AllocVar(ret);
ret->transcriptId = sqlStringComma(&s);
ret->rnaAcc = sqlStringComma(&s);
ret->pepAcc = sqlStringComma(&s);
*pS = s;
return ret;
}

void wgEncodeGencodeRefSeqFree(struct wgEncodeGencodeRefSeq **pEl)
/* Free a single dynamically allocated wgEncodeGencodeRefSeq such as created
 * with wgEncodeGencodeRefSeqLoad(). */
{
struct wgEncodeGencodeRefSeq *el;

if ((el = *pEl) == NULL) return;
freeMem(el->transcriptId);
freeMem(el->rnaAcc);
freeMem(el->pepAcc);
freez(pEl);
}

void wgEncodeGencodeRefSeqFreeList(struct wgEncodeGencodeRefSeq **pList)
/* Free a list of dynamically allocated wgEncodeGencodeRefSeq's */
{
struct wgEncodeGencodeRefSeq *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    wgEncodeGencodeRefSeqFree(&el);
    }
*pList = NULL;
}

void wgEncodeGencodeRefSeqOutput(struct wgEncodeGencodeRefSeq *el, FILE *f, char sep, char lastSep) 
/* Print out wgEncodeGencodeRefSeq.  Separate fields with sep. Follow last field with lastSep. */
{
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->transcriptId);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->rnaAcc);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", el->pepAcc);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

