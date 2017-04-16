/* bigGenePred.c was originally generated by the autoSql program, which also 
 * generated bigGenePred.h and bigGenePred.sql.  This module links the database and
 * the RAM representation of objects. */

#include "common.h"
#include "linefile.h"
#include "dystring.h"
#include "jksql.h"
#include "bigGenePred.h"



char *bigGenePredCommaSepFieldNames = "chrom,chromStart,chromEnd,name,score,strand,thickStart,thickEnd,itemRgb,blockCount,blockSizes,chromStarts,name2,cdsStartStat,cdsEndStat,exonFrames,type,geneName,geneName2,geneType";

struct bigGenePred *bigGenePredLoad(char **row)
/* Load a bigGenePred from row fetched with select * from bigGenePred
 * from database.  Dispose of this with bigGenePredFree(). */
{
struct bigGenePred *ret;

AllocVar(ret);
ret->blockCount = sqlSigned(row[9]);
ret->chrom = cloneString(row[0]);
ret->chromStart = sqlUnsigned(row[1]);
ret->chromEnd = sqlUnsigned(row[2]);
ret->name = cloneString(row[3]);
ret->score = sqlUnsigned(row[4]);
safecpy(ret->strand, sizeof(ret->strand), row[5]);
ret->thickStart = sqlUnsigned(row[6]);
ret->thickEnd = sqlUnsigned(row[7]);
ret->itemRgb = sqlUnsigned(row[8]);
{
int sizeOne;
sqlUnsignedDynamicArray(row[10], &ret->blockSizes, &sizeOne);
assert(sizeOne == ret->blockCount);
}
{
int sizeOne;
sqlUnsignedDynamicArray(row[11], &ret->chromStarts, &sizeOne);
assert(sizeOne == ret->blockCount);
}
ret->name2 = cloneString(row[12]);
ret->cdsStartStat = parseCdsStat(row[13]);
ret->cdsEndStat = parseCdsStat(row[14]);
{
int sizeOne;
sqlSignedDynamicArray(row[15], &ret->exonFrames, &sizeOne);
assert(sizeOne == ret->blockCount);
}
ret->type = cloneString(row[16]);
ret->geneName = cloneString(row[17]);
ret->geneName2 = cloneString(row[18]);
ret->geneType = cloneString(row[19]);
return ret;
}

struct bigGenePred *bigGenePredLoadAll(char *fileName) 
/* Load all bigGenePred from a whitespace-separated file.
 * Dispose of this with bigGenePredFreeList(). */
{
struct bigGenePred *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[20];

while (lineFileRow(lf, row))
    {
    el = bigGenePredLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

struct bigGenePred *bigGenePredLoadAllByChar(char *fileName, char chopper) 
/* Load all bigGenePred from a chopper separated file.
 * Dispose of this with bigGenePredFreeList(). */
{
struct bigGenePred *list = NULL, *el;
struct lineFile *lf = lineFileOpen(fileName, TRUE);
char *row[20];

while (lineFileNextCharRow(lf, chopper, row, ArraySize(row)))
    {
    el = bigGenePredLoad(row);
    slAddHead(&list, el);
    }
lineFileClose(&lf);
slReverse(&list);
return list;
}

#ifdef NOTNOW
struct bigGenePred *bigGenePredCommaIn(char **pS, struct bigGenePred *ret)
/* Create a bigGenePred out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new bigGenePred */
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
ret->itemRgb = sqlUnsignedComma(&s);
ret->blockCount = sqlSignedComma(&s);
{
int i;
s = sqlEatChar(s, '{');
AllocArray(ret->blockSizes, ret->blockCount);
for (i=0; i<ret->blockCount; ++i)
    {
    ret->blockSizes[i] = sqlSignedComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
}
{
int i;
s = sqlEatChar(s, '{');
AllocArray(ret->chromStarts, ret->blockCount);
for (i=0; i<ret->blockCount; ++i)
    {
    ret->chromStarts[i] = sqlSignedComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
}
ret->name2 = sqlStringComma(&s);
ret->cdsStartStat = parseCdsStat(&s);
ret->cdsEndStat = parseCdsStat(&s);
{
int i;
s = sqlEatChar(s, '{');
AllocArray(ret->exonFrames, ret->blockCount);
for (i=0; i<ret->blockCount; ++i)
    {
    ret->exonFrames[i] = sqlSignedComma(&s);
    }
s = sqlEatChar(s, '}');
s = sqlEatChar(s, ',');
}
ret->type = sqlStringComma(&s);
ret->geneName = sqlStringComma(&s);
ret->geneName2 = sqlStringComma(&s);
ret->geneType = sqlStringComma(&s);
*pS = s;
return ret;
}
#endif

void bigGenePredFree(struct bigGenePred **pEl)
/* Free a single dynamically allocated bigGenePred such as created
 * with bigGenePredLoad(). */
{
struct bigGenePred *el;

if ((el = *pEl) == NULL) return;
freeMem(el->chrom);
freeMem(el->name);
freeMem(el->blockSizes);
freeMem(el->chromStarts);
freeMem(el->name2);
freeMem(el->exonFrames);
freeMem(el->type);
freeMem(el->geneName);
freeMem(el->geneName2);
freeMem(el->geneType);
freez(pEl);
}

void bigGenePredFreeList(struct bigGenePred **pList)
/* Free a list of dynamically allocated bigGenePred's */
{
struct bigGenePred *el, *next;

for (el = *pList; el != NULL; el = next)
    {
    next = el->next;
    bigGenePredFree(&el);
    }
*pList = NULL;
}

void bigGenePredOutput(struct bigGenePred *el, FILE *f, char sep, char lastSep) 
/* Print out bigGenePred.  Separate fields with sep. Follow last field with lastSep. */
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
if (el->name != NULL)
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
fprintf(f, "%u", el->itemRgb);
fputc(sep,f);
fprintf(f, "%d", el->blockCount);
fputc(sep,f);
{
int i;
if (sep == ',') fputc('{',f);
for (i=0; i<el->blockCount; ++i)
    {
    fprintf(f, "%d", el->blockSizes[i]);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
}
fputc(sep,f);
{
int i;
if (sep == ',') fputc('{',f);
for (i=0; i<el->blockCount; ++i)
    {
    fprintf(f, "%d", el->chromStarts[i]);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
}
fputc(sep,f);
if (sep == ',') fputc('"',f);
if (el->name2 != NULL)
    fprintf(f, "%s", el->name2);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", genePredCdsStatStr(el->cdsStartStat));
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
fprintf(f, "%s", genePredCdsStatStr(el->cdsEndStat));
if (sep == ',') fputc('"',f);
fputc(sep,f);
{
int i;
if (sep == ',') fputc('{',f);
for (i=0; i<el->blockCount; ++i)
    {
    if (el->exonFrames != NULL)
	fprintf(f, "%d", el->exonFrames[i]);
    else
	fputs("-1", f);
    fputc(',', f);
    }
if (sep == ',') fputc('}',f);
}
fputc(sep,f);
if (sep == ',') fputc('"',f);
if (el->type != NULL)
    fprintf(f, "%s", el->type);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
if (el->geneName != NULL)
    fprintf(f, "%s", el->geneName);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
if (el->geneName2 != NULL)
    fprintf(f, "%s", el->geneName2);
if (sep == ',') fputc('"',f);
fputc(sep,f);
if (sep == ',') fputc('"',f);
if (el->geneType != NULL)
    fprintf(f, "%s", el->geneType);
if (sep == ',') fputc('"',f);
fputc(lastSep,f);
}

/* -------------------------------- End autoSql Generated Code -------------------------------- */

static char *bigGenePredAutoSqlString =
"table bigGenePred\n"
"\"bigGenePred gene models\"\n"
"   (\n"
"   string chrom;       \"Reference sequence chromosome or scaffold\"\n"
"   uint   chromStart;  \"Start position in chromosome\"\n"
"   uint   chromEnd;    \"End position in chromosome\"\n"
"   string name;        \"Name or ID of item, ideally both human readable and unique\"\n"
"   uint score;         \"Score (0-1000)\"\n"
"   char[1] strand;     \"+ or - for strand\"\n"
"   uint thickStart;    \"Start of where display should be thick (start codon)\"\n"
"   uint thickEnd;      \"End of where display should be thick (stop codon)\"\n"
"   uint reserved;       \"RGB value (use R,G,B string in input file)\"\n"
"   int blockCount;     \"Number of blocks\"\n"
"   int[blockCount] blockSizes; \"Comma separated list of block sizes\"\n"
"   int[blockCount] chromStarts; \"Start positions relative to chromStart\"\n"
"   string name2;       \"Alternative/human readable name\"\n"
"   string cdsStartStat; \"enum('none','unk','incmpl','cmpl')\"\n"
"   string cdsEndStat;   \"enum('none','unk','incmpl','cmpl')\"\n"
"   int[blockCount] exonFrames; \"Exon frame {0,1,2}, or -1 if no frame for exon\"\n"
"   string type;        \"Transcript type\"\n"
"   string geneName;    \"Primary identifier for gene\"\n"
"   string geneName2;   \"Alternative/human readable gene name\"\n"
"   string geneType;    \"Gene type\"\n"
"   )\n"
    ;

struct asObject *bigGenePredAsObj()
// Return asObject describing fields of bigGenePred
{
return asParseText(bigGenePredAutoSqlString);
}

