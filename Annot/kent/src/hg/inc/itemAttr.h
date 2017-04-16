/* itemAttr.h was originally generated by the autoSql program, which also 
 * generated itemAttr.c and itemAttr.sql.  This header links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2004 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#ifndef ITEMATTR_H
#define ITEMATTR_H

#define ITEMATTR_NUM_COLS 7

struct itemAttr
/* Relational object of display attributes for individual track items */
    {
    struct itemAttr *next;  /* Next in singly linked list. */
    char *name;	/* name of item */
    char *chrom;	/* chromosome */
    unsigned chromStart;	/* Start position in chromosome */
    unsigned chromEnd;	/* End position in chromosome */
    unsigned char colorR;	/* Color red component 0-255 */
    unsigned char colorG;	/* Color green component 0-255 */
    unsigned char colorB;	/* Color blue component 0-255 */
    };

void itemAttrStaticLoad(char **row, struct itemAttr *ret);
/* Load a row from itemAttr table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct itemAttr *itemAttrLoad(char **row);
/* Load a itemAttr from row fetched with select * from itemAttr
 * from database.  Dispose of this with itemAttrFree(). */

struct itemAttr *itemAttrLoadAll(char *fileName);
/* Load all itemAttr from whitespace-separated file.
 * Dispose of this with itemAttrFreeList(). */

struct itemAttr *itemAttrLoadAllByChar(char *fileName, char chopper);
/* Load all itemAttr from chopper separated file.
 * Dispose of this with itemAttrFreeList(). */

#define itemAttrLoadAllByTab(a) itemAttrLoadAllByChar(a, '\t');
/* Load all itemAttr from tab separated file.
 * Dispose of this with itemAttrFreeList(). */

struct itemAttr *itemAttrCommaIn(char **pS, struct itemAttr *ret);
/* Create a itemAttr out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new itemAttr */

void itemAttrFree(struct itemAttr **pEl);
/* Free a single dynamically allocated itemAttr such as created
 * with itemAttrLoad(). */

void itemAttrFreeList(struct itemAttr **pList);
/* Free a list of dynamically allocated itemAttr's */

void itemAttrOutput(struct itemAttr *el, FILE *f, char sep, char lastSep);
/* Print out itemAttr.  Separate fields with sep. Follow last field with lastSep. */

#define itemAttrTabOut(el,f) itemAttrOutput(el,f,'\t','\n');
/* Print out itemAttr as a line in a tab-separated file. */

#define itemAttrCommaOut(el,f) itemAttrOutput(el,f,',',',');
/* Print out itemAttr as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

struct itemAttrTbl *itemAttrTblNew(char *table);
/* Create a new itemAttr object. This saves the table name, but
 * doesn't load the data. */

void itemAttrTblLoad(struct itemAttrTbl *iat,
                     struct sqlConnection *conn,
                     char *chrom, int start, int end);
/* load itemAttrs for the specified chrom range */

struct itemAttr *itemAttrTblGet(struct itemAttrTbl *iat, char* name,
                                char *chrom, int chromStart, int chromEnd);
/* lookup an itemAttr by name and location */

void itemAttrTblFree(struct itemAttrTbl **iatPtr);
/* free an itemAttrTbl */

#endif /* ITEMATTR_H */

