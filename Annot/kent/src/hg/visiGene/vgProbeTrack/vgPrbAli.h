/* vgPrbAli.h was originally generated by the autoSql program, which also 
 * generated vgPrbAli.c and vgPrbAli.sql.  This header links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2013 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#ifndef VGPRBALI_H
#define VGPRBALI_H

#ifndef JKSQL_H
#include "jksql.h"
#endif

#define VGPRBALI_NUM_COLS 3

struct vgPrbAli
/* Define Probe Alignment status for track vgProbes */
    {
    struct vgPrbAli *next;  /* Next in singly linked list. */
    char *db;	/* assembly */
    int vgPrb;	/* vgPrb id */
    char *status;	/* new,ali,none */
    };

void vgPrbAliStaticLoad(char **row, struct vgPrbAli *ret);
/* Load a row from vgPrbAli table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct vgPrbAli *vgPrbAliLoad(char **row);
/* Load a vgPrbAli from row fetched with select * from vgPrbAli
 * from database.  Dispose of this with vgPrbAliFree(). */

struct vgPrbAli *vgPrbAliLoadAll(char *fileName);
/* Load all vgPrbAli from whitespace-separated file.
 * Dispose of this with vgPrbAliFreeList(). */

struct vgPrbAli *vgPrbAliLoadAllByChar(char *fileName, char chopper);
/* Load all vgPrbAli from chopper separated file.
 * Dispose of this with vgPrbAliFreeList(). */

#define vgPrbAliLoadAllByTab(a) vgPrbAliLoadAllByChar(a, '\t');
/* Load all vgPrbAli from tab separated file.
 * Dispose of this with vgPrbAliFreeList(). */

struct vgPrbAli *vgPrbAliLoadByQuery(struct sqlConnection *conn, char *query);
/* Load all vgPrbAli from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with vgPrbAliFreeList(). */

void vgPrbAliSaveToDb(struct sqlConnection *conn, struct vgPrbAli *el, char *tableName, int updateSize);
/* Save vgPrbAli as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. Strings are automatically escaped to allow insertion into the database. */

struct vgPrbAli *vgPrbAliCommaIn(char **pS, struct vgPrbAli *ret);
/* Create a vgPrbAli out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new vgPrbAli */

void vgPrbAliFree(struct vgPrbAli **pEl);
/* Free a single dynamically allocated vgPrbAli such as created
 * with vgPrbAliLoad(). */

void vgPrbAliFreeList(struct vgPrbAli **pList);
/* Free a list of dynamically allocated vgPrbAli's */

void vgPrbAliOutput(struct vgPrbAli *el, FILE *f, char sep, char lastSep);
/* Print out vgPrbAli.  Separate fields with sep. Follow last field with lastSep. */

#define vgPrbAliTabOut(el,f) vgPrbAliOutput(el,f,'\t','\n');
/* Print out vgPrbAli as a line in a tab-separated file. */

#define vgPrbAliCommaOut(el,f) vgPrbAliOutput(el,f,',',',');
/* Print out vgPrbAli as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* VGPRBALI_H */

