/* snpExceptions.h was originally generated by the autoSql program, which also 
 * generated snpExceptions.c and snpExceptions.sql.  This header links the database and
 * the RAM representation of objects. */

/* Copyright (C) 2013 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#ifndef SNPEXCEPTIONS_H
#define SNPEXCEPTIONS_H

#ifndef JKSQL_H
#include "jksql.h"
#endif

#define SNPEXCEPTIONS_NUM_COLS 5

struct snpExceptions
/* Set of queries to look for snps that appear problematic */
    {
    struct snpExceptions *next;  /* Next in singly linked list. */
    unsigned exceptionId;	/* unique ID for this exception */
    char *query;	/* SQL string to retrieve bad records */
    unsigned num;	/* Count of SNPs that fail this condition */
    char *description;	/* Text string for readability */
    char *resultPath;	/* path for results file */
    };

void snpExceptionsStaticLoad(char **row, struct snpExceptions *ret);
/* Load a row from snpExceptions table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct snpExceptions *snpExceptionsLoad(char **row);
/* Load a snpExceptions from row fetched with select * from snpExceptions
 * from database.  Dispose of this with snpExceptionsFree(). */

struct snpExceptions *snpExceptionsLoadAll(char *fileName);
/* Load all snpExceptions from whitespace-separated file.
 * Dispose of this with snpExceptionsFreeList(). */

struct snpExceptions *snpExceptionsLoadAllByChar(char *fileName, char chopper);
/* Load all snpExceptions from chopper separated file.
 * Dispose of this with snpExceptionsFreeList(). */

#define snpExceptionsLoadAllByTab(a) snpExceptionsLoadAllByChar(a, '\t');
/* Load all snpExceptions from tab separated file.
 * Dispose of this with snpExceptionsFreeList(). */

struct snpExceptions *snpExceptionsLoadByQuery(struct sqlConnection *conn, char *query);
/* Load all snpExceptions from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with snpExceptionsFreeList(). */

void snpExceptionsSaveToDb(struct sqlConnection *conn, struct snpExceptions *el, char *tableName, int updateSize);
/* Save snpExceptions as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. Strings are automatically escaped to allow insertion into the database. */

struct snpExceptions *snpExceptionsCommaIn(char **pS, struct snpExceptions *ret);
/* Create a snpExceptions out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new snpExceptions */

void snpExceptionsFree(struct snpExceptions **pEl);
/* Free a single dynamically allocated snpExceptions such as created
 * with snpExceptionsLoad(). */

void snpExceptionsFreeList(struct snpExceptions **pList);
/* Free a list of dynamically allocated snpExceptions's */

void snpExceptionsOutput(struct snpExceptions *el, FILE *f, char sep, char lastSep);
/* Print out snpExceptions.  Separate fields with sep. Follow last field with lastSep. */

#define snpExceptionsTabOut(el,f) snpExceptionsOutput(el,f,'\t','\n');
/* Print out snpExceptions as a line in a tab-separated file. */

#define snpExceptionsCommaOut(el,f) snpExceptionsOutput(el,f,',',',');
/* Print out snpExceptions as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* SNPEXCEPTIONS_H */

