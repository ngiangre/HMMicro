/* cdwStep.h was originally generated by the autoSql program, which also 
 * generated cdwStep.c and cdwStep.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef CDWSTEP_H
#define CDWSTEP_H

#include "jksql.h"
#define CDWSTEPDEF_NUM_COLS 3

extern char *cdwStepDefCommaSepFieldNames;

struct cdwStepDef
/* A step in an analysis pipeline - produces new files from old. See also cdwStepRun */
    {
    struct cdwStepDef *next;  /* Next in singly linked list. */
    unsigned id;	/* Pipeline id */
    char *name;	/* Name of this analysis step */
    char *description;	/* Description of step, about a sentence. */
    };

void cdwStepDefStaticLoad(char **row, struct cdwStepDef *ret);
/* Load a row from cdwStepDef table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct cdwStepDef *cdwStepDefLoadByQuery(struct sqlConnection *conn, char *query);
/* Load all cdwStepDef from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with cdwStepDefFreeList(). */

void cdwStepDefSaveToDb(struct sqlConnection *conn, struct cdwStepDef *el, char *tableName, int updateSize);
/* Save cdwStepDef as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. This function automatically escapes quoted strings for mysql. */

struct cdwStepDef *cdwStepDefLoad(char **row);
/* Load a cdwStepDef from row fetched with select * from cdwStepDef
 * from database.  Dispose of this with cdwStepDefFree(). */

struct cdwStepDef *cdwStepDefLoadAll(char *fileName);
/* Load all cdwStepDef from whitespace-separated file.
 * Dispose of this with cdwStepDefFreeList(). */

struct cdwStepDef *cdwStepDefLoadAllByChar(char *fileName, char chopper);
/* Load all cdwStepDef from chopper separated file.
 * Dispose of this with cdwStepDefFreeList(). */

#define cdwStepDefLoadAllByTab(a) cdwStepDefLoadAllByChar(a, '\t');
/* Load all cdwStepDef from tab separated file.
 * Dispose of this with cdwStepDefFreeList(). */

struct cdwStepDef *cdwStepDefCommaIn(char **pS, struct cdwStepDef *ret);
/* Create a cdwStepDef out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new cdwStepDef */

void cdwStepDefFree(struct cdwStepDef **pEl);
/* Free a single dynamically allocated cdwStepDef such as created
 * with cdwStepDefLoad(). */

void cdwStepDefFreeList(struct cdwStepDef **pList);
/* Free a list of dynamically allocated cdwStepDef's */

void cdwStepDefOutput(struct cdwStepDef *el, FILE *f, char sep, char lastSep);
/* Print out cdwStepDef.  Separate fields with sep. Follow last field with lastSep. */

#define cdwStepDefTabOut(el,f) cdwStepDefOutput(el,f,'\t','\n');
/* Print out cdwStepDef as a line in a tab-separated file. */

#define cdwStepDefCommaOut(el,f) cdwStepDefOutput(el,f,',',',');
/* Print out cdwStepDef as a comma separated list including final comma. */

#define CDWSTEPRUN_NUM_COLS 3

extern char *cdwStepRunCommaSepFieldNames;

struct cdwStepRun
/* A particular run of a step. */
    {
    struct cdwStepRun *next;  /* Next in singly linked list. */
    unsigned id;	/* Analysis run ID */
    unsigned stepDef;	/* Associated step definition */
    char *stepVersion;	/* Version number of step */
    };

void cdwStepRunStaticLoad(char **row, struct cdwStepRun *ret);
/* Load a row from cdwStepRun table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct cdwStepRun *cdwStepRunLoadByQuery(struct sqlConnection *conn, char *query);
/* Load all cdwStepRun from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with cdwStepRunFreeList(). */

void cdwStepRunSaveToDb(struct sqlConnection *conn, struct cdwStepRun *el, char *tableName, int updateSize);
/* Save cdwStepRun as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. This function automatically escapes quoted strings for mysql. */

struct cdwStepRun *cdwStepRunLoad(char **row);
/* Load a cdwStepRun from row fetched with select * from cdwStepRun
 * from database.  Dispose of this with cdwStepRunFree(). */

struct cdwStepRun *cdwStepRunLoadAll(char *fileName);
/* Load all cdwStepRun from whitespace-separated file.
 * Dispose of this with cdwStepRunFreeList(). */

struct cdwStepRun *cdwStepRunLoadAllByChar(char *fileName, char chopper);
/* Load all cdwStepRun from chopper separated file.
 * Dispose of this with cdwStepRunFreeList(). */

#define cdwStepRunLoadAllByTab(a) cdwStepRunLoadAllByChar(a, '\t');
/* Load all cdwStepRun from tab separated file.
 * Dispose of this with cdwStepRunFreeList(). */

struct cdwStepRun *cdwStepRunCommaIn(char **pS, struct cdwStepRun *ret);
/* Create a cdwStepRun out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new cdwStepRun */

void cdwStepRunFree(struct cdwStepRun **pEl);
/* Free a single dynamically allocated cdwStepRun such as created
 * with cdwStepRunLoad(). */

void cdwStepRunFreeList(struct cdwStepRun **pList);
/* Free a list of dynamically allocated cdwStepRun's */

void cdwStepRunOutput(struct cdwStepRun *el, FILE *f, char sep, char lastSep);
/* Print out cdwStepRun.  Separate fields with sep. Follow last field with lastSep. */

#define cdwStepRunTabOut(el,f) cdwStepRunOutput(el,f,'\t','\n');
/* Print out cdwStepRun as a line in a tab-separated file. */

#define cdwStepRunCommaOut(el,f) cdwStepRunOutput(el,f,',',',');
/* Print out cdwStepRun as a comma separated list including final comma. */

#define CDWSTEPIN_NUM_COLS 5

extern char *cdwStepInCommaSepFieldNames;

struct cdwStepIn
/* Inputs to an cdwStepRun */
    {
    struct cdwStepIn *next;  /* Next in singly linked list. */
    unsigned id;	/* Input table ID */
    unsigned stepRunId;	/* Which cdwStepRun this is associated with */
    char *name;	/* Input name within step */
    unsigned ix;	/* Inputs always potentially vectors.  Have single one with zero ix for scalar input */
    unsigned fileId;	/* Associated file. */
    };

void cdwStepInStaticLoad(char **row, struct cdwStepIn *ret);
/* Load a row from cdwStepIn table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct cdwStepIn *cdwStepInLoadByQuery(struct sqlConnection *conn, char *query);
/* Load all cdwStepIn from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with cdwStepInFreeList(). */

void cdwStepInSaveToDb(struct sqlConnection *conn, struct cdwStepIn *el, char *tableName, int updateSize);
/* Save cdwStepIn as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. This function automatically escapes quoted strings for mysql. */

struct cdwStepIn *cdwStepInLoad(char **row);
/* Load a cdwStepIn from row fetched with select * from cdwStepIn
 * from database.  Dispose of this with cdwStepInFree(). */

struct cdwStepIn *cdwStepInLoadAll(char *fileName);
/* Load all cdwStepIn from whitespace-separated file.
 * Dispose of this with cdwStepInFreeList(). */

struct cdwStepIn *cdwStepInLoadAllByChar(char *fileName, char chopper);
/* Load all cdwStepIn from chopper separated file.
 * Dispose of this with cdwStepInFreeList(). */

#define cdwStepInLoadAllByTab(a) cdwStepInLoadAllByChar(a, '\t');
/* Load all cdwStepIn from tab separated file.
 * Dispose of this with cdwStepInFreeList(). */

struct cdwStepIn *cdwStepInCommaIn(char **pS, struct cdwStepIn *ret);
/* Create a cdwStepIn out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new cdwStepIn */

void cdwStepInFree(struct cdwStepIn **pEl);
/* Free a single dynamically allocated cdwStepIn such as created
 * with cdwStepInLoad(). */

void cdwStepInFreeList(struct cdwStepIn **pList);
/* Free a list of dynamically allocated cdwStepIn's */

void cdwStepInOutput(struct cdwStepIn *el, FILE *f, char sep, char lastSep);
/* Print out cdwStepIn.  Separate fields with sep. Follow last field with lastSep. */

#define cdwStepInTabOut(el,f) cdwStepInOutput(el,f,'\t','\n');
/* Print out cdwStepIn as a line in a tab-separated file. */

#define cdwStepInCommaOut(el,f) cdwStepInOutput(el,f,',',',');
/* Print out cdwStepIn as a comma separated list including final comma. */

#define CDWSTEPOUT_NUM_COLS 5

extern char *cdwStepOutCommaSepFieldNames;

struct cdwStepOut
/* Outputs to an cdwAnalysis */
    {
    struct cdwStepOut *next;  /* Next in singly linked list. */
    unsigned id;	/* Output table ID */
    unsigned stepRunId;	/* Which cdwStepRun this is associated with */
    char *name;	/* Output name within step */
    unsigned ix;	/* Outputs always potentially vectors. Have single one with zero ix for scalar output */
    unsigned fileId;	/* Associated file. */
    };

void cdwStepOutStaticLoad(char **row, struct cdwStepOut *ret);
/* Load a row from cdwStepOut table into ret.  The contents of ret will
 * be replaced at the next call to this function. */

struct cdwStepOut *cdwStepOutLoadByQuery(struct sqlConnection *conn, char *query);
/* Load all cdwStepOut from table that satisfy the query given.  
 * Where query is of the form 'select * from example where something=something'
 * or 'select example.* from example, anotherTable where example.something = 
 * anotherTable.something'.
 * Dispose of this with cdwStepOutFreeList(). */

void cdwStepOutSaveToDb(struct sqlConnection *conn, struct cdwStepOut *el, char *tableName, int updateSize);
/* Save cdwStepOut as a row to the table specified by tableName. 
 * As blob fields may be arbitrary size updateSize specifies the approx size
 * of a string that would contain the entire query. Arrays of native types are
 * converted to comma separated strings and loaded as such, User defined types are
 * inserted as NULL. This function automatically escapes quoted strings for mysql. */

struct cdwStepOut *cdwStepOutLoad(char **row);
/* Load a cdwStepOut from row fetched with select * from cdwStepOut
 * from database.  Dispose of this with cdwStepOutFree(). */

struct cdwStepOut *cdwStepOutLoadAll(char *fileName);
/* Load all cdwStepOut from whitespace-separated file.
 * Dispose of this with cdwStepOutFreeList(). */

struct cdwStepOut *cdwStepOutLoadAllByChar(char *fileName, char chopper);
/* Load all cdwStepOut from chopper separated file.
 * Dispose of this with cdwStepOutFreeList(). */

#define cdwStepOutLoadAllByTab(a) cdwStepOutLoadAllByChar(a, '\t');
/* Load all cdwStepOut from tab separated file.
 * Dispose of this with cdwStepOutFreeList(). */

struct cdwStepOut *cdwStepOutCommaIn(char **pS, struct cdwStepOut *ret);
/* Create a cdwStepOut out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new cdwStepOut */

void cdwStepOutFree(struct cdwStepOut **pEl);
/* Free a single dynamically allocated cdwStepOut such as created
 * with cdwStepOutLoad(). */

void cdwStepOutFreeList(struct cdwStepOut **pList);
/* Free a list of dynamically allocated cdwStepOut's */

void cdwStepOutOutput(struct cdwStepOut *el, FILE *f, char sep, char lastSep);
/* Print out cdwStepOut.  Separate fields with sep. Follow last field with lastSep. */

#define cdwStepOutTabOut(el,f) cdwStepOutOutput(el,f,'\t','\n');
/* Print out cdwStepOut as a line in a tab-separated file. */

#define cdwStepOutCommaOut(el,f) cdwStepOutOutput(el,f,',',',');
/* Print out cdwStepOut as a comma separated list including final comma. */

/* -------------------------------- End autoSql Generated Code -------------------------------- */

#endif /* CDWSTEP_H */

