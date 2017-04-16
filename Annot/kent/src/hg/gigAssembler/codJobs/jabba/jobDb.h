/* jobDb.h was originally generated by the autoSql program, which also 
 * generated jobDb.c and jobDb.sql.  This header links the database and
 * the RAM representation of objects. */

#ifndef JOBDB_H
#define JOBDB_H

struct submission
/* Keeps track of a job submission */
    {
    struct submission *next;  /* Next in singly linked list. */
    char *id;	/* Submission ID from scheduler */
    char *errFile;	/* Error file associated with submission */
    char *outFile;	/* Output file associated with submission */
    float cpuTime;	/* CPU time in seconds */
    char *submitTime;	/* Time submitted */
    char *startTime;	/* Start time of job */
    char *endTime;	/* End time of job */
    int retVal;	/* Return value of job */
    unsigned char gotRetVal;	/* True if got return value */
    unsigned char submitError;	/* An error occurred submitting it */
    unsigned char inQueue;	/* Currently in queuing system */
    unsigned char queueError;	/* In error stat in queue */
    unsigned char trackingError;	/* Have lost track of this somehow - no output, not on queue */
    unsigned char running;	/* Currently running */
    unsigned char crashed;	/* Looks like it ran but crashed */
    unsigned char slow;	/* Run so long we warn user */
    unsigned char hung;	/* Run so long we kill it */
    unsigned char ranOk;	/* Looks like it ran and finished ok */
    };

struct submission *submissionCommaIn(char **pS, struct submission *ret);
/* Create a submission out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new submission */

void submissionFree(struct submission **pEl);
/* Free a single dynamically allocated submission such as created
 * with submissionLoad(). */

void submissionFreeList(struct submission **pList);
/* Free a list of dynamically allocated submission's */

void submissionOutput(struct submission *el, FILE *f, char sep, char lastSep);
/* Print out submission.  Separate fields with sep. Follow last field with lastSep. */

#define submissionTabOut(el,f) submissionOutput(el,f,'\t','\n');
/* Print out submission as a line in a tab-separated file. */

#define submissionCommaOut(el,f) submissionOutput(el,f,',',',');
/* Print out submission as a comma separated list including final comma. */

struct check
/* How to check a job */
    {
    struct check *next;  /* Next in singly linked list. */
    char *when;	/* When to check - currently either 'in' or 'out' */
    char *what;	/* What to check - 'exists' 'nonzero' 'lastLine' */
    char *file;	/* File to check */
    };

struct check *checkCommaIn(char **pS, struct check *ret);
/* Create a check out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new check */

void checkFree(struct check **pEl);
/* Free a single dynamically allocated check such as created
 * with checkLoad(). */

void checkFreeList(struct check **pList);
/* Free a list of dynamically allocated check's */

void checkOutput(struct check *el, FILE *f, char sep, char lastSep);
/* Print out check.  Separate fields with sep. Follow last field with lastSep. */

#define checkTabOut(el,f) checkOutput(el,f,'\t','\n');
/* Print out check as a line in a tab-separated file. */

#define checkCommaOut(el,f) checkOutput(el,f,',',',');
/* Print out check as a comma separated list including final comma. */

struct job
/* Keeps track of a job */
    {
    struct job *next;  /* Next in singly linked list. */
    char *command;	/* Command line for job */
    int checkCount;	/* Count of checks */
    struct check *checkList;	/* Ways to check success of job. */
    int submissionCount;	/* The number of times submitted */
    struct submission *submissionList;	/* List of submissions */
    char *spec;	/* Specification for job */
    };

struct job *jobCommaIn(char **pS, struct job *ret);
/* Create a job out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new job */

void jobFree(struct job **pEl);
/* Free a single dynamically allocated job such as created
 * with jobLoad(). */

void jobFreeList(struct job **pList);
/* Free a list of dynamically allocated job's */

void jobOutput(struct job *el, FILE *f, char sep, char lastSep);
/* Print out job.  Separate fields with sep. Follow last field with lastSep. */

#define jobTabOut(el,f) jobOutput(el,f,'\t','\n');
/* Print out job as a line in a tab-separated file. */

#define jobCommaOut(el,f) jobOutput(el,f,',',',');
/* Print out job as a comma separated list including final comma. */

struct jobDb
/* Keeps track of a batch of jobs.  */
    {
    struct jobDb *next;  /* Next in singly linked list. */
    int jobCount;	/* The number of total jobs */
    struct job *jobList;	/* List of all jobs */
    };

struct jobDb *jobDbCommaIn(char **pS, struct jobDb *ret);
/* Create a jobDb out of a comma separated string. 
 * This will fill in ret if non-null, otherwise will
 * return a new jobDb */

void jobDbFree(struct jobDb **pEl);
/* Free a single dynamically allocated jobDb such as created
 * with jobDbLoad(). */

void jobDbFreeList(struct jobDb **pList);
/* Free a list of dynamically allocated jobDb's */

void jobDbOutput(struct jobDb *el, FILE *f, char sep, char lastSep);
/* Print out jobDb.  Separate fields with sep. Follow last field with lastSep. */

#define jobDbTabOut(el,f) jobDbOutput(el,f,'\t','\n');
/* Print out jobDb as a line in a tab-separated file. */

#define jobDbCommaOut(el,f) jobDbOutput(el,f,',',',');
/* Print out jobDb as a comma separated list including final comma. */

#endif /* JOBDB_H */

