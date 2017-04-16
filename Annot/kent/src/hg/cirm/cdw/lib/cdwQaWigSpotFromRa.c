/* Parser to read in a cdwQaWigSpot from a ra file where tags in ra file correspond to fields in a
 * struct. This program was generated by raToStructGen. */

/* Copyright (C) 2014 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "obscure.h"
#include "sqlNum.h"
#include "sqlList.h"
#include "ra.h"
#include "raToStruct.h"
#include "cdw.h"
#include "cdwQaWigSpotFromRa.h"

struct raToStructReader *cdwQaWigSpotRaReader()
/* Make a raToStructReader for cdwQaWigSpot */
{
static char *fields[] = {
    "spotRatio",
    "enrichment",
    "basesInGenome",
    "basesInSpots",
    "sumSignal",
    "spotSumSignal",
    };
static char *requiredFields[] = {
    "spotRatio",
    "enrichment",
    "basesInGenome",
    "basesInSpots",
    "sumSignal",
    "spotSumSignal",
    };
return raToStructReaderNew("cdwQaWigSpot", ArraySize(fields), fields, ArraySize(requiredFields), requiredFields);
}


struct cdwQaWigSpot *cdwQaWigSpotFromNextRa(struct lineFile *lf, struct raToStructReader *reader)
/* Return next stanza put into an cdwQaWigSpot. */
{
enum fields
    {
    spotRatioField,
    enrichmentField,
    basesInGenomeField,
    basesInSpotsField,
    sumSignalField,
    spotSumSignalField,
    };
if (!raSkipLeadingEmptyLines(lf, NULL))
    return NULL;

struct cdwQaWigSpot *el;
AllocVar(el);

bool *fieldsObserved = reader->fieldsObserved;
bzero(fieldsObserved, reader->fieldCount);

char *tag, *val;
while (raNextTagVal(lf, &tag, &val, NULL))
    {
    struct hashEl *hel = hashLookup(reader->fieldIds, tag);
    if (hel != NULL)
        {
	int id = ptToInt(hel->val);
	if (fieldsObserved[id])
	     errAbort("Duplicate tag %s line %d of %s\n", tag, lf->lineIx, lf->fileName);
	fieldsObserved[id] = TRUE;
	switch (id)
	    {
	    case spotRatioField:
	        {
	        el->spotRatio = sqlDouble(val);
		break;
	        }
	    case enrichmentField:
	        {
	        el->enrichment = sqlDouble(val);
		break;
	        }
	    case basesInGenomeField:
	        {
	        el->basesInGenome = sqlLongLong(val);
		break;
	        }
	    case basesInSpotsField:
	        {
	        el->basesInSpots = sqlLongLong(val);
		break;
	        }
	    case sumSignalField:
	        {
	        el->sumSignal = sqlDouble(val);
		break;
	        }
	    case spotSumSignalField:
	        {
	        el->spotSumSignal = sqlDouble(val);
		break;
	        }
	    default:
	        internalErr();
		break;
	    }
	}
    }

raToStructReaderCheckRequiredFields(reader, lf);
return el;
}

struct cdwQaWigSpot *cdwQaWigSpotLoadRa(char *fileName)
/* Return list of all cdwQaWigSpot in ra file. */
{
struct raToStructReader *reader = cdwQaWigSpotRaReader();
struct lineFile *lf = lineFileOpen(fileName, TRUE);
struct cdwQaWigSpot *el, *list = NULL;
while ((el = cdwQaWigSpotFromNextRa(lf, reader)) != NULL)
    slAddHead(&list, el);
slReverse(&list);
lineFileClose(&lf);
raToStructReaderFree(&reader);
return list;
}

struct cdwQaWigSpot *cdwQaWigSpotOneFromRa(char *fileName)
/* Return cdwQaWigSpot in file and insist there be exactly one record. */
{
struct cdwQaWigSpot *one = cdwQaWigSpotLoadRa(fileName);
if (one == NULL)
    errAbort("No data in %s", fileName);
if (one->next != NULL)
    errAbort("Multiple records in %s", fileName);
return one;
}

