table pslWQueryID
"Summary info about a patSpace alignment with a query ID addition"
    (
    uint    matches;      "Number of bases that match that aren't repeats"
    uint    misMatches;   "Number of bases that don't match"
    uint    repMatches;   "Number of bases that match but are part of repeats"
    uint    nCount;       "Number of 'N' bases"
    uint    qNumInsert;   "Number of inserts in query"
    uint     qBaseInsert;  "Number of bases inserted in query"
    uint    tNumInsert;   "Number of inserts in target"
    uint     tBaseInsert;  "Number of bases inserted in target"
    char[2] strand;       "+ or - for query strand. For mouse second +/- for genomic strand"
    string  qName;        "Query sequence name"
    uint    qSize;        "Query sequence size"
    uint    qStart;       "Alignment start position in query"
    uint    qEnd;         "Alignment end position in query"
    string  tName;        "Target sequence name"
    uint    tSize;        "Target sequence size"
    uint    tStart;       "Alignment start position in target"
    uint    tEnd;         "Alignment end position in target"
    uint    blockCount;   "Number of blocks in alignment"
    uint[blockCount] blockSizes;  "Size of each block"
    uint[blockCount] qStarts;     "Start of each block in query."
    uint[blockCount] tStarts;     "Start of each block in target."
    string    queryID;         "query ID field"
    )
