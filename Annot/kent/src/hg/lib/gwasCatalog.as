table gwasCatalog
"NHGRI's collection of Genome-Wide Association Studies SNPs"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    string name;       "ID of SNP associated with trait"
    uint   pubMedID;   "PubMed ID of publication of the study"
    string author;     "First author of publication"
    string pubDate;    "Date of publication"
    string journal;    "Journal of publication"
    string title;      "Title of publication"
    string trait;      "Disease or trait assessed in study"
    lstring initSample; "Initial sample size"
    lstring replSample; "Replication sample size"
    string region;     "Chromosome band / region of SNP"
    string genes;      "Reported Gene(s)"
    string riskAllele; "Strongest SNP-Risk Allele"
    string riskAlFreq; "Risk Allele Frequency"
    string pValue;     "p-Value"
    string pValueDesc; "p-Value Description"
    string orOrBeta;   "Odds ratio or beta"
    string ci95;       "95% Confidence Interval"
    string platform;   "Platform and [SNPs passing QC]"
    enum(Y,N) cnv;     "Y if Copy Number Variant"
    )
