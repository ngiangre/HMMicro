#!/usr/bin/env bash

ago_bed="$1"
chr="chr22"
ref_bed="../../ref/refseq.hg19.${chr}.tx.bed"

bedtools intersect -a ${ref_bed} -b ${ago_bed} | \
	cut -f 4 | \
	sort | uniq > ${ago_bed%.*}.refseq.${chr}.txt
#bedtools intersect -a ${ref_bed} -b ${ago_bed} | cut -f 4 > .tmp.int

# cut -f 4 ${ref_bed} > .tmp.all

# grep -vF -f .tmp.int .tmp.all > 