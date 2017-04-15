#!/usr/bin/env bash

raw_bed="$1"

out_bed="${raw_bed%.*}.tx.bed"

if [ ! -f ${out_bed} ]; then
	tail -n+2 ${raw_bed} | awk 'BEGIN { OFS="\t" } { print $3,$5,$6,$2,$4 }' > ${out_bed}
fi
