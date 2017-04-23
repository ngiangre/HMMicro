#!/usr/bin/env bash

ago_bed="Hafner.combined_AGO.proc.hg19.bed"
chr="chr22"
ref_bed="../../ref/refseq.hg19.${chr}.tx.bed"
out_bed="${ago_bed%.*}.refseq.${chr}.txt"

tmp_fn1="./.tmp_fn1"
tmp_fn2="./.tmp_fn2"

if [ ! -f ${out_bed} ]; then
	grep "^${chr}	" "${ago_bed}" | while read -r line ; do
		site_id=`echo "${line}" | cut -f 4`
		echo "${line}" | bedtools intersect -a ${ref_bed} -b - | cut -f 4 > ${tmp_fn2}
		n=`cat ${tmp_fn2} | wc -l`
		printf "${site_id}%.0s\n" $( seq 1 ${n} ) > ${tmp_fn1}

		paste "${tmp_fn1}" "${tmp_fn2}" >> "${out_bed}.tmp"
		rm -f ${tmp_fn1} ${tmp_fn2}
	done

	# Get rid of things that didn't match and keep only uniq pairs
	grep -v $'\t$' "${out_bed}.tmp" | sort | uniq > "${out_bed}"
fi
