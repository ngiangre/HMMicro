#!/bin/bash

work_dir="/Users/Felix/Google_Drive/HMM_Project/CLIP_data/HEK293"

bed_all="${work_dir}/Hafner.combined_AGO.proc.bed"
bed_end="${work_dir}/Hafner.combined_AGO.proc.start.end.bed"
bed_t2c="${work_dir}/Hafner.combined_AGO.proc.start.t2c_pos.bed"

out_fn="${bed_all%.*}.hg19.bed"

md5_1=`cut -f 4 ${bed_all} | md5sum`
md5_2=`cut -f 4 ${bed_end} | md5sum`
md5_3=`cut -f 4 ${bed_t2c} | md5sum`

if [ ! "${md5_1}" == "${md5_2}" ] && [ ! "${md5_1}" == "${md5_3}" ]; then
	echo "files are not in same sort order"
	exit 1
fi

paste <( cut -f 1 ${bed_all} ) \
	<( cut -f 2-4 ${bed_end} ) \
	<( cut -f 5 ${bed_all} ) \
	<( cut -f 3 ${bed_t2c} ) \
	<( cut -f 7 ${bed_all} ) > ${out_fn}