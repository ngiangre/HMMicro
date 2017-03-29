#!/bin/bash

#work_dir="/Users/Felix/Google_Drive/HMM_Project/CLIP_data/HEK293"
in_fn="$1"
# include_head="$2"

out_fn1="${in_fn%.*}.start.end.bed"
out_fn2="${in_fn%.*}.start.t2c_pos.bed"

# create bed with start and end fields to be lifted over
awk 'BEGIN { OFS="\t" } { if( NR > 1) { print $2, $3, $4, $1 } } ' ${in_fn} > ${out_fn1}

# create bed with start and t2c_pos fields to be lifted over
awk 'BEGIN { OFS="\t" } { if( NR > 1) { print $2, $3, $5, $1 } } ' ${in_fn} > ${out_fn2}



# if [ -z ${include_head} ] || [ ${include_head} -eq 0 ]; then
# 	awk 'BEGIN { OFS="\t" } { if( NR > 1) { print $2, $3, $4, $1, $6, $5, $7; } } ' ${in_fn} > ${out_fn}
# elif [ ${include_head} > 0 ]; then
# 	printf "chrom\tchromStart\tchromEnd\tname\tstrand\tt2c_pos\ttranscript_alignment\n" > ${out_fn}
# 	awk 'BEGIN { OFS="\t" } { if( NR > 1) { print $2, $3, $4, $1, $6, $5, $7; } } ' ${in_fn} >> ${out_fn}
# fi
