#!/usr/bin/env python

# Import libraries
from multiprocessing.dummy import Pool as ThreadPool
from pandas import DataFrame, read_csv
import pandas as pd
import numpy as np
import sys, os

# bed_dir = '/Users/Felix/ref/hg18/transcript_features'
# bed_fn = os.path.join(bed_dir, 'hg18.' + feat + '.bed')

work_dir = '/Users/Felix/Google_Drive/HMM_Project/CLIP_data/HEK293'
tsv_fn = os.path.join(work_dir, 'Hafner.combined_AGO.raw.txt')

tsv_df = read_csv(tsv_fn, sep='\t', index_col=None, header=None, na_values='nan')

eve = [ i*2   for i in range(tsv_df.shape[0]/2) ]
odd = [ i*2+1 for i in range(tsv_df.shape[0]/2) ]

out_cols = ['ccr_id', 'chr', 'start', 'end', 't2c_pos', 'strand', 'transcript_alignment']
out_df = DataFrame(index=range(len(eve)), columns=out_cols)

for i in range(len(eve)):
    cur_elem = tsv_df.iloc[eve[i],0].split('|')
    out_df.loc[i,'ccr_id'] = cur_elem[1]
    out_df.loc[i,'t2c_pos'] = int(cur_elem[7].split(':')[1])

    coord_elem = cur_elem[6].split(':')
    
    chr_strand = coord_elem[0].split('=')[1]
    out_df.loc[i,'chr'] = chr_strand[:-1]
    out_df.loc[i,'strand'] = chr_strand[-1]

    stt_end = coord_elem[1].split('_')
    out_df.loc[i,'start'] = int(stt_end[0])
    out_df.loc[i,'end'] = int(stt_end[1])

    out_df.loc[i,'transcript_alignment'] = tsv_df.iloc[odd[i],0]

out_path = os.path.join(work_dir, 'Hafner.combined_AGO.proc.tsv')
out_df.to_csv(out_path, index=False, header=True, sep='\t', na_rep='nan')
