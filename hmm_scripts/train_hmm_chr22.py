#!/usr/bin/env python

import os, sys, time, pickle
import numpy as np
from pandas import DataFrame, read_csv
import pandas as pd
import gc
from sklearn.decomposition import PCA
import scipy.stats as stats
from scipy.stats import kstest
# %matplotlib inline
import matplotlib
import matplotlib.pyplot as plt
from pomegranate import *

### Paths, constants, files ###
# work_dir = os.getcwd()
work_dir = '/Users/Felix/GitHub/HMMicro/hmm_scripts'

chromosome = 'chr22'
feature_dir = os.path.join(work_dir,'../features', chromosome) 
pca_path = os.path.join(feature_dir, 'feature_vecs.'+chromosome+'.PCA.tsv') 
if not 'pca_df' in vars():
    pca_df = read_csv(pca_path, sep='\t', header=None, index_col=0)

hek_dir = os.path.join(work_dir, '../CLIP_data/HEK293/')
id_path = os.path.join(hek_dir, 'Hafner.combined_AGO.proc.hg19.refseq.'+chromosome+'.txt')
if not 'id_table' in vars():
    id_table = read_csv(id_path, sep='\t', header=None, index_col=None)
    id_table.columns = ['clip_id', 'refseq_id']

ago_path = os.path.join(hek_dir, 'Hafner.combined_AGO.proc.hg19.bed')
if not 'ago_table' in vars():
    ago_bed = read_csv(ago_path, sep='\t', header=None, index_col=None)
    ago_bed.columns = ['chr', 'start', 'end', 'clip_id', 'strand', 'center', 'sequence']
ago_bed_chr22 = ago_bed.loc[ago_bed['chr'] == chromosome,:]

ref_dir = os.path.join(work_dir, '../ref/')
refseq_path = os.path.join(ref_dir, 'refseq.hg19.'+chromosome+'.bed')
if not 'refseq_bed' in vars():
    refseq_bed = read_csv(refseq_path, sep='\t', header=0, index_col=None)
    
### Helper functions ###

def time_str():
    return time.strftime("%Y%m%d.%H%M%S")

# Pass in comma-sep list of exon starts, ends, get integer positions back
def parse_exons(e_starts, e_ends):
    e_starts = e_starts.split(',')
    e_ends = e_ends.split(',')
    out_list = []
    for i in range(len(e_starts)-1): # last element is empty str ''
        out_list.extend(range(int(e_starts[i])+1, int(e_ends[i])+1))
    return np.array(out_list)

# Pass in a refseq id, return the exonic positions
def get_exon_positions(ref_id, refseq_bed):
    tmp_bool = refseq_bed['name'] == ref_id
    if not np.any(tmp_bool):
        print "id not found"
        return None
    subset_df = refseq_bed.loc[tmp_bool,:]
    out_list = [None] * subset_df.shape[0]
    for i in range(subset_df.shape[0]):
        tmp_exons = parse_exons(subset_df['exonStarts'].iloc[i], subset_df['exonEnds'].iloc[i])
        if subset_df['strand'].iloc[i] == '+':
            out_list[i] = tmp_exons
        else: # Reverse the order of positions if reverse strand
            out_list[i] = np.flipud( tmp_exons )
    return out_list # list of np arrays

# Pass in a clip_id
def get_bind_positions(clip_id, ago_bed):
    tmp_bool = ago_bed['clip_id'] == clip_id
    if not np.any(tmp_bool):
        print "id not found"
        return None
    row = ago_bed.loc[tmp_bool,:]
    out_array = np.array( range(row['start'], row['end'] + 1) )
    if row['strand'].iloc[0] == '-':
        out_array = np.flipud( out_array )
    return out_array

### Get the training and testing data into list of arrays of PC1 value ###
all_ind = []
all_pos = []
all_val = []
all_sts = []
for i, row in id_table.iterrows():
    tmp_pos = get_exon_positions(row['refseq_id'], refseq_bed)
    for arr in tmp_pos:
        cur_val = pca_df.loc[arr,1].values
        clip_pos = get_bind_positions(row['clip_id'], ago_bed_chr22)
        st = np.in1d(arr, clip_pos).astype(int)
        if (not np.any(np.isnan(cur_val))) and (max(st) == 1):
            all_pos.append( arr )
            all_val.append( cur_val )
            all_ind.append( i )
            all_sts.append( st )
all_ind = np.array(all_ind)
all_pos = np.array(all_pos)
all_val = np.array(all_val)
all_sts = np.array(all_sts)
gc.collect()

### Construct k validation partitions ###
k_parts = 5
num_seqs = len(all_ind)
valid_seed = 9012
np.random.seed(valid_seed)
shuf = np.array( range(len(all_ind)) )
np.random.shuffle(shuf)
# shuf_ind = all_ind[shuf]
# shuf_pos = all_pos[shuf]
# shuf_val = all_val[shuf]
# shuf_sts = all_sts[shuf]

partition_list = [None] * k_parts
parts_size = num_seqs / k_parts
for i in range(k_parts):
    if (i+1) == k_parts:
        cur_inds = np.split( shuf, [i*parts_size, num_seqs] )
    else:
        cur_inds = np.split( shuf, [i*parts_size, (i+1)*parts_size] )
    test_inds = cur_inds[1]
    train_inds = np.concatenate( (cur_inds[0],cur_inds[2]) )
    partition_list[i] = {'test':test_inds, 'train':train_inds} 

### Construct HMM ###
# num_fits = 5
model_list = [None] * k_parts

cur_fn = 'fitted.models.chr22.k_{}.seed_{}.p'.format(k_parts, valid_seed)
cur_path = os.path.join(work_dir, cur_fn)

if os.path.exists(cur_path):
    model_list = pickle.load( open( cur_path, "rb" ) )
else:
    for i in range(k_parts):
        # randomly select some starting means and std devs
        cur_mus = np.random.uniform(-1,1, size=2)
        cur_sds = np.random.uniform(0,2, size=2)

        cur_st_trans = np.random.random_sample(3)
        # Build up HMM
        s1 = State( NormalDistribution(cur_mus[0],cur_sds[0]), name='A' ) # non-binding?
        s2 = State( NormalDistribution(cur_mus[1],cur_sds[1]), name='B' ) # binding?
        model = HiddenMarkovModel()
        model.add_states( [s1, s2] )
        model.add_transition( model.start, s1, cur_st_trans[0])
        model.add_transition( model.start, s2, 1-cur_st_trans[0])
        model.add_transition( s1, s1, cur_st_trans[1])
        model.add_transition( s1, s2, 1-cur_st_trans[1])
        model.add_transition( s2, s2, cur_st_trans[2])
        model.add_transition( s2, s1, 1-cur_st_trans[2])
        model.bake()
        # Fit HMM
        cur_train = all_val[partition_list[i]['train']]
        model.fit(cur_train, distribution_inertia=0.1, edge_inertia=0.1)
        model_list[i] = model

    # Save the results
    pickle.dump( model_list, open( cur_path, "wb" ) )
    
### Print transition matrices, emission params ###
for m in model_list:
    print m.dense_transition_matrix()
    
for m in model_list:
    for st in m.states:
        if st.name in ['A','B']:
            print(st.name)
            print(st.distribution)

# Since all models are about the same, we just take one of them
model = model_list[0]
# print model.dense_transition_matrix()
# print model.states

### Viterbi testing
# fb_states = []
test_results = []
for i in range(k_parts):
    model = model_list[i]
    vit_error = []
    vit_states = []
    vit_logP = []
    cur_test_val = all_val[partition_list[i]['test']]
    cur_test_sts = all_sts[partition_list[i]['test']]
    for j in range(len(cur_test_val)):
        val_array = cur_test_val[j]
        sts_array = cur_test_sts[j]
        
        vit = model.viterbi(val_array) # Viterbi
        arr = np.array( [ x[0] for x in vit[1] ] )
        # fb  = np.array( [ x[0] for x in model.predict(val_array, algorithm = 'map')[1] ] ) # Fwd-bwd        
        vit_states.append(arr[1:]) # leave out initialization state
        vit_logP.append(vit[0])
        # fb_states.append(fb)
        cur_error = sum( np.abs(arr[1:] - sts_array) ) # abs error
        # cur_error = sum( (arr[1:] - sts_array) ** 2 ) # sq error
        vit_error.append(cur_error)
    test_results.append({'logP':np.array(vit_logP), 'states':np.array(vit_states),
                         'error':np.array(vit_error)})

gc.collect()

### Messing around
for x in test_results:
    print np.median(x['error'])

