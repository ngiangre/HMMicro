#!/usr/bin/env python

# USAGE:
# ./sim_hmm.py $PWD
#
#
# A quick script to simulate a sequence from a 2-state Markov Model
# with known parameters

# Import libraries
#from multiprocessing.dummy import Pool as ThreadPool
from pandas import DataFrame, read_csv
import pandas as pd
import numpy as np
import sys, os

### Functions ###
def get_next_state(cur_state, state_nms, trans_mat):
    i = np.random.binomial(1, trans_mat.loc[cur_state,cur_state])
    if i == 1:
        next_state = cur_state
    else:
        next_state = state_nms[ np.logical_not(state_nms == cur_state) ][0]
    return next_state

### Get working path ###

try:
    dir_path = sys.argv[1]
except IndexError:
    print 'No directory given, using default'
    dir_path = '/Users/Felix/GitHub/HMMicro/hmm_scripts'
    

### Read parameters from file ###

# Get 2-state transition matrix
trans_mat_path = os.path.join(dir_path, 'trans_mat.tsv')
trans_mat = read_csv(trans_mat_path, sep='\t', header=0, index_col=None)
trans_mat.index = state_nms = trans_mat.columns.tolist() # Get state names as list
state_nms = np.array(state_nms)

# Get emission parameters (mean, sd)
emit_params_path = os.path.join(dir_path, 'emit_params.tsv')
emit_params = read_csv(emit_params_path, sep='\t', header=0, index_col=None)



### Create state sequence ###

np.random.seed(12345) # set seed

seq_len = 1000 # length of sequence
cur_state = state_nms[0] # starting state
state_seq = [''] * seq_len
emit_seq = [0] * seq_len
for i in range(seq_len):
    state_seq[i] = cur_state
    emit_seq[i] = np.random.normal(emit_params[cur_state][0],
                                   emit_params[cur_state][1])
    cur_state = get_next_state(cur_state, state_nms, trans_mat)
    
### Save state and emission sequences ###
state_seq_path = os.path.join(dir_path, 'state_seq.tsv')
np.savetxt(state_seq_path, np.array(state_seq), fmt='%s',delimiter='\n')

emit_seq_path = os.path.join(dir_path, 'emit_seq.tsv')
np.savetxt(emit_seq_path, np.array(emit_seq),delimiter='\n')