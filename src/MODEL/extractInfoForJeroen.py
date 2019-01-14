#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 10:48:57 2018

@author: ngritti
"""

import pickle
import os
import numpy as np

from ACVU_model_analysis_tools_lib import get_run_data_from_sweep_data, get_species_data_from_sim_data, get_T_decision

import matplotlib as mpl

mpl.rc('font',family='Arial')
mpl.rcParams['pdf.fonttype'] = 42

model = 'Model3_5k'

# sweep parameters
sweep_data_dir = model + '_sweep'
sweep_data_file = 'sweep'
sweep_param_lbl = 'DT'

# load sweep data
infile=os.path.join(os.getcwd(), sweep_data_dir, sweep_data_file + '_' + sweep_param_lbl + '.p' )
sweep_data = pickle.load( open( infile, "rb" ) )

# initialize matrix for AC_id and T_dec
data = [[],[],[]]

# get run data for particular parameter value
run_data=get_run_data_from_sweep_data(sweep_data,0)

# extract random seed
data[0] = [ run['random_seed'] for run in run_data ]

#%%
'''
extract AC_id and T_dec
'''
# for each simulation in run_data
for j in range(0,len(run_data)):
    # get species data
    (t,ts)=get_species_data_from_sim_data(run_data[j]['sim_data'],'D',['M','D'])
    # print(t)
    # break
    # get T_decision
    (AC_id,T_dec, ind)=get_T_decision(t,ts,0.3)
    # print(AC_id, T_dec)
    
    data[1].append(AC_id)
    data[2].append(T_dec)
    
data = np.array(data)

#%%
'''
filter only long trajectories
'''

thr_Tdec=360.
seeds = data[0,data[2]>thr_Tdec].astype(np.uint16)

long_data = [{'param_value':0.,'run_data':[]}]
for seed in seeds:
    long_data[0]['run_data'].append( sweep_data[0]['run_data'][seed] )

pickle.dump(long_data,open(model+'_longTrajData.p','wb'))







