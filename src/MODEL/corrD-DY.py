#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 10:03:09 2018

@author: ngritti
"""

import numpy as np
import pickle
import matplotlib.pyplot as plt
from ACVU_model_analysis_tools_lib import get_run_data_from_sweep_data, get_species_data_from_sim_data, get_T_decision
from ACVU_model_analysis_tools_lib import plot_single_trajectory, setup_trajectory_figure
import warnings
from tqdm import tqdm
warnings.filterwarnings('ignore')

# load sweep_data
model = 'Model8'
sweep_data = pickle.load(open(model+'_sweep/sweep_DT.p','rb'))
save_species_list = ['M','D','MY','DY']
thr_Tdec = -1

# get run data for particular parameter value
Dt=10
run_data=get_run_data_from_sweep_data(sweep_data,Dt)

# setup figure
fig,ax=plt.subplots(figsize=(5,5),nrows=1,ncols=1)
fig.suptitle(model + '- T_dec correlation')

# select single (or multiple) trajectory index
#np.random.seed(20)
idx = np.random.randint(0,len(run_data),1)
idx = np.arange(len(run_data))

# for the selected trajectories
T_dec=[[],[]]
for i in tqdm(idx):
    # get species data
    (t,ts)=get_species_data_from_sim_data(run_data[i]['sim_data'],'D',save_species_list)
    # get T_decision
    (AC_id,td, ind)=get_T_decision(t,ts,0.3)
    T_dec[0].append(td)
    # get species data
    (t,ts)=get_species_data_from_sim_data(run_data[i]['sim_data'],'DY',save_species_list)
    # get T_decision
    (AC_id,td, ind)=get_T_decision(t,ts,0.3)
    T_dec[1].append(td)
    
ax.plot(T_dec[0],T_dec[1],'ok')
ax.set_ylim(0,360)
ax.set_xlim(0,360)
ax.set_xlabel('T_dec with LAG2 (min)')
ax.set_ylabel('T_dec with YFP (min)')

plt.show()
