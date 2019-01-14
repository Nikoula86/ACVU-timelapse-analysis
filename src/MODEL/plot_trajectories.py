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
model = 'Model7'
sweep_data = pickle.load(open(model+'_sweep/sweep_DT.p','rb'))
save_species_list = ['M','D','MY','DY']
species = 'D'
thr_Tdec = -1

# get run data for particular parameter value
Dt=90
run_data=get_run_data_from_sweep_data(sweep_data,Dt)

# setup figure
ax = setup_trajectory_figure(Dt/60.+2)
ax[0].set_title(species+' - '+model+', Dt: %02d min, T_dec>%02d min'%(Dt,thr_Tdec))

# select single (or multiple) trajectory index
#np.random.seed(20)
idx = np.random.randint(0,len(run_data),1)
idx = np.arange(len(run_data))

# for the selected trajectories
for i in tqdm(idx[::25]):
    # get species data
    (t,ts)=get_species_data_from_sim_data(run_data[i]['sim_data'],species,save_species_list)

    # get T_decision
    (AC_id,T_dec, ind)=get_T_decision(t,ts,0.3)
    if T_dec>thr_Tdec:

        # calculate and plot difference in D_i(t) for all three stages     
        plot_single_trajectory(t,ts,ax,thr=0.3)
#        print('T_decision for trajectory %03d: %.4f'%(i,T_dec/60.))


plt.show()
