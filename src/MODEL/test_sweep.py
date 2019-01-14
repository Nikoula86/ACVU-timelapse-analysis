import matplotlib.pyplot as plt
import pickle as pickle
import os

from ACVU_model_analysis_tools_lib import get_run_data_from_sweep_data, get_species_data_from_sim_data, get_T_decision

model = 'Model6_5k'

if model == 'Model4':
    from ACVU_model_gillespie_lib_M4 import ACVU_simulation
    save_species_list = ['M','D']
elif (model == 'Model0') or (model=='Model1'):
    from ACVU_model_gillespie_lib_noYFP import ACVU_simulation
    save_species_list = ['M','D']    
else:
    from ACVU_model_gillespie_lib import ACVU_simulation
    save_species_list = ['M','D','MY','DY']

# sweep parameters
sweep_data_dir = model + '_sweep'
sweep_data_file = 'sweep'
sweep_param_lbl = 'DT'
sweep_param_list = [0.]#,10.,20.,30.,40.,50.,60.,70.,80.,90.]
N_run = 5000

### run sweep

# initialize simulation from parameter file    
sim = ACVU_simulation(model + '.param')
# and do sweep
sim.do_sweep(sweep_param_lbl,sweep_param_list,N_run,sweep_data_dir,sweep_data_file,save_species_list,coreFraction=0.8)


### analyze sweep data

# load sweep data
infile=os.path.join(os.getcwd(), sweep_data_dir, sweep_data_file + '_' + sweep_param_lbl + '.p' )
sweep_data = pickle.load( open( infile, "rb" ) )

# get run data for particular parameter value
run_data=get_run_data_from_sweep_data(sweep_data,20)

# plot dynamics and T_dec
col='kcm'
# for each simulation in run_data
for i in range(0,len(run_data))[::10]:
    # get species data
    (t,D_vs_t)=get_species_data_from_sim_data(run_data[i]['sim_data'],'D',save_species_list)
    # calculate and plot difference in D_i(t) for all three stages 
    diff_D=[]        
    for j in [0,1,2]:
        diff_D.append( (D_vs_t[j][:,0] - D_vs_t[j][:,1]) / (D_vs_t[j][:,0] + D_vs_t[j][:,1]) )
        plt.plot(t[j]/60.,diff_D[j],'-'+col[j],lw=2)
    # get T_decision
    (AC_id,T_dec, ind)=get_T_decision(t,D_vs_t,0.5)
    # and plot in figure
    if len(ind)>0:
        plt.plot(t[ ind[0] ][ ind[1] ]/60., diff_D[ ind[0] ][ ind[1] ], 'ok' )

plt.plot([0,12],[0,0],'--k')    
plt.ylim([-1.1,1.1])
plt.xlabel('Time (hours)')
plt.ylabel('Normalized difference')

plt.show()
