import matplotlib.pyplot as plt
import numpy as np

from ACVU_model_gillespie_lib_noYFP import ACVU_simulation
from ACVU_model_analysis_tools_lib import get_species_data_from_sim_data#, get_T_decision

save_species_list = ['M','D']#,'MY','DY']

# initialize simulation from parameter file    
#sim = ACVU_simulation('Basic_model.param')
#sim = ACVU_simulation('Notch_dependent_degradation_model.param')
sim = ACVU_simulation('Model0.param')

# initialize reactant levels N, using initial conditions <IC> for each cell
N=sim.initialize_reactants([{'Oc':1},{'Oc':1}])

np.random.seed(10)
sim_data=sim.run_full_simulation(N,save_species_list)

# plot sim_data

(t,D_vs_t)=get_species_data_from_sim_data(sim_data,'D',save_species_list)

plt.subplot(211)
for j in [0,1,2]:
    plt.plot(t[j],D_vs_t[j][:,0],'--b')
    plt.plot(t[j],D_vs_t[j][:,1],'--r')

plt.subplot(212)
# calculate and plot difference in D_i(t) for all three stages 
col='kcm'
diff_D=[]        
for j in [0,1,2]:
    diff_D.append( (D_vs_t[j][:,0] - D_vs_t[j][:,1]) / (D_vs_t[j][:,0] + D_vs_t[j][:,1]) )
    plt.plot(t[j]/60.,diff_D[j],'--'+col[j],lw=2)

plt.show()