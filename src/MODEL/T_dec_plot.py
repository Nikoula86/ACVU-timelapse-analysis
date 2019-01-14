import matplotlib.pyplot as plt
import pickle as pickle
import os
import numpy as np

from ACVU_model_analysis_tools_lib import get_run_data_from_sweep_data, get_species_data_from_sim_data, get_T_decision

import matplotlib as mpl

mpl.rc('font',family='Arial')
mpl.rcParams['pdf.fonttype'] = 42

def plot_Tdec(model, color, fig, species):

    # sweep parameters
    save_species_list = ['M','D','MY','DY']
    sweep_data_dir = model + '_sweep'
    sweep_data_file = 'sweep'
    sweep_param_lbl = 'DT'
    sweep_param_list = [0.,10.,20.,30.,40.,50.,60.,70.,80.,90.]

    ### analyze sweep data

    # load sweep data
    infile=os.path.join(os.getcwd(), sweep_data_dir, sweep_data_file + '_' + sweep_param_lbl + '.p' )
    sweep_data = pickle.load( open( infile, "rb" ) )

    # initialize matrix for AC_id and T_dec
    AC_ids = np.array( [ np.zeros((len(i['run_data']))) for i in sweep_data ] )
    T_decs = np.array( [ np.zeros((len(i['run_data']))) for i in sweep_data ] )

    # for each sweep_parameter
    for i, sweep_param in enumerate( sweep_param_list ):
        print(sweep_param)

        # get run data for particular parameter value
        run_data=get_run_data_from_sweep_data(sweep_data,sweep_param)

        # extract AC_id and T_dec
        # for each simulation in run_data
        for j in range(0,len(run_data)):
            # get species data
            (t,ts)=get_species_data_from_sim_data(run_data[j]['sim_data'],species,save_species_list)
            # print(t)
            # break
            # get T_decision
            (AC_id,T_dec, ind)=get_T_decision(t,ts,0.3)
            # print(AC_id, T_dec)

            AC_ids[i,j] = AC_id
            T_decs[i,j] = T_dec

    frac = np.array( [ np.sum(AC_id)/AC_id.shape[0] for AC_id in AC_ids ] )
    mean_T_dec = np.array( [ np.mean(T_dec) for T_dec in T_decs ] )
    std_T_dec = np.array( [ np.std(T_dec) for T_dec in T_decs ] )

    print(mean_T_dec)
    print(np.sum(np.array(T_decs[0])>30)/len(T_decs[0]))
    print(T_decs.shape)

    [ ax1, ax2, ax3 ] = fig.get_axes()
    ax1line = ax1.plot(sweep_param_list, frac, '-o', ms = 3, color = color, lw = 0.5)
    ax2line = ax2.plot(sweep_param_list, mean_T_dec,'-o', ms = 3, color = color, lw = 0.5)
    ax3line = ax3.hist(T_decs[0], 
                        bins = 48,
                        range = [0,480],
                        orientation='horizontal', 
                        weights = [ 1./len(T_decs[0]) for i in np.arange(len(T_decs[0])) ],
                        color = color, alpha=.5)

    return ax1line, ax2line, ax3line

def setup_figure():

    # plt results
    fig = plt.figure( figsize = (6,3) )
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    # birth order bias plot
    ax1.set_xlim(-10,100)
    ax1.set_ylim(0.05,1.05)
    ax1.set_yticks([0.5,0.6,0.7,0.8,0.9,1.0])

    ax1.plot([-10,100],[0.5,0.5],'--k',lw=0.25)    
    ax1.plot([-10,100],[1.,1.],'--k',lw=0.25)    

    ax1.set_ylabel('Fraction of animals (%)', fontsize = 7, fontname = 'Arial',labelpad=0)
    ax1.set_xlabel('Time between\n birth (min)', fontsize = 7, fontname = 'Arial',labelpad=0)

    ax1.set_xticks([0,30,60,90])

    # time to decision plot
    ax2.set_xlim(-10,100)
    ax2.set_ylim(-10.,450)

    ax2.set_ylabel('Time to decision (min)', fontsize = 7, fontname = 'Arial',labelpad=0)
    ax2.set_xlabel('Time between\n birth (min)', fontsize = 7, fontname = 'Arial',labelpad=0)

    ax2.set_xticks([0,30,60,90])
    ax2.set_yticks([0,120,240,360])

    # histogram of Tdec for dT=0 plot
    ax3.set_ylim(0,450)
    ax3.set_xlim(0.,.5)

    ax3.set_xlabel('Fraction of animals', fontsize = 7, fontname = 'Arial',labelpad=0)
    ax3.set_ylabel('Time to decision (min)', fontsize = 7, fontname = 'Arial',labelpad=0)
    
    ax3.set_xticks([0,.5])
    ax3.set_yticks([0,120,240,360])

    fig.subplots_adjust(left=0.15, right=.95, top=.9, bottom=0.2, wspace=.5)

    for ax in fig.get_axes():
        for tl in ax.get_xticklabels():
            tl.set_fontsize(7)
        for tl in ax.get_yticklabels():
            tl.set_fontsize(7)

        
        # remove label axis and so on
        ax.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left='off',      # ticks along the bottom edge are off
            right='off',         # ticks along the top edge are off
            labelleft='off', # labels along the bottom edge are off
            length = 1.5,
            pad = 0 )
        ax.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off', # labels along the bottom edge are off
            length = 1.5,
            pad = 0 )
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(0.25)
        
        ax.xaxis.set_tick_params(width=.25)
        ax.yaxis.set_tick_params(width=.25)
        
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        
        # ax.set_yticks([0,5,10])
        
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontname('Arial')
            label.set_fontsize(7)

    return fig

models = [ 
#            [ 'Model0', 'M0', 'black', 'D' ],
#            [ 'Model0_original', 'M0_orig', 'green', 'D' ],
#            [ 'Model1', 'M1-daughters', 'gray', 'D' ],
#            [ 'Model1_original', 'M1origin-daughters', 'green', 'D' ],
            [ 'Model2_5k', 'M2-NotchDep', 'orange', 'D' ],
#            [ 'Model2', 'M2-NotchDep', 'orange', 'DY' ],
#            [ 'Model3', 'M3-lowHigh', 'green', 'D' ],
#            [ 'Model3', 'M3-lowHigh', 'green', 'DY' ],
#            [ 'Model4', 'M4-NotchInBothDaughters', 'magenta', 'D' ],
#            [ 'Model5', 'M5-NotchDep-NOBURSTS', 'red', 'D' ],
#            [ 'Model5', 'M5-NotchDep-NOBURSTS', 'red', 'DY' ],
#             [ 'Model6', 'M6-lowHIGH-NOBURSTS', 'blue', 'D' ],
#             [ 'Model6', 'M6-lowHIGH-NOBURSTS', 'blue', 'DY' ],
#            ['Model7', 'M7', 'black', 'D']
	  ]
fig = setup_figure()
#fig.suptitle('Model7 comparison - YFP dynamics')
for mod in models:
    l1,l2,l3 = plot_Tdec(mod[0],  mod[2], fig, mod[3])
    
fig.get_axes()[1].legend([m[1] for m in models], fontsize = 5)
#fig.savefig('compareStats06_bDN,H1_bDN,Nobur,H1.pdf')
plt.show()



