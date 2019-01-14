import numpy as np

def get_run_data_from_sweep_data(sweep_data,sweep_param_value):
    for s in sweep_data:
        if s['param_value']==sweep_param_value:
            return(s['run_data'])
    return([])

def get_species_data_from_sim_data(sim_data,species_lbl,species_list):
    ind=1+[x for (x,y) in enumerate(species_list) if y==species_lbl][0]
    t=[]
    S=[]
    for i in [0,1,2]:
        t.append(sim_data[i]['trajectory_data'][0])
        S.append(sim_data[i]['trajectory_data'][ind])
    return(t,S)

def get_T_decision(t,D_vs_t,h=0.3):
    # filter only last part of the simulation
    t = t[-1] - t[-1][0]
    D_vs_t = D_vs_t[-1]
    # compute ratio, always with ending in the positive region
    (AC_id,f)=(0,1) if ( D_vs_t[-1,0] - D_vs_t[-1,1] ) > 0 else (1,-1)
    diff_D = f*((D_vs_t[:,0] - D_vs_t[:,1]) / (D_vs_t[:,0] + D_vs_t[:,1]))
    # compute when threshold is passed
    r=np.where(diff_D < h)[0]
    if len(r)>0:
        ind = r[-1]
        return ((AC_id, t[ind], [-1,ind]))
    else:
        return((AC_id,t[0],[-1,0]))
    
    return None
    
def get_decision_data(sweep_data, sweep_param_value, species_list):
    run_data=get_run_data_from_sweep_data(sweep_data, sweep_param_value)

    AC_id_list=[]
    T_dec_list=[]
    for i in range(0,len(run_data)):
        # get species data
        (t,D_vs_t)=get_species_data_from_sim_data(run_data[i]['sim_data'],'D',species_list)
        # get T_decision
        (AC_id,T_dec, ind)=get_T_decision(t,D_vs_t,0.5)
        AC_id_list.append (AC_id)
        T_dec_list.append (T_dec)
        
    return( (np.array(AC_id_list),np.array(T_dec_list)) )

def plot_single_trajectory(t,ts,ax,thr=0.3):
    # plot absolute
    from matplotlib import cm
    ax[0].set_ylim(-5,5000)#np.round((np.max(ts[2])/1000)+1)*1000)
    lineages = [['1.p','1.pp','1.ppp'],['4.a','4.aa','4.aaa']]
    col = np.array( [ np.array( [cm.Blues(int(i))[:3] for i in (np.arange(len(lineages[0]))*254./(len(lineages[0])-1) / 254. ) * 127 + 127] ),
                np.array( [cm.Reds(int(i))[:3] for i in (np.arange(len(lineages[0]))*254./(len(lineages[0])-1) / 254. ) * 127 + 127 ] ) ] )
    for c in [0,1]:
        for j in [0,1,2]:
            ax[0].plot(t[j]/60.,ts[j][:,c],'-',color=col[c,j],lw=2)
            
    # plot ratio
    col='kcm'
    diff_D=[]        
    for j in [0,1,2]:
        diff_D.append( (ts[j][:,0] - ts[j][:,1]) / (ts[j][:,0] + ts[j][:,1]) )
        ax[1].plot(t[j]/60.,diff_D[j],'-'+col[j],lw=2)

    
    # get T_decision
    (AC_id,T_dec, ind)=get_T_decision(t,ts,0.3)
    # and plot in figure
    if len(ind)>0:
        ax[1].plot(t[ ind[0] ][ ind[1] ]/60., diff_D[ ind[0] ][ ind[1] ], 'ok' )
        
def setup_trajectory_figure(Dt):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize = (6,4),nrows=2,ncols=1)
    ax[0].set_ylim(-2,1000)
    ax[0].plot([0,12],[0,0],'--k')
    ax[0].set_xlabel('Time (hours)')
    ax[0].set_ylabel('Protein level')
    ax[1].set_ylim(-1.05,1.05)
    ax[1].plot([0,12],[0,0],'--k')    
    ax[1].plot([Dt,Dt],[-1,1],'--k')
    ax[1].set_xlabel('Time (hours)')
    ax[1].set_ylabel('Normalized difference')
    return ax
