import numpy as np
import pickle as pickle
import os as os
import time
from multiprocessing import Pool



class ACVU_simulation:
    
    def __init__(self,param_file):
        
        # number of species and reactions per cell
        self.N_S=10
        self.N_R=18

        # species labels and indeces
        self.species_lbl=['O','Oc','OH','OcH','M','D','OY','OYH','MY','DY']
        self.species_ind={'O':0,'Oc':1,'OH':2,'OcH':3,'M':4,'D':5,'OY':6,'OYH':7,'MY':8,'DY':9}

        self.param_list=self.read_param_file(param_file)
        
        if self.param_list['propensity_function']=='basic':
            self.calculate_propensity=self.calculate_propensity_basic
            self.dN=self.construct_dN() 
        elif self.param_list['propensity_function']=='Notch_dependent_degradation':
            self.calculate_propensity=self.calculate_propensity_Notch_dependent_degradation
            self.dN=self.construct_dN() 
        else:
            print("No propensity function assigned!!!")

    def read_param_file(self,param_file):
        param_list={}
        with open(param_file) as f:
            lines = f.readlines()
            for l in lines:
                try:
                    (p,val)=l.split(":")
                    try:
                        # save value as an integer
                        param_list[p]=np.float(val)
                    except ValueError:
                        # otherwise save as text with \n and \t stripped
                        param_list[p]=val.strip()
                except ValueError:
                    pass
        return (param_list)
        
    def S_ind(self,S_lbl,cell_lbl):
        return ( cell_lbl*self.N_S + self.species_ind[S_lbl])
 
    def construct_dN(self):
        # construct the stoichiometry matrix dN for the two cells
        dN=np.zeros((2*self.N_R,2*self.N_S))
        
        # for each cell i
        for i in [0,1]:
            # r0: Oc --> O
            dN[ i*self.N_R+0, self.S_ind('O',i) ] = +1
            dN[ i*self.N_R+0, self.S_ind('Oc',i) ] = -1
            # r1: O --> Oc
            dN[ i*self.N_R+1, self.S_ind('O',i) ] = -1
            dN[ i*self.N_R+1, self.S_ind('Oc',i) ] = +1
            # r2: OcH --> OH
            dN[ i*self.N_R+2, self.S_ind('OH',i) ] = +1
            dN[ i*self.N_R+2, self.S_ind('OcH',i) ] = -1
            # r3: OH --> OcH
            dN[ i*self.N_R+3, self.S_ind('OH',i) ] = -1
            dN[ i*self.N_R+3, self.S_ind('OcH',i) ] = +1
            # r4: O --> OH
            dN[ i*self.N_R+4, self.S_ind('OH',i) ] = +1
            dN[ i*self.N_R+4, self.S_ind('O',i) ] = -1
            # r5: OH --> O
            dN[ i*self.N_R+5, self.S_ind('OH',i) ] = -1
            dN[ i*self.N_R+5, self.S_ind('O',i) ] = +1
            # r6: Oc --> OcH
            dN[ i*self.N_R+6, self.S_ind('OcH',i) ] = +1
            dN[ i*self.N_R+6, self.S_ind('Oc',i) ] = -1
            # r7: OcH --> Oc
            dN[ i*self.N_R+7, self.S_ind('OcH',i) ] = -1
            dN[ i*self.N_R+7, self.S_ind('Oc',i) ] = +1
            # r8: OH --> OH + M
            dN[ i*self.N_R+8, self.S_ind('M',i) ] = +1
            # r9: M --> 0
            dN[ i*self.N_R+9, self.S_ind('M',i) ] = -1
            # r10: M --> D + M
            dN[ i*self.N_R+10, self.S_ind('D',i) ] = +1
            # r11: D --> 0
            dN[ i*self.N_R+11, self.S_ind('D',i) ] = -1
            # r12: OY --> OYH
            dN[ i*self.N_R+12, self.S_ind('OYH',i) ] = +1
            dN[ i*self.N_R+12, self.S_ind('OY',i) ] = -1
            # r13: OYH --> OY
            dN[ i*self.N_R+13, self.S_ind('OYH',i) ] = -1
            dN[ i*self.N_R+13, self.S_ind('OY',i) ] = +1
            # r14: OYH --> OYH + MY
            dN[ i*self.N_R+14, self.S_ind('MY',i) ] = +1
            # r15: MY --> 0
            dN[ i*self.N_R+15, self.S_ind('MY',i) ] = -1
            # r16: MY --> DY + MY
            dN[ i*self.N_R+16, self.S_ind('DY',i) ] = +1
            # r17: DY --> 0
            dN[ i*self.N_R+17, self.S_ind('DY',i) ] = -1
            
        return(dN)
            
    def initialize_reactants(self,init_data):
        N=np.zeros(2*self.N_S,dtype=int)
        # for each cell
        for i in [0,1]:
            # get initial conditions for each species
            init_list=init_data[i]
            # for each species label in init_list
            for lbl in init_list.keys():
                # add reactant particle number to array of reactant numbers
                N[ i*self.N_S + self.S_ind(lbl,0) ] = init_list[lbl]
        return(N)
        
    def calculate_propensity_basic(self,N,cell_type):
        a=np.zeros(2*self.N_R,dtype=float)
        
        # i represents the Cell label
        for i in [0,1]:
            # r0: Oc --> O
            a[ i*self.N_R+0 ] = self.param_list['k0_'+cell_type[i] ] * N[ self.S_ind('Oc',i) ]
            # r1: O --> Oc
            a[ i*self.N_R+1 ] = self.param_list['k1_'+cell_type[i] ] * N[ self.S_ind('O',i) ]
            # r2: OcH --> OH
            a[ i*self.N_R+2 ] = self.param_list['k0_'+cell_type[i] ] * N[ self.S_ind('OcH',i) ]
            # r3: OH --> OcH
            a[ i*self.N_R+3 ] = self.param_list['k1_'+cell_type[i] ] * N[ self.S_ind('OH',i) ]
            
            if i==0:
                S=self.param_list['phi_'+cell_type[i] ] * N[ self.S_ind('D',1) ]
            else:
                S=self.param_list['phi_'+cell_type[i] ] * N[ self.S_ind('D',0) ]
                
            K=self.param_list['K']
            n=self.param_list['n']
            H=self.param_list['KH']*(pow(K,n)+pow(S,n))/(pow(K,n)+(1+self.param_list['beta'])*pow(S,n))
            
            # r4: O --> OH
            a[ i*self.N_R+4 ] = self.param_list['fOH'] * H * N[ self.S_ind('O',i) ]
            # r5: OH --> O
            a[ i*self.N_R+5 ] = self.param_list['bOH'] * N[ self.S_ind('OH',i) ]
            # r6: Oc --> OcH
            a[ i*self.N_R+6 ] = self.param_list['fOH'] * H * N[ self.S_ind('Oc',i) ]
            # r7: OcH --> Oc
            a[ i*self.N_R+7 ] = self.param_list['bOH'] * N[ self.S_ind('OcH',i) ]
            # r8: OH --> OH + M
            a[ i*self.N_R+8 ] = self.param_list['v0_'+cell_type[i] ] * N[ self.S_ind('OH',i) ]
            # r9: M --> 0
            a[ i*self.N_R+9 ] = self.param_list['d0'] * N[ self.S_ind('M',i) ]
            # r10: M --> D + M
            a[ i*self.N_R+10 ] = self.param_list['fD_'+cell_type[i] ] * N[ self.S_ind('M',i) ]
            # r11: D --> 0
            a[ i*self.N_R+11 ] = self.param_list['bD_'+cell_type[i] ] * N[ self.S_ind('D',i) ]
            # r12: OY --> OYH
            a[ i*self.N_R+12 ] = self.param_list['fOYH'] * H * N[ self.S_ind('OY',i) ]
            # r13: OYH --> OY
            a[ i*self.N_R+13 ] = self.param_list['bOYH'] * N[ self.S_ind('OYH',i) ]
            # r14: OYH --> OYH + MY
            a[ i*self.N_R+14 ] = self.param_list['fMY_'+cell_type[i] ] * N[ self.S_ind('OYH',i) ]
            # r15: MY --> 0
            a[ i*self.N_R+15 ] = self.param_list['bMY_'+cell_type[i] ] * N[ self.S_ind('MY',i) ]
            # r16: MY --> DY + MY
            a[ i*self.N_R+16 ] = self.param_list['fDY_'+cell_type[i] ] * N[ self.S_ind('MY',i) ]
            # r17: DY --> 0
            a[ i*self.N_R+17 ] = self.param_list['bDY_'+cell_type[i] ] * N[ self.S_ind('DY',i) ]
        
        return(a)

    def calculate_propensity_Notch_dependent_degradation(self,N,cell_type):
        a=np.zeros(2*self.N_R,dtype=float)
        
        for i in [0,1]:
            # r0: Oc --> O
            a[ i*self.N_R+0 ] = self.param_list['k0_'+cell_type[i] ] * N[ self.S_ind('Oc',i) ]
            # r1: O --> Oc
            a[ i*self.N_R+1 ] = self.param_list['k1_'+cell_type[i] ] * N[ self.S_ind('O',i) ]
            # r2: OcH --> OH
            a[ i*self.N_R+2 ] = self.param_list['k0_'+cell_type[i] ] * N[ self.S_ind('OcH',i) ]
            # r3: OH --> OcH
            a[ i*self.N_R+3 ] = self.param_list['k1_'+cell_type[i] ] * N[ self.S_ind('OH',i) ]
            
            if i==0:
                S=self.param_list['phi_'+cell_type[i] ] * N[ self.S_ind('D',1) ]
            else:
                S=self.param_list['phi_'+cell_type[i] ] * N[ self.S_ind('D',0) ]
                
            K=self.param_list['K']
            n=self.param_list['n']
            H=self.param_list['KH']*(pow(K,n)+pow(S,n))/(pow(K,n)+(1+self.param_list['beta'])*pow(S,n))
            
            # r4: O --> OH
            a[ i*self.N_R+4 ] = self.param_list['fOH'] * H * N[ self.S_ind('O',i) ]
            # r5: OH --> O
            a[ i*self.N_R+5 ] = self.param_list['bOH'] * N[ self.S_ind('OH',i) ]
            # r6: Oc --> OcH
            a[ i*self.N_R+6 ] = self.param_list['fOH'] * H * N[ self.S_ind('Oc',i) ]
            # r7: OcH --> Oc
            a[ i*self.N_R+7 ] = self.param_list['bOH'] * N[ self.S_ind('OcH',i) ]
            # r8: OH --> OH + M
            a[ i*self.N_R+8 ] = self.param_list['v0_'+cell_type[i] ] * N[ self.S_ind('OH',i) ]
            # r9: M --> 0
            a[ i*self.N_R+9 ] = self.param_list['d0'] * N[ self.S_ind('M',i) ]
            # r10: M --> D + M
            a[ i*self.N_R+10 ] = self.param_list['fD'] * N[ self.S_ind('M',i) ]
            # r11: D --> 0
            K1=self.param_list['K1']
            n=1
            bD=(pow(K1,n)*self.param_list['bD_min']+pow(S,n)*self.param_list['bD_max'])/(pow(K1,n)+pow(S,n))
            a[ i*self.N_R+11 ] =  bD * N[ self.S_ind('D',i) ]
            # r12: OY --> OYH
            a[ i*self.N_R+12 ] = self.param_list['fOYH'] * H * N[ self.S_ind('OY',i) ]
            # r13: OYH --> OY
            a[ i*self.N_R+13 ] = self.param_list['bOYH'] * N[ self.S_ind('OYH',i) ]
            # r14: OYH --> OYH + MY
            a[ i*self.N_R+14 ] = self.param_list['fMY_'+cell_type[i] ] * N[ self.S_ind('OYH',i) ]
            # r15: MY --> 0
            a[ i*self.N_R+15 ] = self.param_list['bMY_'+cell_type[i] ] * N[ self.S_ind('MY',i) ]
            # r16: MY --> DY + MY
            a[ i*self.N_R+16 ] = self.param_list['fDY_'+cell_type[i] ] * N[ self.S_ind('MY',i) ]
            # r17: DY --> 0
            a[ i*self.N_R+17 ] = self.param_list['bDY_'+cell_type[i] ] * N[ self.S_ind('DY',i) ]
        
        return(a)

    def print_reactions(self,a,i):
        print("R[%d]:"%i,end=' ')
        r=np.where(self.dN[i,:]==-1)[0]
        if len(r)>0:
            if r[0]>=self.N_S:
                ind=r[0]-self.N_S
                print(self.species_lbl[ind]+'_p',end='')
            else:
                ind=r[0]
                print(self.species_lbl[ind]+'_a',end='')
        else:
            print('0',end='')
        print("-->",end=''),
        r=np.where(self.dN[i,:]==1)[0]
        if len(r)>0:
            if r[0]>=self.N_S:
                ind=r[0]-self.N_S
                print(self.species_lbl[ind]+'_p',end='')
            else:
                ind=r[0]
                print(self.species_lbl[ind]+'_a',end='')
        else:
            print('0',end='')
        print(", a=%f"%a[i])

    def print_molecules(self,t,N):
        cell_lbl=['_a','_p']
        print("t=%f"%t,end=' ')
        for i in [0,1]:
            for j in range(0,self.N_S):
                print(self.species_lbl[j]+cell_lbl[i]+":%d"%N[i*self.N_S+j],end=' ')
        print()

    def init_save_data(self,t,N,species_list):
        sim_data=[]
        # first, save time
        sim_data.append(np.array([t]))
        for S_lbl in species_list:
            if S_lbl=='open_complex':
                N_save = N[ [self.S_ind('O',0),self.S_ind('O',1)] ] + N[ [self.S_ind('OH',0),self.S_ind('OH',1)] ]
            else:
                N_save = N[ [self.S_ind(S_lbl,0),self.S_ind(S_lbl,1)] ]
            sim_data.append( np.array([N_save]) )
        return(sim_data)

    def add_save_data(self,sim_data,t,N,species_list):
        # first, save time
        sim_data[0]=np.hstack((sim_data[0],t))
        c=1
        for S_lbl in species_list:
            if S_lbl=='open_complex':
                N_save = N[ [self.S_ind('O',0),self.S_ind('O',1)] ] + N[ [self.S_ind('OH',0),self.S_ind('OH',1)] ]
            else:
                N_save = N[ [self.S_ind(S_lbl,0),self.S_ind(S_lbl,1)] ]
            sim_data[c]=np.vstack((sim_data[c],N_save))
            c=c+1
        return(sim_data)

    def propagate_reactions(self,t,N,cell_types,T_sim,DT_out,save_species_list,generate_moments=False):

        # calculate the time when the next data is saved
        t_out_new = np.ceil(t/DT_out)*DT_out
        
        # initialize stored data
        trajectory_data=self.init_save_data(t,N,save_species_list)

        # moments
        if generate_moments:
            T=0
            N_m=np.zeros((2,2*self.N_S),dtype=float)
        
        # while time is smaller than t_sim    
        while t<T_sim:
            # calculate propensity function
            a=self.calculate_propensity(N,cell_types)
            # as well as cumulative propensity        
            A=a.cumsum()
            # and total propensity
            A_tot=A[-1]
            
            # draw time to next reaction
            dt=np.random.exponential(1/A_tot,1)[0]
    
            # draw reaction
            # first, draw random number r        
            r=A_tot*np.random.rand()
            # find which reaction occurs according to cumulative propensity
            cont=True
            i=0
            while cont:
                if r<=A[i]:
                    react_ind=i
                    cont=False
                i=i+1

            if generate_moments:   
                T += dt
                N_m[0,:] += dt*N
                N_m[1,:] += dt*N**2     
            # add time        
            t=t+dt
            
            # save data
            if t>t_out_new:
                trajectory_data=self.add_save_data(trajectory_data,t_out_new,N,save_species_list)
                t_out_new += DT_out
    
            # change number of reactants according to reaction
            N=N+self.dN[react_ind,:]

        save_data={'trajectory_data':trajectory_data}

        if generate_moments:               
            N_m=N_m/T
            save_data['moment_data']=N_m
                        
        return(t,N,save_data)

    def divide_cell(self,N,div_cell,div_species_lbl_list):
        # for each reactant species in cell <div_cell>
        for i in range(0,self.N_S):
            # check is species in <div_species_lbl_list>
            if (self.species_lbl[i] in div_species_lbl_list):
                # if so, segregate molecules into daughter cell in a binomial manner
                ind = self.S_ind( self.species_lbl[i], div_cell)
                N[ind] = np.random.binomial( N[ind], 0.5 )
            # if not, leave like it is. This means that daughter cell will inherit the
            # same level or state as the mother cell
        return(N)

    def run_full_simulation(self,N,save_species_list):

        t=0    
        run_data=[]
        
        # 1) simulate two mother cells
        (t,N,sim_data) = self.propagate_reactions(t,N,['m','m'],self.param_list['T0'],self.param_list['DT_out'],save_species_list)
        run_data.append(sim_data)
        # 2) execute first division
        N=self.divide_cell(N,0,['M', 'D','MY','DY'])
        # 3) simulate daughter + mother cell
        (t,N,sim_data) = self.propagate_reactions(t,N,['d','m'],self.param_list['T0']+self.param_list['DT'],self.param_list['DT_out'],save_species_list)
        run_data.append(sim_data)
        # 4) execute second division
        N=self.divide_cell(N,1,['M', 'D','MY','DY'])
        # 5) simulate wo daughter cells 
        (t,N,sim_data) = self.propagate_reactions(t,N,['d','d'],self.param_list['T0']+self.param_list['T1'],self.param_list['DT_out'],save_species_list)
        run_data.append(sim_data)
     
        return(run_data)
        
    #### sweep functions        
        
    def get_sweep_file_name(self,sweep_data_dir,sweep_data_file,sweep_param_lbl):
        cwd = os.getcwd()
        return os.path.join(cwd, sweep_data_dir, sweep_data_file + '_' + sweep_param_lbl + '.p')
        
    def init_sweep(self,sweep_data_dir,sweep_data_file,sweep_param_lbl):
    
        data_file=self.get_sweep_file_name(sweep_data_dir,sweep_data_file,sweep_param_lbl)
        
        try:
            # try loading file
            tmp = pickle.load( open( data_file, "rb" ) )
            print("Loading existing sweep data output file <%s>..."%data_file)
            sweep_data = tmp
            
        except FileNotFoundError:
            # if not found, first try to create sweep_data_dir
            print("Starting new sweep data output file <%s>..."%data_file)
            try:
                os.makedirs(sweep_data_dir)
                sweep_data=[]
                print("\tCreating new sweep data directory <%s>..."%sweep_data_dir)
            except FileExistsError:
                sweep_data=[]
        print('\n')
        return (sweep_data)
    
    def parallel_sweep( self, zipped_data ):
        seed = zipped_data[0]
        save_species_list = zipped_data[1]
        sweep_param_lbl = zipped_data[2]
        sweep_param = zipped_data[3]

        print("Running simulation %d for %s=%f"%(seed,sweep_param_lbl,sweep_param))

        # set random seed
        np.random.seed(seed)        

        # initialize reactant levels N, using initial conditions <IC> for each cell
        N=self.initialize_reactants([{'O':1,'OY':1},{'O':1,'OY':1}])
        # and calculate stoichiometry matrix
#        dN=self.construct_dN()

        sim_data=self.run_full_simulation(N,save_species_list)

        data_entry = {'random_seed':seed,'sim_data':sim_data}
        return data_entry
 

    def do_sweep(self,sweep_param_lbl,sweep_param_list,N_run,sweep_data_dir,sweep_data_file,save_species_list,seedStart=0,coreFraction=1.):
    
        sweep_data=self.init_sweep(sweep_data_dir,sweep_data_file,sweep_param_lbl)
        
        for sweep_param in sweep_param_list:
            tStart = time.time()
            # set sweep parameter to current value
            self.param_list[sweep_param_lbl] = sweep_param
        
            sweep_data_ind=-1
            for i in range(0,len(sweep_data)):
                if sweep_data[i]['param_value']==sweep_param:
                    sweep_data_ind=i
            if sweep_data_ind != -1:
                run_data= sweep_data[ sweep_data_ind ]['run_data']
            else:
                run_data=[]
                
            # print( [ (i['param_value'], len(i['run_data'])) for i in sweep_data ] )
            num_new_sim = int( N_run - len(run_data) )
#            new_run_data = [ None for i in range(num_new_sim) ]
            zipped_data = zip( np.arange(num_new_sim) + seedStart, 
                                [save_species_list for i in range(num_new_sim)],
                                [sweep_param_lbl for i in range(num_new_sim)],
                                [sweep_param for i in range(num_new_sim)]
                                )

            if os.cpu_count()>1:
                print('Multiple cores found. parallelize over ', int( os.cpu_count() * 1. ),' cores.')
                with Pool( int( os.cpu_count() * coreFraction ) ) as pool:
                    new_run_data = pool.map( self.parallel_sweep, zipped_data )

            for data_entry in new_run_data:
                run_data.append(data_entry)

            if sweep_data_ind != -1:
                # if there were already sim with same sweep_param, redefine only run_data
                sweep_data[ sweep_data_ind ]['run_data'] = run_data
            else:
                # if it's a new sweep_param, append a new dict to sweep_data
                sweep_data.append({'param_value':sweep_param,'run_data':run_data})

            pickle.dump( sweep_data, open( self.get_sweep_file_name(sweep_data_dir,sweep_data_file,sweep_param_lbl), "wb" ) )
            print('Time elapsed to run ' + str(num_new_sim) + ' simulations: ' + str(time.time()-tStart) + '\n' )
