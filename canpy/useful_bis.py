# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 17:14:42 2016

@author: Adrien
"""

from canpy.dumplog import dumplog
from canpy.point_defects import point_defects
from canpy.mixing import mixing,plot_mixing
from canpy.useful import time_to_snap
import matplotlib.pyplot as plt

import time



def complete_analysis(dump_name,log_name,equi=True):
    '''
    The idea of this function is to implement a basic analyse routin on the results
    given by LAMMPS. 
    The idea is basically to print graphs of the thermal spike and of the final
    configuration.
    The function also aims to print the residual defects (at the last snapshot)
    as well as the interstitial and vacancy clustering final configuration.
    '''
    
    t=time.clock()
    
    
    # Declaration of a dumplog object
    a=dumplog(dump_name,log_name,equi=equi)
    
    # Save a movie of the evolution of the defects
    a.defects_movie()
    
    # Save a plot of the defects during the thermal spike
    time_spike=a.highest_disorder()['time']
    spike=point_defects(a,time_spike)
    spike.plot_centered(name='plot_highest_disorder')
    
    
    # Save a plot of the defects at the enfd of the simulation
    final=point_defects(a,a.time_real)
    final.plot_centered(name='plot_final_defects_configuration')
    
    # Videos of the final defects configuration with rotation around the z axis
    final.rot_frame_anim()
    final.rot_fixframe_anim()
    
    
    # Write some details about the simulation
    file=open('Results.txt','w')
    file.write('Number of atoms: '+str(a.nb_atoms))
    file.write('\nType of lattice: '+a.lattice)
    file.write('\nLattice parameter: '+str(a.lattice)+' A')
    file.write('\nLength of the box: '+str(a.box_dim[0])+' A')
    file.write('\nTime step: '+str(a.time_step)+' ps')
    file.write('\nTotal number of time steps run: '+str(a.step_max))
    file.write('\nTotal time of the simulation: '+str(a.time_real)+' ps')
    
    
    
    # Write the results (defects and clusters) in a text file
    file.write('\n\nNumber of Frenkel Pairs: '+str(final.defects['FP_count']))
    
    NN1_clust=final.clusters(NN=1)
    NN2_clust=final.clusters(NN=2)
    NN3_clust=final.clusters(NN=3)
    NN4_clust=final.clusters(NN=4)
    file.write('\nInterstitial clustering configuration: ')
    file.write('\n\t1st nearest neighbour criterium: '+str(NN1_clust['inter']['shape']))
    file.write('\n\t2nd nearest neighbour criterium: '+str(NN2_clust['inter']['shape']))
    file.write('\n\t3rd nearest neighbour criterium: '+str(NN3_clust['inter']['shape']))
    file.write('\n\t4th nearest neighbour criterium: '+str(NN4_clust['inter']['shape']))
    file.write('\nVacancy clustering configuration: ')
    file.write('\n\t1st nearest neighbour criterium: '+str(NN1_clust['vac']['shape']))
    file.write('\n\t2nd nearest neighbour criterium: '+str(NN2_clust['vac']['shape']))
    file.write('\n\t3rd nearest neighbour criterium: '+str(NN3_clust['vac']['shape']))
    file.write('\n\t4th nearest neighbour criterium: '+str(NN4_clust['vac']['shape']))
    
    # PLot a graph of the mixing
    plot_mixing(a,a.time_real)
    
    # Measure of the execution time
    t=time.clock()-t
    file.write('\nExecution time of basic_analysis: '+str(t)+' s')
    
    file.close()
    
    
    
    
    
    
    
    
    
    
    
def graph_defects(dl, t_start=0, t_end='standard'):
    '''
    This function allows to plot a graph of the defects as a function of time.
    The user has to provide a dumplog object. He can also select the departure
    time as well as the ending time in ps through the variables t_start and
    t_end. If the user doesn't change the value of t_start and t_end they will
    automatically be set such that the graph is made on the whole range of time
    of the simulation.
    '''
    
    if t_end == 'standard':
        t_end = dl.time_real
    
    interval = dl.inter * dl.time_step
    
    # Definition of the arrays used to plot the graph
    defects_array = []
    time_array = []
    max_FP = 0
    
    current_time = t_start
    while abs(current_time - t_end) > 1e-6:
        alpha = point_defects(dl,current_time)
        defects_array += [alpha.defects['FP_count']]
        print(alpha.defects['FP_count'])
        
        if alpha.defects['FP_count'] > max_FP:
            max_FP = alpha.defects['FP_count']
        
        time_array += [current_time]
        
        current_time += interval
    
    
    #Ploting of defects as a function of time
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(time_array,defects_array)
    ax.set_xlim(0,t_end*1.01)
    ax.set_ylim(0,1.1*max_FP)
    ax.set_xlabel('Time [ps]')
    ax.set_ylabel('Number of Frenkel pairs')
    ax.set_title('Number of Frenkel pairs from '+str(t_start)+' ps to '+str(t_end)+' ps')
    fig.savefig('plot_FPcount_'+str(t_start)+'ps_to_'+str(t_end)+'ps.png')
    
    plt.close(fig)

    
    
    
    
    
    
    
def graph_semilog_defects(dl, t_start=0, t_end='standard'):
    '''
    This function allows to plot a graph of the defects as a function of time.
    The user has to provide a dumplog object. He can also select the departure
    time as well as the ending time in ps through the variables t_start and
    t_end. If the user doesn't change the value of t_start and t_end they will
    automatically be set such that the graph is made on the whole range of time
    of the simulation.
    '''
    
    if t_end == 'standard':
        t_end = dl.time_real
    
    interval = dl.inter * dl.time_step
    
    # Definition of the arrays used to plot the graph
    defects_array = []
    time_array = []
    max_FP = 0
    
    current_time = t_start
    while abs(current_time - t_end) > 1e-6:
        alpha = point_defects(dl,current_time)
        defects_array += [alpha.defects['FP_count']]
        print(alpha.defects['FP_count'])
        
        if alpha.defects['FP_count'] > max_FP:
            max_FP = alpha.defects['FP_count']
        
        time_array += [current_time]
        
        current_time += interval
    
    
    
    # Definition of the maximum value for the y-axis of the graph
    test = 1
    while max_FP/test > 1:
        test *= 10
    sup_y = test

    
    #Ploting of defects as a function of time
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time_array,defects_array)
    ax.set_yscale('log')
    ax.set_xlabel('Time [ps]')
    ax.set_ylabel('Number of Frenkel pairs')
    ax.set_title('Number of Frenkel pairs from '+str(t_start)+' ps to '+str(t_end)+' ps')
    ax.set_xlim(0,t_end*1.01)
    ax.set_ylim(1,sup_y)
    
    fig.savefig('plot_FPcount_'+str(t_start)+'ps_to_'+str(t_end)+'ps.png')
    
    plt.close(fig)

    
    
    








def graph_log_defects(dl, t_start=0, t_end='standard'):
    '''
    This function allows to plot a graph of the defects as a function of time.
    The user has to provide a dumplog object. He can also select the departure
    time as well as the ending time in ps through the variables t_start and
    t_end. If the user doesn't change the value of t_start and t_end they will
    automatically be set such that the graph is made on the whole range of time
    of the simulation.
    '''
    
    if t_end == 'standard':
        t_end = dl.time_real
    
    interval = dl.inter * dl.time_step
    
    # Definition of the arrays used to plot the graph
    defects_array = []
    time_array = []
    max_FP = 0
    
    current_time = t_start
    while abs(current_time - t_end) > 1e-6:
        alpha = point_defects(dl,current_time)
        defects_array += [alpha.defects['FP_count']]
        print(alpha.defects['FP_count'])
        
        if alpha.defects['FP_count'] > max_FP:
            max_FP = alpha.defects['FP_count']
        
        time_array += [current_time]
        
        current_time += interval
    
    
    
    # Definition of the maximum value for the y-axis of the graph
    test = 1
    while max_FP/test > 1:
        test *= 10
    sup_y = test
    
    test = 1
    while t_end/test > 1:
        test *= 10
    sup_x = test
    
    test = 1
    while interval/test < 1:
        test = test/10
    inf_x = test
    
    
    
    #Ploting of defects as a function of time
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(time_array,defects_array)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Time [ps]')
    ax.set_ylabel('Number of Frenkel pairs')
    ax.set_title('Number of Frenkel pairs from '+str(t_start)+' ps to '+str(t_end)+' ps')
    ax.set_xlim(inf_x,sup_x)
    ax.set_ylim(1,sup_y)
    
    fig.savefig('plot_FPcount_'+str(t_start)+'ps_to_'+str(t_end)+'ps.png')
    
    plt.close(fig)




    
    
    
    
    
    
    
    
    

def T_infile(nb_cells, temp):
    '''
    This function aims to create a text file to initialise the electronic 
    temperatures of the cells created using fix ttm in LAMMPS.
    The user has to provide the number of cells in which he wants to divide
    each side of the simulation cell. We advise to cut the simulation box in
    cells such that each of them contains approximately 50 to 500 atoms depending
    on the level of accuracy desired. This information must be provided using
    the variable nb_cells.
    The user also has to provide the temperature to initialise the cells with
    using the variable temp.
    The simulation box must be a cube
    '''
    
    # Opening of the file
    file = open('T'+str(temp)+'_Cell'+str(nb_cells)+'_infile.in','w')
    
    i=0
    
    
    while i < nb_cells:
        
        j=0
        while j < nb_cells:
            
            k=0
            while k < nb_cells:
                file.write(str(i)+' '+str(j)+' '+str(k)+' '+str(temp)+'\n')
                
                k+=1
            
            j+=1
        
        i+=1
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    