# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 09:57:19 2016

@author: Adrien
"""


from canpy.useful import time_to_snap

import numpy as np
import matplotlib.pyplot as plt
import math





def mixing(dl,time,dist_thresh='standard'):
    '''
    The mixing function computes the atomic mixing at the time provided by the
    user. The threshold provided is a distance in Angström.
    dump_array has to be a numpy array returned by the array method of the
    dumplog class. The mixing function return a list which contains the results.
    The first value of the list is the calculation of the value which appraise
    the mixing and the second element of the list provide the number of atoms
    which have been taken accounted in the calculation.
    The time has to be given in ps to get consistent results.
    If the distance between the original position of an atom and its position
    at the time provided by the user is smaller distance than 
    threshold*lattice_parameter then the contribution to the mixing of this
    atom is ignored. By default the threshold value is set at 0% of the 
    lattice parameter but the user is free to change it to obtain more
    significant results. The main idea behind the use of such a threshold is
    to avoid the contribution of the vibrational movements of the atoms to
    the mixing value.
    '''
    
    # Determination of the threshold
    if dist_thresh == 'standard':
        if dl.lattice == 'bcc':
            dist_thresh = (math.sqrt(3)/4) * dl.lat_param
        
        elif dl.lattice == 'fcc':
            dist_thresh = (math.sqrt(2)/4) * dl.lat_param
        else:
            print('Error in determination of the threshold of the mixing function')
    
    
    
    # Here we define a variable giving the coresponding snapshot for the time provided
    n_snap=round((time/dl.time_step)/dl.inter)
    
    # Sort of the rows of the numpy array corresponding to the 0th snapshot using a quicksort algorithm according to increasing atom id order
    to_sort_0=np.copy(dl.data[0])
    sorted_array_0 = to_sort_0[to_sort_0[:,0].argsort()]
    
    # Sort of the rows of the numpy array corresponding to the n_snapth snapshot using a quicksort algorithm according to increasing atom id order
    to_sort_n=np.copy(dl.data[n_snap])   
    sorted_array_n = to_sort_n[to_sort_n[:,0].argsort()]
    
    # Attribution of the dimensions of the box using
    length=dl.box_dim[0]
    
    # atom_count is used to count the number of atoms that we consider involved in the mixing
    atom_count=0
    
    # tot is used to sum all the contributions to the mixing. It is actually a sum of norms
    tot=0
    
    # Initialisation of the counter for the first while loop
    i=0
    
    while i<dl.nb_atoms:
        
        # Use of variable with more explicit meaning for the clearness of the code
        x0=sorted_array_0[i][2]
        y0=sorted_array_0[i][3]
        z0=sorted_array_0[i][4]
        x1=sorted_array_n[i][2]
        y1=sorted_array_n[i][3]
        z1=sorted_array_n[i][4]
        
        # Coordinates of the vector steming from the soustraction of the final position vector with the initial position vector
        dx=x1-x0
        dy=y1-y0
        dz=z1-z0
        
        if dx > length/2:
            x1=x1-length
        elif dx < -length/2:
            x1=x1+length
    
        if dy > length/2:
            y1=y1-length
        elif dy < -length/2:
            y1=y1+length
    
        if dz > length/2:
            z1=z1-length
        elif dz < -length/2:
            z1=z1+length
    
        dx=x1-x0
        dy=y1-y0
        dz=z1-z0
    
        norm=math.sqrt(dx**2 + dy**2 + dz**2)
        
        # Not finshed here to the extend that the threshold is supposed to be a percentage of the lattice parameter
        if norm >= dist_thresh:
            tot+=norm
            atom_count+=1
        
        i+=1
    
    
    # Declaration of a list which contains the results
    if atom_count != 0:
        result=[tot/atom_count,atom_count,n_snap]
    else:
        result=[0,atom_count,n_snap]
    
    return result










def plot_mixing(dl,dist_thresh=0,t0=0,tf='standard',name='plot_mixing.png',style='bo'):
    '''
    The idea of this function is to build a plot of the mixing in function of
    the time in ps.
    
    '''
    
    if tf == 'standard':
        tf = dl.time_real
    
    
    # Convert t0 and tf in terms of number of snapshot
    snap_0=time_to_snap(t0,dl)
    snap_f=time_to_snap(tf,dl)
    
    # Initialisation of the lists which will contain the x and y axis
    list_mix=[]
    list_time=[]
    
    if t0==0:
        i=snap_0+1
        list_mix+=[0]
        list_time+=[0]
    else:
        i=snap_0
    
    while i <= snap_f:
        
        curr_time=i*dl.inter*dl.time_step
        transi=mixing(dl,curr_time,dist_thresh=dist_thresh)
        
        list_mix+=[transi[0]]
        list_time+=[curr_time]
        
        i+=1
    

    # Plot the mixing according to the time
    fig=plt.figure()
    plt.subplot(111)
    plt.plot(list_time,list_mix,style)
    plt.xlabel('Time [ps]')
    plt.ylabel('Mixing [Angström]')
    plt.title('Mixing from '+str(t0)+' ps to '+str(tf)+' ps')
    plt.savefig('plot_mixing_'+str(t0)+'ps_to_'+str(tf)+'ps.png')
    
    plt.close(fig)
    
    
    