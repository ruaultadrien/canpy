# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 17:46:57 2016

@author: Adrien
"""



import scipy.constants
import math




    
def index_heading(dl,name):
    '''
    This function is made to provide the index of the element corresponding to
    the name (string) provided in the list called heading which is a variable 
    output of the dumplog class.
    Therefore the user has to provide a dumplog object and the name of one of
    the columns of the dump file and the function returns the index for which
    this name (string) is stored in the heading list.
    '''
    
    # Definition of the index which starts at 0
    i=0
    
    while dl.heading[i] != name:
        i+=1
    # We have to ignore the two first elements : 'ITEM:' and 'ATOMS'
    i=i-2
    return i












def time_to_snap(time,dl):
    '''
    This function is used to convert a given time in its corresponding snapshot.
    It returns the index of the snapshot related to the time provided.
    '''
    
    return round((time/dl.time_step)/dl.inter)
    
    
    







def snap_to_time(snap,dl):
    '''
    This function is used to convert a snap index in its corresponding time in
    ps.
    '''
    
    return round((snap*dl.time_step*dl.inter),8)







def PKA_vector(u,v,w,pka_energy,amu):
    '''
    This function gives the velocity vector to apply to the PKA at the very beginning
    of the simulation cascade according to different parameters. The returned
    velocity vector is given in Angstrom/ps.
    The parameters that must be provided are the direction of the vector through
    the (u, v, w) variables, the PKA energy in keV through the variable pka_energy and
    finally the atomic mass unit of the atom in kg through amu.
    '''
    
    # Conversion of the PKA energy from eV to joules
    eV_to_J = scipy.constants.physical_constants['electron volt-joule relationship'][0]
    
    energy = pka_energy * 1e3 * eV_to_J
    
    mass = amu * scipy.constants.physical_constants['atomic mass constant'][0]
    
    a = math.sqrt((2*energy)/(mass*(u**2 + v**2 + w**2)))
    
    #Conversion of a from m/s to Angstrom/ps
    a = a*1e-2

    # Construction of the initial velocity vector
    velocity = [a*u,a*v,a*w]
    
    return velocity
    







