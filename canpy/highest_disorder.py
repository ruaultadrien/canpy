# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 09:55:44 2016

@author: Adrien
"""



import scipy.constants
from canpy.useful import index_heading 
# necessary to import starting with 'package.' because the current directory is the one of the main.py module

# Check if it's compulsory to import this function as importing dumplog doesn't seem to be compulsory
#from index_heading import index_heading





def highest_disorder(dl,peng_name,threshold='npt_standard'):
    '''
    The purpose of this function is to return the snapshot for which the 
    disorder is at its maximum.
    The notion of disorder is not obvious and not well define. That is why we
    need to provide a clear definition of how the disorder is measured through
    this function. Knowing that in a perfect crystal all the atoms should be
    at a minimum level of energy (the cohesive energy), more important is the
    disorder in a crystal higher should be the potential energy of the atoms.
    Thus the idea here is to count the number of atoms with a potential energy
    above a given threshold (standard or provided by the user). Indeed we can
    reasonably consider that if the threshold is well chosen, this atom count
    can be a good indicator of the disorder. Therefore basically the function
    count the number of atoms in a high disorder state for each snapshot and
    return the one with the highest number of atoms counted.
    The threshold can be chosen as standard if an equilibration simulation in
    npt environment has been run before the cascade. In this case the temperature
    will be directly extracted through the dl variable which has to be a dumplog
    class. If this condition is not respected, the user is invited to provide
    his own threshold.
    The treshold in the standard condition is calculated using the following
    formula : threshold=10*((3/2)*k*T) with k, the Boltzmann constant and T the 
    temperature extracted from dl (dumplog class). It is thus advised to the
    user to estimate the temperature and to use a similar kind of threshold 
    to get consistent results. 
    The cohesive energy is estimated using the average of the potential energy
    of the first snapshot (before that the cascade starts really) even if it
    actually estimates the cohesive energy plus the average kinetic energy.
    dl is a dumplog variable (see the first class definition).
    threshold must be provided in eV if it is not let as 'npt_standard'.
    Finally note that the user must provide the name which is used to
    characterise the column in the dump file where the potential energy of the
    atoms is stored.
    '''
    
    
    # Attribution of the Boltzmann constant
    k=scipy.constants.k
    
    T=dl.temperature
    
    # Joule -> eV conversion
    J_to_eV=scipy.constants.physical_constants['joule-electron volt relationship'][0]
    
    
    if threshold == 'npt_standard':
        # Be careful to the temperature: extract it in dumplog class
        threshold = 5*((3/2)*k*T)
        # Conversion in eV
        threshold = J_to_eV * threshold
    
    # We need an average of the (cohesive energy + average kinetic energy) using the first snapshot
    # First we extract the index of the potential energy information
    index_peng=index_heading(dl,peng_name)
    
    # Now we compute the average of the (cohesive energy + average kinetic energy)
    tot=0
    j=0
    while j<dl.nb_atoms:
        epot=dl.data[0][j][index_peng]
        tot+=epot
        j+=1
    
    average=tot/dl.nb_atoms
    
    
    # Final is the number max of atoms with potential energy greater than average + threshold
    final=0    
    
    k=0
    while k < dl.nb_snapshot:
        l=0
        
        # Declaration of the counter to count the number of atoms with potential energy greater than the average + threshold
        count=0
        
        while l < dl.nb_atoms:
            epot=dl.data[k][l][index_peng]
            if epot > average + threshold:
                count+=1
            l+=1
        
        if count > final:
            final=count
            snap=k
        
        k+=1
                
    return [snap,final]
    
    
    
    
    
    
    
    '''
    
    
    
    
def index_heading(dl,name):
'''
    '''
    This function is made to provide the index of the element corresponding to
    the name (string) provided in the list called heading which is a variable 
    output of the dumplog class.
    Therefore the user has to provide a dumplog object and the name of one of
    the columns of the dump file and the function returns the index for which
    this name (string) is stored in the heading list.
    '''
    '''
    # Definition of the index which starts at 0
    i=0
    
    while dl.heading[i] != name:
        i+=1
    # We have to ignore the two first elements : 'ITEM:' and 'ATOMS'
    i=i-2
    return i
    
    '''
    
    