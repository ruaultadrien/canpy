# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 09:52:18 2016

@author: Adrien
"""



import numpy as np
import scipy.constants




class former_dumplog:
    '''
    To call this class in a programm code one has to provide: the name of the 
    dump files in the same way as it has been called in the LAMMPS's input file
    (with the * character) as well as the name of the log file generated.
    The idea of this class is to extract some important informations needed to 
    study the outcome of the simulation. This class prevent the user to input
    manually in the python programm a lot of data such as the results produced
    by LAMMPS but also some informations that have already been written in the
    input file which would be boring and time consuming to provide again.
    For the moment the use of this class has some requirements on the LAMMPS
    input: you have to use the units metal, you have to decompose your dump
    file in separeted files using the * character (this problem could be
    solved easily by inserting a default argument which could be used to tell
    that the output informations are in a unique dump file and the latter
    could be then decomposed in several textfile using python)
    '''
    
    def __init__(self,dump_name,log_name,equi=True):
        '''
        Initialisation of the different attributes of the dumplog class
        Basically here we extract the different important data provided
        through the input text file. However we do not use the text file. We
        only use the two output files of a LAMMPS simuation which are the 
        dump file and the log file.
        The variable called equi is used to tell if an equilibration stage
        has been run before that the before the real cascade stage to be run.
        '''
        # Extraction of dump informations
        a=dump_name.split('*')
        first_dump=a[0]+'0'+a[1]    # Give the name of the first dump file so that it can be read
        fdump=open(first_dump,'r')
        
        i=0
        while i<9:
            line=fdump.readline()
            if i==3:
                nb_atoms=int(line)
            
            elif i==5:
                # Extraction of the length of the box in Angström using lx
                transi=line.split()
                box_dim=[float(transi[1]) - float(transi[0])]
                
                # Add the beginning and the end of the lx segment
                box_dim+=[float(transi[0])]
                box_dim+=[float(transi[1])]                
                
            elif i==8:
                heading=line.split()
                
            i+=1
            
        fdump.close()
        
        #--------------------------------------------------
        # Extraction of log informations
        
        temperature='non fixed'
        
        # Test used in the following if statements to not take the number of time_step of the equilibration phase
        test_timestep=equi
        
        # In case equi=True, test used to take the temperature at the first fix statement and avoid the second one
        test_fix=equi        
        
        flog=open(log_name,'r')
        
        line='a'
        while line!='':
            line=flog.readline()
            if line[0:3]=='run':
                a=line.split()
                step_max=int(a[1])
                
            elif line[0:4]=='dump':
                a=line.split()
                inter=int(a[4])
            
            elif line[0:7]=='lattice':
                a=line.split()
                lattice=a[1]
            
            
            # Extraction of the temperature only if an equilibration phase has been run in an npt ensemble
            elif line[0:3]=='fix' and test_fix==True:
                a=line.split()
                temperature=int(a[6])
                test_fix=False
                
            
            # Be careful for the extraction of the timestep we need to state two elif statement in case of an equilibration simulation has been run
            elif line[0:8]=='timestep' and test_timestep==False:
                a=line.split()
                time_step=eval(a[1])
                
            elif line[0:8]=='timestep' and test_timestep==True:
                test_timestep=False
                
        flog.close()
        
        nb_snapshot=int(step_max//inter + 1)
        time_real=round((step_max*time_step),12)
        
        
        '''
        The idea here is to extract the numerical data contained in the dump
        file using a numpy array. We will use the loadtxt function contained
        in the numpy package to extract this data and store them in an numpy
        array.
        '''
        # This part is given to build a numpy array which contains the data present in the dump file
        cut_name=dump_name.split('*')
        
        
        # Initialisation of the array which will host the data contained in the dump file
        #data=numpy.empty(nb_snapshot,numpy.ndarray)
        data=np.empty((nb_snapshot,nb_atoms,len(heading)-2))
        
        # Find the index of the position x in the list called heading defined earlier
        index_pos=0
        while heading[index_pos] != 'x':
            index_pos+=1
        # We have to ignore the two first elements : 'ITEM:' and 'ATOMS'
        index_pos = index_pos - 2
        
        # prog is a variable used to display the percentage of the following loop that has been run        
        prog=0        
        
        i=0
        while i < nb_snapshot:
            current_dump=cut_name[0]+str(i*inter)+cut_name[1]
            f=open(current_dump,'r')
            data[i]=np.loadtxt(current_dump,skiprows=9)
            f.close()
            
            # We now apply a translation of all the points so that the origin be 0
            k=0
            while k < nb_atoms:
                data[i][k][index_pos] = data[i][k][index_pos] - box_dim[1]
                data[i][k][index_pos+1] = data[i][k][index_pos+1] - box_dim[1]
                data[i][k][index_pos+2] = data[i][k][index_pos+2] - box_dim[1]
                
                k+=1
                
            i+=1
            
            if i%int(0.1*nb_snapshot) == 0:
                prog+=1
                #print('Progression of dumplog:',prog*10,'\b%')
        
        
        
        # Attributes of the dumplog class
        self.heading=heading[:]     # Is there a difference with self.heading=heading ?
        self.nb_atoms=nb_atoms     # Number of atoms present in the simulation
        self.step_max=step_max     # Total number of steps in the simulations
        self.time_step=time_step     # Time step of the simulation in ps
        self.temperature=temperature    # Temperature of the simulation in K. If no equilibration phase has been run in npt, the output is the string: 'not fixed'
        self.time_real=time_real     # Real time which the simulation has last in ps
        self.inter=inter     # Number of time steps between each snapshot
        self.lattice=lattice    # Gives the nature of the lattice
        self.nb_snapshot=nb_snapshot     # Number of snapshots in the dump file
        self.box_dim=box_dim[:]     # Contains [length of the box,departure lx, end lx]
        self.dump_name=dump_name
        self.log_name=log_name
        self.data=data     # Do we need to write self.data=numpy.copy(data) instead ?
        










    def highest_disorder(self,peng_name,threshold='npt_standard'):
        '''
        The purpose of this function is to return the index of the snapshot for
        which the disorder is at its maximum.
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
        will be directly extracted through the self variable which has to be a dumplog
        class. If this condition is not respected, the user is invited to provide
        his own threshold.
        The treshold in the standard condition is calculated using the following
        formula : threshold=10*((3/2)*k*T) with k, the Boltzmann constant and T the 
        temperature extracted from self (dumplog class). It is thus advised to the
        user to estimate the temperature and to use a similar kind of threshold 
        to get consistent results. 
        The cohesive energy is estimated using the average of the potential energy
        of the first snapshot (before that the cascade starts really) even if it
        actually estimates the cohesive energy plus the average kinetic energy.
        self is a dumplog variable (see the first class definition).
        threshold must be provided in eV if it is not let as 'npt_standard'.
        Finally note that the user must provide the name which is used to
        characterise the column in the dump file where the potential energy of the
        atoms is stored.
        '''
        
        
        # Attribution of the Boltzmann constant
        k=scipy.constants.k
        
        T=self.temperature
        
        # Joule -> eV conversion
        J_to_eV=scipy.constants.physical_constants['joule-electron volt relationship'][0]
        
        
        if threshold == 'npt_standard':
            # Be careful to the temperature: extract it in dumplog class
            threshold = 5*((3/2)*k*T)
            # Conversion in eV
            threshold = J_to_eV * threshold
        
        # We need an average of the (cohesive energy + average kinetic energy) using the first snapshot
        # First we extract the index of the potential energy information
        index_peng=index_heading(self,peng_name)
        
        # Now we compute the average of the (cohesive energy + average kinetic energy)
        tot=0
        j=0
        while j<self.nb_atoms:
            epot=self.data[0][j][index_peng]
            tot+=epot
            j+=1
        
        average=tot/self.nb_atoms
        
        
        # Number max of atoms with potential energy greater than average + threshold
        above_thresh=0
        
        # Declaration of the variable which will be used to contain the index of the snapshot with the highest disorder
        snap=0
        
        k=0
        while k < self.nb_snapshot:
            l=0
            
            # Declaration of the counter to count the number of atoms with potential energy greater than the average + threshold
            count=0
            
            while l < self.nb_atoms:
                epot=self.data[k][l][index_peng]
                if epot > average + threshold:
                    count+=1
                l+=1
            
            if count > above_thresh:
                above_thresh=count
                snap=k
            
            k+=1
        
        # We now want to calculate the corresponding time using the snap given by the previous operations
        
        time = round(snap * self.inter * self.time_step,12)
        
        
        
        return [time,snap,above_thresh]






















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