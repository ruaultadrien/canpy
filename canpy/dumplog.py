# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 14:49:26 2016

@author: Adrien
"""

import numpy as np
import scipy.constants


from canpy.point_defects import point_defects
from canpy.useful import index_heading,time_to_snap,snap_to_time

from mpl_toolkits.mplot3d import axes3d
from matplotlib import animation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import juggle_axes







class dumplog:
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
        fdump=open(dump_name,'r')
        
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
        # Declaration of the vector used to store the data read from the text file
        data_copy = np.empty((nb_snapshot*nb_atoms,len(heading)-2))
        
        f=open(dump_name,'r')
        
        
        k=0
        line='0'
        while line:
            line=f.readline()
            if line == 'ITEM: TIMESTEP\n':
                i=0
                while i<9:
                    line=f.readline()
                    i+=1
                    
            list_string=line.split()
            i=0
            for a in list_string:
                data_copy[k,i]=float(a)
                i+=1
            
            # Following statement used to display the progression of the reading of the dump file
            if k == 0:
                print('Progression of dumplog:')
            if k % int(0.2*nb_atoms*nb_snapshot)==0:
                print(str(round(100*k/(nb_atoms*nb_snapshot)))+'%')
            
            k+=1
            
        
        f.close()
        
        
        # Declaration of a numpy array which will contain the same elements as data_copy but with another dimension to rank the data according to the snapshots
        data = np.empty((nb_snapshot,nb_atoms,len(heading)-2))
        
        # Reorganisation of the data_copy array through the use of a loop and using the new array called data
        data[0,:,:] = data_copy[0:nb_atoms,:]
        
        i=1
        while i < nb_snapshot:
            data[i,:,:] = data_copy[i*(nb_atoms):(i+1)*nb_atoms,:]
            i+=1
        
        
        
        # Find the index of the position x in the list called heading defined earlier
        index_pos=0
        while heading[index_pos] != 'x':
            index_pos+=1
        # We have to ignore the two first elements : 'ITEM:' and 'ATOMS'
        index_pos = index_pos - 2         
        
        # We now apply a translation of all the points for the origin of the frame to be set to 0
        i=0
        while i < nb_snapshot:
            k=0
            while k < nb_atoms:
                data[i][k][index_pos] = data[i][k][index_pos] - box_dim[1]
                data[i][k][index_pos+1] = data[i][k][index_pos+1] - box_dim[1]
                data[i][k][index_pos+2] = data[i][k][index_pos+2] - box_dim[1]
                
                k+=1 
            i+=1
        
        
        # Calculation of the lattice parameter
        if lattice == 'bcc':
            # Length of the box
            length=box_dim[0]        
                
            # Determination of the number of cubic cells per edge of the simulation box
            num_edge = round((nb_atoms/2)**(1/3))
            #print('num_edge:',num_edge)
                
            # Determination of the lattice parameter
            lat_param = length / num_edge
        
        elif lattice == 'fcc':
            # Length of the box
            length=box_dim[0]        
                
            # Determination of the number of cubic cells per edge of the simulation box
            num_edge = round((nb_atoms/4)**(1/3))
            #print('num_edge:',num_edge)
                
            # Determination of the lattice parameter
            lat_param = length / num_edge
        
        else:
            print('Error: the lattice variable in __init__ of dumplog is neither bcc nor fcc')
        
        
        
        
        
        
        
        
        # Attributes of the dumplog class
        self.heading=heading[:]     # Is there a difference with self.heading=heading ?
        self.nb_atoms=nb_atoms     # Number of atoms present in the simulation
        self.step_max=step_max     # Total number of steps in the simulations
        self.time_step=time_step     # Time step of the simulation in ps
        self.temperature=temperature    # Temperature of the simulation in K. If no equilibration phase has been run in npt, the output is the string: 'not fixed'
        self.time_real=time_real     # Real time which the simulation has last in ps
        self.inter=inter     # Number of time steps between each snapshot
        self.lattice=lattice    # Gives the nature of the lattice
        self.lat_param=lat_param     # lattice parameter in Angström
        self.nb_snapshot=nb_snapshot     # Number of snapshots in the dump file
        self.box_dim=box_dim[:]     # Contains [length of the box,departure lx, end lx]
        self.dump_name=dump_name
        self.log_name=log_name
        self.data=data     # Do we need to write self.data=numpy.copy(data) instead ?
        










    def highest_disorder(self,peng_name='c_peng',epot_thresh='npt_standard'):
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
        
        
        if epot_thresh == 'npt_standard':
            # Be careful to the temperature: extract it in dumplog class
            epot_thresh = 10*((3/2)*k*T)
            # Conversion in eV
            epot_thresh = J_to_eV * epot_thresh
        
        # We need an average of the (cohesive energy + average kinetic energy) using the first snapshot
        # First we extract the index of the potential energy information
        index_peng=index_heading(self,peng_name)
        #print('index_peng:',index_peng)
        
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
                if k==100 and l%1000==0:
                    print('epot:',epot)
                if epot > average + epot_thresh:
                    count+=1
                l+=1
            
            if count > above_thresh:
                above_thresh=count
                snap=k
            
            k+=1
        
        # We now want to calculate the corresponding time using the snap given by the previous operations
        time = snap_to_time(snap,self)
        
        
        # return [time in ps corresponding to highest disorder, the corresponding snapshot, the number of atoms found above the energy threshold for this snapshot]
        return {'time':time,'snap':snap,'atoms_above_thresh':above_thresh}

    
    
    
    
    
    
    
    def final_defects(self):
        '''
        Gives the configuration of the defects using the point_defects object.
        To find the final configuration we use the last snapshot available.
        '''
        
        return point_defects(self,self.time_real)







    def defects_movie(self,t_start=0,t_end='end',speed=1,peng_name='c_peng'):
        '''
        This method aims to create a movie of the defects produced during the 
        cascade from time t_start to t_end (in ps). t_start and t_end can be
        provided by the user or let as default value: from the beginning to 
        the end.
        Note that the user can control the speed of the movie by using the 
        variable speed. The standard speed is set to 1 and the user can increase
        or decrease it by respectively setting speed to a factor greater or
        lower than 1.
        '''
        
        
        if t_end == 'end':
            snap_end = self.nb_snapshot-1
        else:
            snap_end = time_to_snap(t_end,self)
        
        snap_start = time_to_snap(t_start,self)
        
        # The method that we are using in animates has a problem when there is no defects (with the first step typically) so we need to avoid it for the moment
        if snap_start == 0:
            snap_start+=1
            
        # We now want to find the snap with the highest disorder to then extract one relevant transition vector to center the defects in the middle of the graph
        high_dis = self.highest_disorder(peng_name=peng_name)
        time_high_dis=high_dis['time']
        defects = point_defects(self,time_high_dis)
        tran = defects.transi_center()['tran']
        
        
        fig = plt.figure()
        #ax = fig.gca(projection='3d')
        ax = fig.add_subplot(111, projection='3d')
            
        
        pi = ax.scatter([],[],[], label='Interstitials', c='r', marker='^')
        pv = ax.scatter([],[],[], label='Vacancies', c='b', marker='o')
            
        ax.set_title('Movie of the defects from '+str(t_start)+' ps to '+str(snap_to_time(snap_end,self))+' ps')       
        
        ax.set_xlabel('x axis[Angström]')
        ax.set_ylabel('y axis[Angström]')
        ax.set_zlabel('z axis[Angström]')
        
        ax.set_xlim3d(0,self.box_dim[0])
        ax.set_ylim3d(0,self.box_dim[0])
        ax.set_zlim3d(0,self.box_dim[0])
        
        ax.legend(frameon = True, fancybox = True, ncol = 1, fontsize = 'x-small', loc = 'lower right')
        
        
        def animate(i):
            
            gamma = point_defects(self,snap_to_time(snap_start+i,self))
            
            gamma_transited = gamma.transi(tran)
            
            inter = gamma_transited['inter'][:]
            vac = gamma_transited['vac'][:]
            
            
            # 3 arrays containing the coordinates of the interstitials
            xi,yi,zi = np.transpose(inter)
            
            # 3 arrays containing the coordinates of the vacancies
            xv,yv,zv = np.transpose(vac)
            
            pi._offsets3d = juggle_axes(xi,yi,zi,'z')
            pv._offsets3d = juggle_axes(xv,yv,zv,'z')
            
            
            '''
            ax.clear()
            ax.scatter(xi,yi,zi, label='Interstitials', c='r', marker='^')
            ax.scatter(xv,yv,zv, label='Vacancies', c='b', marker='o')
    
            
            ax.set_title('Movie of the defects from '+str(t_start)+' ps to '+str(snap_to_time(snap_end,self))+' ps')       
            
            ax.set_xlabel('x axis[Angström]')
            ax.set_ylabel('y axis[Angström]')
            ax.set_zlabel('z axis[Angström]')
            
            ax.set_xlim3d(0,self.box_dim[0])
            ax.set_ylim3d(0,self.box_dim[0])
            ax.set_zlim3d(0,self.box_dim[0])
            
            ax.legend(frameon = True, fancybox = True, ncol = 1, fontsize = 'x-small', loc = 'lower right')
            '''
            
            
            
            # Display of the progression
            if i % round(0.2*(snap_end-snap_start+1))==0 and i // round(0.2*(snap_end-snap_start+1)) != 0:
                if i // round(0.2*(snap_end-snap_start+1)) == 1:
                    print('Progression of defects_movie:')
                print(str(round(100*i/(snap_end-snap_start+1)))+'%')
                
            # We want to close the window to let the program continue working (bug when we run the program in cmd)
            #if i == snap_end-snap_start:
            #    plt.close(fig)
            
        # Animate
        anim = animation.FuncAnimation(fig, animate, interval = int(12000/speed),
                                       frames=snap_end-snap_start+1)
        #init_func=init
        # blit=True
                                       
        # Save
        anim.save('defects_movie.mp4', dpi=150, extra_args=['-vcodec', 'libx264'])
        # !! Think to find a way to change the dpi to make a quality video in full screen
        
        print('wowowo')
        
        
        















'''
def index_heading(dl,name):
    This function is made to provide the index of the element corresponding to
    the name (string) provided in the list called heading which is a variable 
    output of the dumplog class.
    Therefore the user has to provide a dumplog object and the name of one of
    the columns of the dump file and the function returns the index for which
    this name (string) is stored in the heading list.
    
    
    # Definition of the index which starts at 0
    i=0
    
    while dl.heading[i] != name:
        i+=1
    # We have to ignore the two first elements : 'ITEM:' and 'ATOMS'
    i=i-2
    return i
    
    '''