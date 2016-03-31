# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 18:02:59 2016

@author: Adrien
"""


from canpy.useful import index_heading,time_to_snap

import numpy as np
import scipy.constants
#from index_heading import index_heading



class mapping:
    '''
    This class divides the simulation box in cells and returns a numpy array
    with the coordinates of the centre of each cell and some informations about
    the atoms contained in the cells.
    The informations contained about each cell in the numpz arraz are: 1) the
    sum of the potential energy of all the atoms contained in the cell and
    2) the number of atom contained in the cell
    The user must provide a dumplog class object called dl
    The user provides also nb_cut which is used to define in which number the 
    edges of the box will be divides. If the user defines nb_cut as 3 then the 
    edges are divided in 3 and the whole simulation box is divided into 27 cells.
    Finally the user has to provide as argument the time in ps at which he wants
    the mapping to be created.
    Note that this function has been written such that the simulation box is
    cubic, some slight modifications may be required to apply it to non cubic
    simulation boxes.
    Note also that to run this function we need the potential energy to be computed
    in the LAMMPS simulation and we need the user to provide the name used to
    store the potential energy in the dumplog file.
    '''

    def __init__(self,dl,nb_cut,time,peng_name):
        
        # We calculate the side length of the cells using the length of the simulation box and nb_cut
        a_cell= dl.box_dim[0] / nb_cut
        
        # We now want to find the snapshot corresponding to the time provided by the user
        snap=time_to_snap(time,dl)
        
        # Total number of cells
        nb_cells=(nb_cut)**3
        
        # Initialisation of a numpy array with 0 which will contain all the informations about the cells (coordinates of centre / sum of potential energy / nb of atoms in the cell)
        cell=np.zeros((nb_cells,6))
        
        # Index which will be used to identify the cells
        m=0
        
        # 3 loops to fill the cell numpy array with the coordinates of the centre of each of the cells
        i=0
        while i < nb_cut:
            
            j=0
            while j < nb_cut:
                
                k=0
                while k < nb_cut:
                    cell[m][0] = (a_cell/2) + i*a_cell
                    cell[m][1] = (a_cell/2) + j*a_cell
                    cell[m][2] = (a_cell/2) + k*a_cell
                    
                    m+=1
                    k+=1
                
                j+=1
            
            i+=1
        
        # We will now going through all the atoms using the numpy array data provided by the dl object and complete the remaining informations in the cell numpy array (sum of potential energy and number of atoms)
        # First we need the index of the x position in the heading list contained in the dumplog object dl
        index_pos=index_heading(dl,'x')
        
        # Then we need the index of the potential energy
        index_peng=index_heading(dl,peng_name)
        
        i=0
        while i < dl.nb_atoms:
            # For each atom we extract xc, yc and zc which will then be used to find the index of the cell where the atom is located
            # Three if statements are added to be sure that if a_cell is not precise enough: xc, yc, zc will never overcome nb_cut which would involve an IndexError
            xc = int(dl.data[snap][i][index_pos] // a_cell)
            if xc >= nb_cut:
                xc=nb_cut-1
            elif xc < 0:
                xc=0
            
            yc = int(dl.data[snap][i][index_pos + 1] // a_cell)
            if yc>=nb_cut:
                yc=nb_cut-1
            elif yc < 0:
                yc=0
            
            zc = int(dl.data[snap][i][index_pos + 2] // a_cell)
            if zc>=nb_cut:
                zc=nb_cut-1
            elif zc < 0:
                zc=0
            
            # Using the 3 previous varibles (xc, yc and zc) we now find the cell index (m) where the atom is found
            m = int(xc*(nb_cut**2) + yc*(nb_cut) + zc)
            
            # Now that we have the m index we can complete the informations in the cell numpy array
            cell[m][3]+=dl.data[snap][i][index_peng]
            cell[m][4]+=1
            
            i+=1
        
        # Compute the average potential energy for each of the cells
        i=0
        while i < nb_cells:
            cell[i][5] = cell[i][3] / cell[i][4]
            i+=1
        
        # Control if every atom has been assigned to a cell exactly once
        counter=0
        i=0
        while i < nb_cells:
            counter += cell[i][4]
            i+=1
            
        if abs(counter-dl.nb_atoms) > 0.1:
            print('Problem in mapping() function, number of atoms attributed to cells is wrong')
        
        self.cell = cell
        
        # Declararation of a varibale containing the side length of the box to be passed more easily in the method called volume and declared below
        self.length = dl.box_dim[0]
        
        # Temperature of the npt simulation if an equilibration phase is run before
        self.temperature = dl.temperature









    def volume(self,epot_thresh='npt_standard'):
        '''
        The idea of this function is to return the volume of the cascade in Angstroms, 
        using the output of the mapping() function. The mapping function divides 
        the simulation box in cells defined by the user and it returns a numpy 
        array with the average potential energy of each of the cells. All the cells
        with an average potential energy greater than the threshold defined by the 
        user are considered as beeing part of the cascade and all the others are considered as outside
        the cascade. We can thus easily calculate the volume of the cascade, knowing
        the volume of a cell and knowing the number of cells that constitute the
        cascade.
        dl is a dumplog object
        map_out is the output of the mapping() function
        epot_thresh is the average potential energy threshold above which a cell is
        considered as beeing part of the volume.
        '''
        
        # In this part we define the threshold if it has been set as the standard one
        
        if epot_thresh == 'npt_standard':
            # Attribution of the Boltzmann constant
            k=scipy.constants.k
            
            T=self.temperature
            
            # Joule -> eV conversion
            J_to_eV=scipy.constants.physical_constants['joule-electron volt relationship'][0]

            # Be careful to the temperature: extract it in dumplog class
            epot_thresh = 5*((3/2)*k*T)
            # Conversion in eV
            epot_thresh = J_to_eV * epot_thresh
        
        
        
        # Calculation of the volume of a cell
        # Number of cells in which the simulation boy has been divided
        nb_cells = np.shape(self.cell)[0]
        
        # Volume of the simulation box
        vol_box = (self.length)**3
        
        # Volume of an individual cell
        vol_cell = vol_box / nb_cells
        
        # The following loop is made to count the number of cells with an average potential energy greater than the threshold (epot_thresh)
        # Counter of cells respecting the previous stated condition
        cell_count=0
        
        i=0
        while i < nb_cells:
            
            if self.cell[i][5] > epot_thresh:
                cell_count+=1
            
            i+=1
        
        # Computation of the cascade volume
        vol_cascade = cell_count * vol_cell
        
        return vol_cascade




    
    
    
    
    
    
    