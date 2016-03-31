# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 09:53:50 2016

@author: Adrien
"""




def split(name):
    '''
    This function can be used by the user if the dump output of LAMMPS is 
    contained in only one file. This kind of file cannot be treated by the
    dumplog class so the split function can be used to divide the dump file
    in as many files as the snapshots asked by the user in the LAMMPS input
    file when he called the dump command.
    The file name can't have more than one dot otherwise the programme fails
    '''
    
    
    # This part aims to extract the left and right part of the name of the dump file provided
    index=name.find('.')
    
    left=name[:index]
    right=name[index:]
    
    
    # Now looking for the number of atoms in the simulation to have the information about the number of lines in each snapshots
    fdump=open(name,'r')
    
    i=0
    while i<4:
        line=fdump.readline()
        i+=1
    
    nb_atoms=int(line)
    
    # We now return to the beginning of the file for breaking it down
    fdump.seek(0,0)
    
    # Copying of the different snapshots from the big unique file in different file with their own snapshot
    # test variable allows us to go out of the while
    test=True
    
    while test:
        
        # The if/else statements test the first line of the snapshot, if we reach EOF we directly exit the while-loop        
        line=fdump.readline()
        if line != '':
            line=fdump.readline()
            
            # Extraction of the name of the new file + opening of the file which will contain the text related to the current snapshot
            step=int(line)
            new_name=left+str(step)+right
            fnew=open(new_name,'w')
            
            # Write the two first lines of the snapshot description on the new file
            fnew.write('ITEM: TIMESTEP\n')
            fnew.write(line)
            
            # Write the rest of the lines in the new file
            j=0
            while j < nb_atoms+7:
                line=fdump.readline()
                fnew.write(line)
                j+=1
            
            # Close the new file
            fnew.close()
            
    
        else:
            test=False
    
    # Close the initial dump file to end the function processings
    fdump.close()
    