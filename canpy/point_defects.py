# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 16:42:41 2016

@author: Adrien
"""

import numpy as np

import math

from canpy.useful import index_heading,time_to_snap

from mpl_toolkits.mplot3d import axes3d
from matplotlib import animation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import juggle_axes





class point_defects:
    '''
    point_defect is a class with a set of method to analyse the defects produced
    during the cascade at a time t.
    '''

    def __init__(self,dl,time):
        '''
        The idea here is to identify the point defects produced during the simulation
        cascade. By point defects we mean interstitial atoms and vacancies. 
        For doing this we need the user to provide a dumplog object as argument as
        well as the time in ps at which we want to know the defect configuration.
        The idea is to generate a network with the theoretical coordinates of each
        atomic position in a perfect lattice. Then we want to check each atom and
        identify the closest point to this atom in the network of theoretical points
        that we have just generated. if the distance between the atom and the closest
        corresponding point of the perfect lattice is greater than half of the minimal
        distance between two atoms in the lattice then this atom is considered as 
        interstitial and its coordinates are stored in a numpy array. 
        If the distance between the theoretical point and the atom checked is lower
        than half of the minimal distance between two atoms then the theoretical point
        of the lattice is considered as occupied by the atom.
        If the atom checked is considered as located at a theoretical point of the
        lattice where another atom is already supposed to be then we consider the
        closest atom to the point of the lattice as occupying the latter and the
        other one is considered as an interstitials and again its coordinates are
        stored in a numpy array.
        At the end of the process we determine the vacancies as the points of the
        theoretical lattice which haven't been associated with another atom.
        This function is based on the assumption that because of the thermal activation
        of the atoms the latter are never at the exact emplacement of their
        theoretical site in the lattice. 
        This function present the weakness to not work for systems which haven't been
        equilibriated before the cascade simulation to be run. Indeed it would lead
        for the atoms exact     lat_param
        '''
        
        self.lat_param = dl.lat_param
        self.length = dl.box_dim[0]
        self.lattice = dl.lattice
        self.time = time
        
        if dl.lattice == 'bcc':
            # bcc_point_defects() is a function defined below in this script
            self.defects = bcc_point_defects(dl,time)
                
        elif dl.lattice == 'fcc':
            # fcc_point_defects() is a function defined below in this script
            self.defects = fcc_point_defects(dl,time)
            
        else:
            print('Error')
            return 0


    
    
    
    
    
    
    
    def transi(self,tran):
        '''
        This function allows us to apply a linear transition to the interstitials
        and vacancies in the simulation box. The transition must be applied on
        a point_defects object. The user is free to provide any transition vector.
        The periodic boundaries are taken into account and a defects translated
        outside the boundaries of the box will be reflected through the periodic
        boundaries.
        tran must be a one dimension numpy array containing the coordinates
        of the transition in the following order: x,y,z
        '''
        
        length = self.length
        
        inter = self.defects['inter'][:]
        vac = self.defects['vac'][:]
        
        inter = inter + tran
        m=0
        while m < self.defects['FP_count']:
            n=0
            while n < 3:
                a = inter[m,n]
                if a > length:
                    a = a - length
                elif a < 0:
                    a = a + length
                            
                inter[m,n]=a
                n+=1
                        
            m+=1
        
        
        
        vac = vac + tran
        m=0
        while m < self.defects['FP_count']:
            n=0
            while n < 3:
                a = vac[m,n]
                if a > length:
                    a = a - length
                elif a < 0:
                    a = a + length
                            
                vac[m,n]=a
                n+=1
                        
            m+=1
        
        
        return {'inter':inter,'vac':vac}
    
    
    
    
    







    def transi_center(self):
        '''
        The idea of this method is to return the same kind of output as the
        point_defects class but with a transition of the coordinates of the
        vacancies and interstitials such that their center of mass is centered
        in the middle of the box. It would allow us to have a representation
        of the defects without any reflections of some of them through the
        periodic boundaries of the simulation box. 
        To do that we will carry out a minimisation of the standard deviation
        of the position of the vacancies for different transition vectors. This
        will give us a set of transition vectors for which the standard deviation 
        of the position is minimised. The transition vector that we will take 
        is the one for which the sum the CDM vector gives the closest point to 
        the middle of the box. 
        We choose the vacancies because the interstitials are spread out on
        a large volume which means that an interstitial reflected through a
        periodic boundary could participate to reduce the standard deviation
        and would lead to the failure of the method. The vacancies are supposed
        to be gathered enough so that a vacancies at the wrong place due to
        a reflection through a periodic boundary would nearly necesseraly lead
        to an increase of the standard deviation.
        '''
        
        if abs(self.time) < 1e-10:
            return {'inter':np.array([0,0,0]),'vac':np.array([0,0,0]),'tran':np.array([0,0,0])}
        
        
        length = self.length
        
        # Definition of an increment for constructing the different transition vectors that we will try
        incr = 0.1*length
        
        # Storage of the numpy array containing the coordinates of the vacancies
        vac = self.defects['vac'][:]         
        
        # Initialisation of the transition vector
        tran = np.zeros(3)
        
        # Creation of a numpy array supposed to store the vac vector after transition
        transited = np.zeros_like(vac)
        
        # Creationof a numpy array which will store all the tran and com vectors
        store = np.zeros((2,8000,3))
        
        # Creation of a numpy array to store the coordinates of the center of mass of the transited array
        com = np.zeros(3)
        
        # Initialisation of a variable used to check which tran vectors are considered as giving the least standard deviation from the com
        mini = length**4 # initialisation to be sure that mini is not too small
        
        # Initialisation of an index list which will contain all the index of the array called store for which the standard deviation is minimised       
        store_indexlist=[]        
        
        # Initialisation of the counters for the following triple loop
        i=-4
        j=-4
        k=-4
        
        p=0
        
        while i <= 5:
            tran[0] = i*incr
            
            j=0
            while j <= 5:
                tran[1] = j*incr
                
                k=0
                while k <= 5:
                    tran[2] = k*incr
                    store[0,p]=tran[:]
                    
                    # We apply the transition
                    transited = vac + tran
                    
                    # We want to be sure that the translated vacancies are not outside the boundaries of the box, if so we reflect the vacancies through the periodic boundaries
                    m=0
                    while m < self.defects['FP_count']:
                        n=0
                        while n < 3:
                            a = transited[m,n]
                            if a > length:
                                a = a - length
                            elif a < 0:
                                a = a + length
                            
                            transited[m,n]=a
                            n+=1
                        
                        m+=1
                            
                    # Calculation of the center of mass of the transited array
                    com[0]=np.average(transited[:,0])
                    com[1]=np.average(transited[:,1])
                    com[2]=np.average(transited[:,2])
                    
                    store[1,p]=com[:]
                    
                    # We want to calculate the distance between every transited vacancy and the com
                    coor_dist = transited - com
                    
                    stalidev=0
                    for a in coor_dist[:]:
                        stalidev += np.linalg.norm(a)**2
                        
                    stalidev = stalidev / self.defects['FP_count']
                    
                    if abs(stalidev - mini) < 1e-4:
                        store_indexlist+=[p]
                    if abs(stalidev - mini) >= 1e-4 and stalidev < mini:
                        mini = stalidev
                        store_indexlist=[p]
                    
                    
                    p+=1
                    k+=1
                
                j+=1
            
            i+=1
        
        
        # Definition of a numpy array which contains the coordinates of the center of the box
        center = np.array([length/2,length/2,length/2])
        
        # Initialisation of a variable to find the minimal distance in the following for loop
        mini_dist= length**3 # Initialised with an arbitrary high value        
        
        for a in store_indexlist:
            com=store[1][a]
            dist = np.linalg.norm(com-center)
            
            if dist < mini_dist:
                mini_dist=dist
                tran_index = a
                
        # Knowing the index corresponding to the right tran vector in the store array we now want to transform the arrays containing the vac and inter
        tran = store[0][tran_index]
        
        
        
        vac = self.transi(tran)['vac']
        inter = self.transi(tran)['inter']
        
        
        '''
        vac = vac + tran
        
        # We want to be sure that the translated vacancies are not outside the boundaries of the box, if so we reflect the vacancies through the periodic boundaries
        m=0
        while m < self.defects['FP_count']:
            n=0
            while n < 3:
                a = vac[m,n]
                if a > length:
                    a = a - length
                elif a < 0:
                    a = a + length
                            
                vac[m,n]=a
                n+=1
                        
            m+=1
        
        
        inter = self.defects['inter'][:] + tran
        
        # We want to be sure that the translated vacancies are not outside the boundaries of the box, if so we reflect the vacancies through the periodic boundaries
        m=0
        while m < self.defects['FP_count']:
            n=0
            while n < 3:
                a = inter[m,n]
                if a > length:
                    a = a - length
                elif a < 0:
                    a = a + length
                            
                inter[m,n]=a
                n+=1
                        
            m+=1
        '''
        
        
        return {'inter':inter,'vac':vac,'tran':tran}
            
            
        
    
        
        
        
        
        


    def plot_original_frame(self,name='plot_original_frame'):
        '''
        This method is used to plot the vacancies and interstitials using the
        original frame, that is to say the one that we have considered during
        the LAMMPS simulation. This method of plotting make the morphology of 
        the cascade hard to visualise because the frame is not adapted and
        the defects often appear far from each other due to the reflections
        through the periodic boundaries.
        '''    
        
        xi = self.defects['inter'][:,0]
        yi = self.defects['inter'][:,1]
        zi = self.defects['inter'][:,2]
        
        xv = self.defects['vac'][:,0]
        yv = self.defects['vac'][:,1]
        zv = self.defects['vac'][:,2]
        

        fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        ax = fig.gca(projection='3d')
    
        ax.scatter(xi, yi, zi, c='r', marker='^')
        ax.scatter(xv, yv, zv, c='b', marker='o')
    
        ax.set_xlabel('x axis[Angström]')
        ax.set_ylabel('y axis[Angström]')
        ax.set_zlabel('z axis[Angström]')
        
        ax.set_xlim3d(0,self.length)
        ax.set_ylim3d(0,self.length)
        ax.set_zlim3d(0,self.length)
        
        plt.show()
        fig.savefig(name)
        plt.close(fig)
    
    
    
    
    
    
    
    
    def plot_centered(self,name='plot_centered'):
        '''
        The idea of this function is to plot the defects with a translation of
        the frame. Indeed we want to use the center of mass of the vacancies
        as reference and put it at the center of the graph. We thus want to
        translate all the defects shuch that the center of mass of vacancies
        is found at the center of the box (1/2 of box length, 1/2 of box length,
        1/2 of box length)
        We here make the assumption that the vacancies which have coordinates
        greater than half 
        '''

        inter = self.transi_center()['inter'][:]
        vac = self.transi_center()['vac'][:]
        
        
        # Plotting of the translated vacancies and interstitials
        xi = inter[:,0]
        yi = inter[:,1]
        zi = inter[:,2]
        
        
        xv = vac[:,0]
        yv = vac[:,1]
        zv = vac[:,2]


        fig = plt.figure()
        ax = fig.gca(projection='3d')
    
        ax.scatter(xi, yi, zi, label='Interstitials', c='r', marker='^')
        ax.scatter(xv, yv, zv, label='Vacancies', c='b', marker='o')
        
        #ax.scatter(x_clus,y_clus,z_clus, c='r', marker='^')
        
        ax.set_title('Interstitials and vacancies at '+str(self.time)+' ps')       
        
        ax.set_xlabel('x axis[Angström]')
        ax.set_ylabel('y axis[Angström]')
        ax.set_zlabel('z axis[Angström]')
        
        ax.set_xlim3d(0,self.length)
        ax.set_ylim3d(0,self.length)
        ax.set_zlim3d(0,self.length)
        
        ax.legend(frameon = True, fancybox = True, ncol = 1, fontsize = 'x-small', loc = 'lower right')
        
        plt.show()
        fig.savefig(name)
        plt.close(fig)











    def rot_frame_anim(self,speed=1):
        '''
        The idea of this method is to create an animated plot of the interstitials
        and vacancies. The animation involves a 360° rotation of the camera on 
        the (x,y) plane and so around the z-axis.
        rot_animation() is based on the same plot as the one displayed by the
        method plot_vac_centered().
        The user can change the speed of rotation by modifying the value of the
        speed argument. A value of 1 is set as standard. To slow down the rotation
        the user must set a value lower than 1 and to speed up: a value greater
        than 1 must be set.
        '''
        
        
        inter = self.transi_center()['inter'][:]
        vac = self.transi_center()['vac'][:]
        
        
        xi,yi,zi = np.transpose(inter)
        
        xv,yv,zv = np.transpose(vac)
        
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        
        def init():
            ax.scatter(xi, yi, zi, label='Interstitials', c='r', marker='^')
            ax.scatter(xv, yv, zv, label='Vacancies', c='b', marker='o')
    
            
            ax.set_title('Interstitials and vacancies at '+str(self.time)+' ps')       
            
            ax.set_xlabel('x axis[Angström]')
            ax.set_ylabel('y axis[Angström]')
            ax.set_zlabel('z axis[Angström]')
            
            ax.set_xlim3d(0,self.length)
            ax.set_ylim3d(0,self.length)
            ax.set_zlim3d(0,self.length)
            
            ax.legend(frameon = True, fancybox = True, ncol = 1, fontsize = 'x-small', loc = 'lower right')
            
        
        
        def animate(i):
            ax.view_init(elev=0,azim=i)
            
        # Animate
        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=360, blit=True)
        # Save
        anim.save('rot_frame_anim_.mp4', fps=30*speed, dpi=150, extra_args=['-vcodec', 'libx264'])
        
        plt.close(fig)







    def rot_fixframe_anim(self,speed=1):
        '''
        Produce a movie diplaying the defects rotating at 360° in a fix frame
        '''
            
        transited = self.transi_center()
            
        inter = transited['inter'][:]
        vac = transited['vac'][:]
        
        # numpy array used to center the coordinates around the origine so that we can apply a rotation
        length = self.length
        back_center = np.array([0.5*length,0.5*length,0.5*length])
        inter = inter - back_center
        vac = vac - back_center
        
        # Transpose inter and vac so that we can apply the rotation matrix defined later
        inter = inter.transpose()
        vac = vac.transpose()
        
        
        fig = plt.figure()
        #ax = fig.gca(projection='3d')
        ax = fig.add_subplot(111, projection='3d')

        pi = ax.scatter([],[],[], label='Interstitials', c='r', marker='^')
        pv = ax.scatter([],[],[], label='Vacancies', c='b', marker='o')
        
        ax.set_title('Interstitials and vacancies at '+str(self.time)+' ps')       
        
        ax.set_xlabel('x axis[Angström]')
        ax.set_ylabel('y axis[Angström]')
        ax.set_zlabel('z axis[Angström]')
        
        ax.set_xlim3d(0,length)
        ax.set_ylim3d(0,length)
        ax.set_zlim3d(0,length)
        
        ax.view_init(elev=0,azim=70)
        
        ax.legend(frameon = True, fancybox = True, ncol = 1, fontsize = 'x-small', loc = 'lower right')
        
        
        def animate(i,inter,vac,back_center,length):

            # declaration of the rotation matrix
            angle= i * ((2*math.pi)/360)
            rotation = np.matrix([[math.cos(angle),-math.sin(angle),0],
                                  [math.sin(angle),math.cos(angle),0],
                                  [0,0,1]])
            
            # Application of the matrix to rotate the defects and application of back_center to put them in the center of the frame 
            inter = np.array(rotation*inter)
            vac = np.array(rotation*vac)
            
            inter = inter.transpose() + back_center
            vac = vac.transpose() + back_center
            
            # We want to be sure that the rotated interstitials are not outside the boundaries of the box, if so we reflect the interstitials through the periodic boundaries
            m=0
            while m < self.defects['FP_count']:
                n=0
                while n < 3:
                    a = inter[m,n]
                    if a > length:
                        a = a - length
                    elif a < 0:
                        a = a + length
                            
                    inter[m,n]=a
                    n+=1
                        
                m+=1
            
            # We want to be sure that the rotated vacancies are not outside the boundaries of the box, if so we reflect the vacancies through the periodic boundaries
            m=0
            while m < self.defects['FP_count']:
                n=0
                while n < 3:
                    a = vac[m,n]
                    if a > length:
                        a = a - length
                    elif a < 0:
                        a = a + length
                            
                    vac[m,n]=a
                    n+=1
                        
                m+=1
            
            
            
            # 3 arrays containing the coordinates of the interstitials
            xi,yi,zi = inter.transpose()
            
            
            # 3 arrays containing the coordinates of the vacancies
            xv,yv,zv = vac.transpose()
            
            pi._offsets3d = juggle_axes(xi,yi,zi,'z')
            pv._offsets3d = juggle_axes(xv,yv,zv,'z')
        
        
        # Animate
        anim = animation.FuncAnimation(fig, animate, fargs=(inter,vac,back_center,length), frames=360)
        #init_func=init
        # blit=True
                                       
        # Save
        anim.save('rot_fixframe_movie.mp4', fps=30*speed, dpi=150, extra_args=['-vcodec', 'libx264'])
        # !! Think to find a way to change the dpi to make a quality video in full screen
        
        plt.close(fig)
        
        
        





    def clusters(self,NN=1):
        
        if self.lattice == 'bcc':
            
            return {'inter':bcc_inter_clusters(self,NN),'vac':bcc_vac_clusters(self,NN)}
        
        elif self.lattice == 'fcc':
            
        
            return {'inter':fcc_inter_clusters(self,NN),'vac':fcc_vac_clusters(self,NN)}            
            
        else:
            print('Error: lattice is neither bcc nor fcc')
            return 0












def bcc_point_defects(dl,time):
    '''
    This function is just a subpart of the main one called point_defects.
    If the lattice is bcc then this function is called through the main one via
    the use of some if statements.
    Another subpart of the main function (point_defects) called fcc_point_defects
    is used if the lattice is found as beeing fcc. This function is define later
    in the script.
    '''
    
    # Here we define a variable giving the coresponding snapshot for the time provided
    n_snap=time_to_snap(time,dl)
    #print('n_snap:',n_snap)

    # Extract the heading indexes of the x position
    index_pos=index_heading(dl,'x')
    #print('index_pos:',index_pos)
    
    # Length of the box
    length=dl.box_dim[0]
        
    # Determination of the number of cubic cells per edge of the simulation box
    num_edge = round((dl.nb_atoms/2)**(1/3))
    #print('num_edge:',num_edge)
        
    # Determination of the lattice parameter
    lat_param = length / num_edge
    #print('lat_param',lat_param)
    
    # Determination of the minimal distance between 2 atoms
    mini_dist = (math.sqrt(3)/4)*lat_param
    
    # Counter to fill the informations about the vacancies and the interstitials in the defects numpy array
    inter_count = 0
    vac_count = 0
    
    # A bcc lattice is constituted of two interpenatrated cubic lattices, the 2nd one beeing a translation of (1/2a,1/2a,1/2a) from the first one
    # Initialisation of the 1st cubic lattice
    cub1=np.zeros((round(dl.nb_atoms/2),7))
    
    # Initialisation of the 2nd interpenetrating cubic lattice
    cub2=np.zeros((round(dl.nb_atoms/2),7))
        
    # Initialisation of the index used to identify the points of the lattice that we intend to generate        
    m=0        
    
    # 3 Loops to generate the points of the lattices by defining the 2 interpenetrating networks cub1 and cub2        
    i=0
    while i < num_edge:
    
        j=0
        while j < num_edge:
        
            k=0
            while k < num_edge:
                cub1[m][0] = i*lat_param
                cub1[m][1] = j*lat_param
                cub1[m][2] = k*lat_param
                
                cub2[m][0] = (lat_param/2) + i*lat_param
                cub2[m][1] = (lat_param/2) + j*lat_param
                cub2[m][2] = (lat_param/2) + k*lat_param
            
                m+=1
                k+=1
        
            j+=1
        
        i+=1
        

    
    # Initialisation of the numpy arrays which will contain the coordinates of interstitials and vacancies
    inter = np.zeros((int(dl.nb_atoms/2),3))
    vac = np.zeros((int(dl.nb_atoms/2),3))
    
    # Now we need to check all the atom through the dl.data array to attribute them to a lattice point or to define them as interstitials
    i=0
    while i < dl.nb_atoms:
        
        # We define the position variables of the atom to clarify the code
        x_at = dl.data[n_snap][i][index_pos]
        y_at = dl.data[n_snap][i][index_pos+1]
        z_at = dl.data[n_snap][i][index_pos+2]
        
        i+=1
            
        # We find the closest lattice point from this atom in cub1
        a1 = int(x_at // lat_param)
        b1 = int(y_at // lat_param)
        c1 = int(z_at // lat_param)
            
            
        if abs(a1*lat_param - x_at) > abs((a1+1)*lat_param - x_at):
            a1+=1
            
        if abs(b1*lat_param - y_at) > abs((b1+1)*lat_param - y_at):
            b1+=1
            
        if abs(c1*lat_param - z_at) > abs((c1+1)*lat_param - z_at):
            c1+=1
            
        a1 = a1%num_edge
        b1 = b1%num_edge
        c1 = c1%num_edge
            
        # Using a, b and c we find the index in the numpy array cub1 of the closest lattice point to the atom
        index_cub1 = int(a1*(num_edge**2) + b1*(num_edge) + c1)
        

        # We find the closest lattice point from this atom in cub2
        a2 = int((x_at - lat_param/2) // lat_param)
        b2 = int((y_at - lat_param/2) // lat_param)
        c2 = int((z_at - lat_param/2) // lat_param)
            
            
        if abs((a2 + 1/2)*lat_param - x_at) > abs((a2+1 + 1/2)*lat_param - x_at):
            a2+=1
        
        if abs((b2 + 1/2)*lat_param - y_at) > abs((b2+1 + 1/2)*lat_param - y_at):
            b2+=1
        
        if abs((c2 + 1/2)*lat_param - z_at) > abs((c2+1 + 1/2)*lat_param - z_at):
            c2+=1                
            
        a2 = a2%num_edge
        b2 = b2%num_edge
        c2 = c2%num_edge
            
        # Using a, b and c we find the index in the numpy array cub1 of the closest lattice point to the atom
        index_cub2 = int(a2*(num_edge**2) + b2*(num_edge) + c2)
            
            
        # Determination of which of the theoretical lattice points is the closest to the atom checked
        # First we need to calculate the distances between the theoretical point and the atom associated
        # We need to be careful that this distance is not greater than half of the length of the box
        
        # dist1 (use of the intermediate variables x_dist, y_dist and z_dist)
        x_dist=x_at
        y_dist=y_at
        z_dist=z_at
            
        if x_at-cub1[index_cub1][0] < (-length/2):
            x_dist = x_dist + length
        elif x_at-cub1[index_cub1][0] > length/2:
            x_dist = x_dist - length
        
        if y_at-cub1[index_cub1][1] < (-length/2):
            y_dist = y_dist + length
        elif y_at-cub1[index_cub1][1] > length/2:
            y_dist = y_dist - length
        
        if z_at-cub1[index_cub1][2] < (-length/2):
            z_dist = z_dist + length
        elif z_at-cub1[index_cub1][2] > length/2:
            z_dist = z_dist - length
        
        dist1 = math.sqrt( (x_dist-cub1[index_cub1][0])**2 + (y_dist-cub1[index_cub1][1])**2 + (z_dist-cub1[index_cub1][2])**2 )
        
        
        # dist2 (use of the intermediate variables x_dist, y_dist and z_dist)
        x_dist=x_at
        y_dist=y_at
        z_dist=z_at            
        
        if x_at-cub2[index_cub2][0] < (-length/2):
            x_dist = x_dist + length
        elif x_at-cub2[index_cub2][0] > length/2:
            x_dist = x_dist - length
        
        if y_at-cub2[index_cub2][1] < (-length/2):
            y_dist = y_dist + length
        elif y_at-cub2[index_cub2][1] > length/2:
            y_dist = y_dist - length
        
        if z_at-cub2[index_cub2][2] < (-length/2):
            z_dist = z_dist + length
        elif z_at-cub2[index_cub2][2] > length/2:
            z_dist = z_dist - length            
        
        dist2 = math.sqrt( (x_dist-cub2[index_cub2][0])**2 + (y_dist-cub2[index_cub2][1])**2 + (z_dist-cub2[index_cub2][2])**2 )
        
        lourdeur=1
        if lourdeur%100==0:
            print(lourdeur)
        
        # We determine which of the 2 theoretical lattices gives the closest point to the atom
        if dist1 < dist2:
                
            # We record the atom checked as an interstitial or not
            if (dist1 < mini_dist) and (abs(cub1[index_cub1][6]) < 1e-8):
                cub1[index_cub1][3] = x_at
                cub1[index_cub1][4] = y_at
                cub1[index_cub1][5] = z_at
                cub1[index_cub1][6] = dist1
                lourdeur+=1
                
            elif dist1 < mini_dist and (dist1 < cub1[index_cub1][6]):
                inter[inter_count][0] = cub1[index_cub1][3]
                inter[inter_count][1] = cub1[index_cub1][4]
                inter[inter_count][2] = cub1[index_cub1][5]
            
                inter_count+=1                            
            
                cub1[index_cub1][3] = x_at
                cub1[index_cub1][4] = y_at
                cub1[index_cub1][5] = z_at
                cub1[index_cub1][6] = dist1
                    
                
            else:
                inter[inter_count][0] = x_at
                inter[inter_count][1] = y_at
                inter[inter_count][2] = z_at
                       
                inter_count+=1 
          
        else:
            
            # We record the atom checked as an interstitial or not
            if dist2 < mini_dist and (abs(cub2[index_cub2][6]) < 1e-8):
                cub2[index_cub2][3] = x_at
                cub2[index_cub2][4] = y_at
                cub2[index_cub2][5] = z_at
                cub2[index_cub2][6] = dist2
                
                lourdeur+=1
                
            elif dist2 < mini_dist and (dist2 < cub2[index_cub2][6]):
                inter[inter_count][0] = cub2[index_cub2][3]
                inter[inter_count][1] = cub2[index_cub2][4]
                inter[inter_count][2] = cub2[index_cub2][5]
                    
                inter_count+=1                            
                
                cub2[index_cub2][3] = x_at
                cub2[index_cub2][4] = y_at
                cub2[index_cub2][5] = z_at
                cub2[index_cub2][6] = dist2
               
            else:
                inter[inter_count][0] = x_at
                inter[inter_count][1] = y_at
                inter[inter_count][2] = z_at
                       
                inter_count+=1 
            
        
    # Now we want to know the vacancies, to do that we just need to identify the points of the theoretical lattice to which no atom has been attributed
    i=0
    while i < (dl.nb_atoms/2):
            
        if abs(cub1[i][3]+cub1[i][4]+cub1[i][5]+cub1[i][6]) < 1e-10:
            vac[vac_count][0] = cub1[i][0]
            vac[vac_count][1] = cub1[i][1]
            vac[vac_count][2] = cub1[i][2]
                
            vac_count+=1
            
            
        if abs(cub2[i][3]+cub2[i][4]+cub2[i][5]+cub2[i][6]) < 1e-10:
            vac[vac_count][0] = cub2[i][0]
            vac[vac_count][1] = cub2[i][1]
            vac[vac_count][2] = cub2[i][2]
            
            vac_count+=1
            
        i+=1
        
    # Adaptation of the size of the numpy arrays that the function will return to not have too large ones
    inter_fitted = np.zeros((inter_count,3))
    vac_fitted = np.zeros((vac_count,3))
    
    inter_fitted = inter[0:inter_count]
    vac_fitted = vac[0:vac_count]
    
    if inter_count != vac_count:
        print('Error: the number of interstitials and vacancies is different')
    
    # Number of Frenkel pairs
    FP_count = vac_count
    
        
    return {'inter':inter_fitted,'vac':vac_fitted,'FP_count':FP_count,'lat_param':lat_param}
    
    









def fcc_point_defects(dl,time):
    '''
    This function is a subpart of the main function called point_defects.
    If the lattice is fcc then this function is called through the main one via
    the use of some if statements.
    Another subpart of the main function (point_defects) called bcc_point_defects
    is used if the lattice is found as beeing bcc. This function is define earlier
    in the script.
    '''
    # Here we define a variable giving the coresponding snapshot for the time provided
    n_snap=time_to_snap(time,dl)

    # Extract the heading indexes of the x position
    index_pos=index_heading(dl,'x')
    #print('index_pos:',index_pos)
    
    # Length of the box
    length=dl.box_dim[0]        
        
    # Determination of the number of cubic cells per edge of the simulation box
    num_edge = round((dl.nb_atoms/4)**(1/3))
    #print('num_edge:',num_edge)
        
    # Determination of the lattice parameter
    lat_param = length / num_edge
    #print('lat_param',lat_param)
    
    # Determination of the minimal distance between 2 atoms
    mini_dist = (math.sqrt(2)/4)*lat_param
    
    # Counter to fill the informations about the vacancies and the interstitials in the defects numpy array
    inter_count = 0
    vac_count = 0
    
    # A fcc lattice is constituted of four interpenatrated cubic lattices
    # From the first cubic lattice we generate the entire fcc lattice by adding three interpenetrating lattices 
    # The three extra lattices are translated from the first one by (1/2a,1/2a,0) or (1/2a,0,1/2a) or (0,1/2a,1/2a)
    # Initialisation of the 1st cubic lattice
    cub1=np.zeros((round(dl.nb_atoms/4),7))
    
    # Initialisation of the 2nd interpenetrating cubic lattice
    cub2=np.zeros((round(dl.nb_atoms/4),7))
    
    # Initialisation of the 3rd interpenetrating cubic lattice
    cub3=np.zeros((round(dl.nb_atoms/4),7))
    
    # Initialisation of the 4th interpenetrating cubic lattice
    cub4=np.zeros((round(dl.nb_atoms/4),7))
    
    
    # Initialisation of the index used to identify the points of the lattice that we intend to generate        
    m=0        
    
    # 3 Loops to generate the points of the lattices by defining the 2 interpenetrating networks cub1 and cub2        
    i=0
    while i < num_edge:
    
        j=0
        while j < num_edge:
        
            k=0
            while k < num_edge:
                cub1[m][0] = i*lat_param
                cub1[m][1] = j*lat_param
                cub1[m][2] = k*lat_param
                
                cub2[m][0] = i*lat_param
                cub2[m][1] = (lat_param/2) + j*lat_param
                cub2[m][2] = (lat_param/2) + k*lat_param
                
                cub3[m][0] = (lat_param/2) + i*lat_param
                cub3[m][1] = j*lat_param
                cub3[m][2] = (lat_param/2) + k*lat_param
                
                cub4[m][0] = (lat_param/2) + i*lat_param
                cub4[m][1] = (lat_param/2) + j*lat_param
                cub4[m][2] = k*lat_param                
                
                m+=1
                k+=1
        
            j+=1
        
        i+=1
    
    
    # Initialisation of the numpy arrays which will contain the coordinates of interstitials and vacancies
    inter = np.zeros((int(dl.nb_atoms),3))
    vac = np.zeros((int(dl.nb_atoms),3))
    
    # Now we need to check all the atom through the dl.data array to attribute them to a lattice point or to define them as interstitials
    i=0
    while i < dl.nb_atoms:
        
        # We define the position variables of the atom to clarify the code
        x_at = dl.data[n_snap][i][index_pos]
        y_at = dl.data[n_snap][i][index_pos+1]
        z_at = dl.data[n_snap][i][index_pos+2]
        
        i+=1
        
        
        # We find the closest lattice point from this atom in cub1
        a1 = int(x_at // lat_param)
        b1 = int(y_at // lat_param)
        c1 = int(z_at // lat_param)
            
            
        if abs(a1*lat_param - x_at) > abs((a1+1)*lat_param - x_at):
            a1+=1
            
        if abs(b1*lat_param - y_at) > abs((b1+1)*lat_param - y_at):
            b1+=1
            
        if abs(c1*lat_param - z_at) > abs((c1+1)*lat_param - z_at):
            c1+=1
            
        a1 = a1%num_edge
        b1 = b1%num_edge
        c1 = c1%num_edge
            
        # Using a, b and c we find the index in the numpy array cub1 of the closest lattice point to the atom
        index_cub1 = int(a1*(num_edge**2) + b1*(num_edge) + c1)
    
    
        # We find the closest lattice point from this atom in cub2
        a2 = int(x_at // lat_param)
        b2 = int((y_at - lat_param/2) // lat_param)
        c2 = int((z_at - lat_param/2) // lat_param)
            
            
        if abs(a2*lat_param - x_at) > abs((a2+1)*lat_param - x_at):
            a2+=1
        
        if abs((b2 + 1/2)*lat_param - y_at) > abs((b2+1 + 1/2)*lat_param - y_at):
            b2+=1
        
        if abs((c2 + 1/2)*lat_param - z_at) > abs((c2+1 + 1/2)*lat_param - z_at):
            c2+=1                
            
        a2 = a2%num_edge
        b2 = b2%num_edge
        c2 = c2%num_edge
            
        # Using a, b and c we find the index in the numpy array cub1 of the closest lattice point to the atom
        index_cub2 = int(a2*(num_edge**2) + b2*(num_edge) + c2)
    
        
        # We find the closest lattice point from this atom in cub3
        a3 = int((x_at - lat_param/2) // lat_param)
        b3 = int(y_at // lat_param)
        c3 = int((z_at - lat_param/2) // lat_param)
            
            
        if abs((a3 + 1/2)*lat_param - x_at) > abs((a3+1 + 1/2)*lat_param - x_at):
            a3+=1
        
        if abs(b3*lat_param - y_at) > abs((b3+1)*lat_param - y_at):
            b3+=1
        
        if abs((c3 + 1/2)*lat_param - z_at) > abs((c3+1 + 1/2)*lat_param - z_at):
            c3+=1                
            
        a3 = a3%num_edge
        b3 = b3%num_edge
        c3 = c3%num_edge
            
        # Using a, b and c we find the index in the numpy array cub1 of the closest lattice point to the atom
        index_cub3 = int(a3*(num_edge**2) + b3*(num_edge) + c3)
        
        
        # We find the closest lattice point from this atom in cub4
        a4 = int((x_at - lat_param/2) // lat_param)
        b4 = int((y_at - lat_param/2) // lat_param)
        c4 = int(z_at // lat_param)
            
            
        if abs((a4 + 1/2)*lat_param - x_at) > abs((a4+1 + 1/2)*lat_param - x_at):
            a4+=1
        
        if abs((b4 + 1/2)*lat_param - y_at) > abs((b4+1 + 1/2)*lat_param - y_at):
            b4+=1
        
        if abs(c4*lat_param - z_at) > abs((c4+1)*lat_param - z_at):
            c4+=1                
            
        a4 = a4%num_edge
        b4 = b4%num_edge
        c4 = c4%num_edge
            
        # Using a, b and c we find the index in the numpy array cub1 of the closest lattice point to the atom
        index_cub4 = int(a4*(num_edge**2) + b4*(num_edge) + c4)
    
        
        # Determination of which of the theoretical lattice points is the closest to the atom checked
        # First we need to calculate the distances between the theoretical point and the atom associated
        # We need to be careful that this distance is not greater than half of the length of the box
        
        # dist1 (use of the intermediate variables x_dist, y_dist and z_dist)
        x_dist=x_at
        y_dist=y_at
        z_dist=z_at
            
        if x_at-cub1[index_cub1][0] < (-length/2):
            x_dist = x_dist + length
        elif x_at-cub1[index_cub1][0] > length/2:
            x_dist = x_dist - length
        
        if y_at-cub1[index_cub1][1] < (-length/2):
            y_dist = y_dist + length
        elif y_at-cub1[index_cub1][1] > length/2:
            y_dist = y_dist - length
        
        if z_at-cub1[index_cub1][2] < (-length/2):
            z_dist = z_dist + length
        elif z_at-cub1[index_cub1][2] > length/2:
            z_dist = z_dist - length
        
        dist1 = math.sqrt( (x_dist-cub1[index_cub1][0])**2 + (y_dist-cub1[index_cub1][1])**2 + (z_dist-cub1[index_cub1][2])**2 )
    
        
        # dist2 (use of the intermediate variables x_dist, y_dist and z_dist)
        x_dist=x_at
        y_dist=y_at
        z_dist=z_at            
        
        if x_at-cub2[index_cub2][0] < (-length/2):
            x_dist = x_dist + length
        elif x_at-cub2[index_cub2][0] > length/2:
            x_dist = x_dist - length
        
        if y_at-cub2[index_cub2][1] < (-length/2):
            y_dist = y_dist + length
        elif y_at-cub2[index_cub2][1] > length/2:
            y_dist = y_dist - length
        
        if z_at-cub2[index_cub2][2] < (-length/2):
            z_dist = z_dist + length
        elif z_at-cub2[index_cub2][2] > length/2:
            z_dist = z_dist - length            
        
        dist2 = math.sqrt( (x_dist-cub2[index_cub2][0])**2 + (y_dist-cub2[index_cub2][1])**2 + (z_dist-cub2[index_cub2][2])**2 )

        
        # dist3 (use of the intermediate variables x_dist, y_dist and z_dist)
        x_dist=x_at
        y_dist=y_at
        z_dist=z_at
            
        if x_at-cub3[index_cub3][0] < (-length/2):
            x_dist = x_dist + length
        elif x_at-cub3[index_cub3][0] > length/2:
            x_dist = x_dist - length
        
        if y_at-cub3[index_cub3][1] < (-length/2):
            y_dist = y_dist + length
        elif y_at-cub3[index_cub3][1] > length/2:
            y_dist = y_dist - length
        
        if z_at-cub3[index_cub3][2] < (-length/2):
            z_dist = z_dist + length
        elif z_at-cub3[index_cub3][2] > length/2:
            z_dist = z_dist - length
        
        dist3 = math.sqrt( (x_dist-cub3[index_cub3][0])**2 + (y_dist-cub3[index_cub3][1])**2 + (z_dist-cub3[index_cub3][2])**2 )

        
        # dist4 (use of the intermediate variables x_dist, y_dist and z_dist)
        x_dist=x_at
        y_dist=y_at
        z_dist=z_at
            
        if x_at-cub4[index_cub4][0] < (-length/2):
            x_dist = x_dist + length
        elif x_at-cub4[index_cub4][0] > length/2:
            x_dist = x_dist - length
        
        if y_at-cub4[index_cub4][1] < (-length/2):
            y_dist = y_dist + length
        elif y_at-cub4[index_cub4][1] > length/2:
            y_dist = y_dist - length
        
        if z_at-cub4[index_cub4][2] < (-length/2):
            z_dist = z_dist + length
        elif z_at-cub4[index_cub4][2] > length/2:
            z_dist = z_dist - length
        
        dist4 = math.sqrt( (x_dist-cub4[index_cub4][0])**2 + (y_dist-cub4[index_cub4][1])**2 + (z_dist-cub4[index_cub4][2])**2 )

        
        # We determine the smallest distance between the 4 defined earlier
        smallest_dist = min(dist1,dist2,dist3,dist4)
        
        # We determine which of the 4 theoretical lattices gives the closest point to the atom
        if abs(dist1 - smallest_dist) < 1e-10:
                
            # We record the atom checked as an interstitial or not
            if (dist1 < mini_dist) and (abs(cub1[index_cub1][6]) < 1e-8):
                cub1[index_cub1][3] = x_at
                cub1[index_cub1][4] = y_at
                cub1[index_cub1][5] = z_at
                cub1[index_cub1][6] = dist1
                
            elif dist1 < mini_dist and (dist1 < cub1[index_cub1][6]):
                inter[inter_count][0] = cub1[index_cub1][3]
                inter[inter_count][1] = cub1[index_cub1][4]
                inter[inter_count][2] = cub1[index_cub1][5]
            
                inter_count+=1                            
            
                cub1[index_cub1][3] = x_at
                cub1[index_cub1][4] = y_at
                cub1[index_cub1][5] = z_at
                cub1[index_cub1][6] = dist1
                    
                
            else:
                inter[inter_count][0] = x_at
                inter[inter_count][1] = y_at
                inter[inter_count][2] = z_at
                       
                inter_count+=1 
            
        elif abs(dist2 - smallest_dist) < 1e-10:
            
            # We record the atom checked as an interstitial or not
            if dist2 < mini_dist and (abs(cub2[index_cub2][6]) < 1e-8):
                cub2[index_cub2][3] = x_at
                cub2[index_cub2][4] = y_at
                cub2[index_cub2][5] = z_at
                cub2[index_cub2][6] = dist2
                
            elif dist2 < mini_dist and (dist2 < cub2[index_cub2][6]):
                inter[inter_count][0] = cub2[index_cub2][3]
                inter[inter_count][1] = cub2[index_cub2][4]
                inter[inter_count][2] = cub2[index_cub2][5]
                    
                inter_count+=1                            
                
                cub2[index_cub2][3] = x_at
                cub2[index_cub2][4] = y_at
                cub2[index_cub2][5] = z_at
                cub2[index_cub2][6] = dist2
               
            else:
                inter[inter_count][0] = x_at
                inter[inter_count][1] = y_at
                inter[inter_count][2] = z_at
                       
                inter_count+=1 
        
        elif abs(dist3 - smallest_dist) < 1e-10:
            
            # We record the atom checked as an interstitial or not
            if dist3 < mini_dist and (abs(cub3[index_cub3][6]) < 1e-8):
                cub3[index_cub3][3] = x_at
                cub3[index_cub3][4] = y_at
                cub3[index_cub3][5] = z_at
                cub3[index_cub3][6] = dist3
                
            elif dist3 < mini_dist and (dist3 < cub3[index_cub3][6]):
                inter[inter_count][0] = cub3[index_cub3][3]
                inter[inter_count][1] = cub3[index_cub3][4]
                inter[inter_count][2] = cub3[index_cub3][5]
                    
                inter_count+=1                            
                
                cub3[index_cub3][3] = x_at
                cub3[index_cub3][4] = y_at
                cub3[index_cub3][5] = z_at
                cub3[index_cub3][6] = dist3
               
            else:
                inter[inter_count][0] = x_at
                inter[inter_count][1] = y_at
                inter[inter_count][2] = z_at
                       
                inter_count+=1 
        
        else:
            
            # We record the atom checked as an interstitial or not
            if dist4 < mini_dist and (abs(cub4[index_cub4][6]) < 1e-8):
                cub4[index_cub4][3] = x_at
                cub4[index_cub4][4] = y_at
                cub4[index_cub4][5] = z_at
                cub4[index_cub4][6] = dist4
                
            elif dist4 < mini_dist and (dist4 < cub4[index_cub4][6]):
                inter[inter_count][0] = cub4[index_cub4][3]
                inter[inter_count][1] = cub4[index_cub4][4]
                inter[inter_count][2] = cub4[index_cub4][5]
                    
                inter_count+=1                            
                
                cub4[index_cub4][3] = x_at
                cub4[index_cub4][4] = y_at
                cub4[index_cub4][5] = z_at
                cub4[index_cub4][6] = dist4
               
            else:
                inter[inter_count][0] = x_at
                inter[inter_count][1] = y_at
                inter[inter_count][2] = z_at
                       
                inter_count+=1
        
    # Now we want to know the vacancies, to do that we just need to identify the points of the theoretical lattice to which no atom has been attributed
    i=0
    while i < int(dl.nb_atoms/4):
            
        if abs(cub1[i][3]+cub1[i][4]+cub1[i][5]+cub1[i][6]) < 1e-10:
            vac[vac_count][0] = cub1[i][0]
            vac[vac_count][1] = cub1[i][1]
            vac[vac_count][2] = cub1[i][2]
                
            vac_count+=1
               
               
        if abs(cub2[i][3]+cub2[i][4]+cub2[i][5]+cub2[i][6]) < 1e-10:
            vac[vac_count][0] = cub2[i][0]
            vac[vac_count][1] = cub2[i][1]
            vac[vac_count][2] = cub2[i][2]
            
            vac_count+=1
        
        
        if abs(cub3[i][3]+cub3[i][4]+cub3[i][5]+cub3[i][6]) < 1e-10:
            vac[vac_count][0] = cub3[i][0]
            vac[vac_count][1] = cub3[i][1]
            vac[vac_count][2] = cub3[i][2]
            
            vac_count+=1
        
        
        if abs(cub4[i][3]+cub4[i][4]+cub4[i][5]+cub4[i][6]) < 1e-10:
            vac[vac_count][0] = cub4[i][0]
            vac[vac_count][1] = cub4[i][1]
            vac[vac_count][2] = cub4[i][2]
            
            vac_count+=1
              
        i+=1
               
        
    # Adaptation of the size of the numpy arrays that the function will return to not have too large ones
    inter_fitted = np.zeros((inter_count,3))
    vac_fitted = np.zeros((vac_count,3))
    
    inter_fitted = inter[0:inter_count]
    vac_fitted = vac[0:vac_count]
    
    
    if inter_count != vac_count:
        print('Error: the number of interstitials and vacancies is different')
    
    # Number of Frenkel pairs
    FP_count = vac_count
    
        
    return {'inter':inter_fitted,'vac':vac_fitted,'FP_count':FP_count,'lat_param':lat_param}
    
    
    
    
    
    
    
    
    
    
    

def bcc_inter_clusters(self,NN):
    '''
    This is a subpart of the clusters() method of the point_defects class
    '''
    
    # recovery of the number of vacancies
    inter_count=self.defects['FP_count']    
    
    # recovery of the lattice parameter
    lat_param = self.defects['lat_param'] 
    
    # initialisation of a list to store the indexes of the clusters and the atoms of which they are constituted
    stor_index=np.zeros((inter_count,inter_count))
    
    # Initialisation of the shape list which contain one element per cluster and each element indicates the number of atoms included in the related cluster
    shape_list=[]
    
    # Distance cut-off between 2 atoms below which they are considered as neighbours
    if NN == 1:
        cut_off = ((2+math.sqrt(3))/4) * lat_param
    elif NN == 2:
        cut_off = ((1+math.sqrt(2))/2) * lat_param
    elif NN == 3:
        cut_off = ((math.sqrt(8)+math.sqrt(11))/(2*math.sqrt(4))) * lat_param
    elif NN == 4:
        cut_off = ((math.sqrt(12)+math.sqrt(11))/(2*math.sqrt(4))) * lat_param
    else:
        print('Nearest neighbour variable NN is invalid in bcc_inter_clusters')
        
    
    # initialisation of two lists which will be useful in the following double loop
    forbidden=[]
    index_list=[]
    
    test=1
    memory=0
    nb_clusters=0
    highest_length=0
    i=0
    while i < inter_count:
        
        x = self.defects['inter'][i,0]
        y = self.defects['inter'][i,1]
        z = self.defects['inter'][i,2]
        
        if i not in forbidden:
            forbidden+=[i]
            index_list+=[i]        
        
        if test == 1:
            memory=i
        
        j=i+1
        while j in forbidden:
            j+=1
        
        while j < inter_count:
            
            dist = math.sqrt((x - self.defects['inter'][j,0])**2 + (y - self.defects['inter'][j,1])**2 + (z - self.defects['inter'][j,2])**2)
            
            if dist < cut_off:
                
                forbidden+=[j]
                index_list+=[j]
            
            j+=1
            while j in forbidden:
                j+=1
        
        if test < len(index_list):
            i=index_list[test]
            test+=1
        else:
            test=1
            
            # Iteration of the number of clusters
            if len(index_list) > 1:
                shape_list+=[len(index_list)]
                stor_index[nb_clusters][0:len(index_list)]=index_list
                nb_clusters+=1
            
            # Determination of the number of atoms in the biggest cluster to fit the shape of the following numpy array called clusters
            if len(index_list) > highest_length:
                highest_length=len(index_list)
            index_list=[]
            
            i=memory+1
            while i in forbidden:
                i+=1
        
    
    
    clusters=np.zeros((nb_clusters,highest_length,3))
    
    i=0
    while i < nb_clusters:
        inclustered_atoms = shape_list[i]
        j=0
        while j < inclustered_atoms:
            # Index of the numpy array that contains the coordinates of the vacancies, steming from the point_defects class
            index = round(stor_index[i,j])
            clusters[i,j,0:3] = self.defects['vac'][index,0:3]
            
            j+=1
        
        i+=1
    
    #print('forbidden''s length',len(forbidden),'forbidden',forbidden)
    
    return {'clusters':clusters,'shape':shape_list}
    
    
    
    
    
    
    
    
    
    
    

def bcc_vac_clusters(self,NN):
    '''
    This is a subpart of the clusters() method of the point_defects class
    '''
    
    # recovery of the number of vacancies
    vac_count=self.defects['FP_count']    
    
    # recovery of the lattice parameter
    lat_param = self.defects['lat_param'] 
    
    # initialisation of a list to store the indexes of the clusters and the atoms of which they are constituted
    stor_index=np.zeros((vac_count,vac_count))
    
    # Initialisation of the shape list which contain one element per cluster and each element indicates the number of atoms included in the related cluster
    shape_list=[]
    
    # Distance cut-off between 2 atoms below which they are considered as neighbours
    if NN == 1:
        cut_off = ((2+math.sqrt(3))/4) * lat_param
    elif NN == 2:
        cut_off = ((1+math.sqrt(2))/2) * lat_param
    elif NN == 3:
        cut_off = ((math.sqrt(8)+math.sqrt(11))/(2*math.sqrt(4))) * lat_param
    elif NN == 4:
        cut_off = ((math.sqrt(12)+math.sqrt(11))/(2*math.sqrt(4))) * lat_param
    else:
        print('Nearest neighbour variable NN is invalid in bcc_inter_clusters')
        
    
    # initialisation of two lists which will be useful in the following double loop
    forbidden=[]
    index_list=[]
    
    test=1
    memory=0
    nb_clusters=0
    highest_length=0
    i=0
    while i < vac_count:
        
        x = self.defects['vac'][i,0]
        y = self.defects['vac'][i,1]
        z = self.defects['vac'][i,2]
        
        if i not in forbidden:
            forbidden+=[i]
            index_list+=[i]        
        
        if test == 1:
            memory=i
        
        j=i+1
        while j in forbidden:
            j+=1
        
        while j < vac_count:
            
            dist = math.sqrt((x - self.defects['vac'][j,0])**2 + (y - self.defects['vac'][j,1])**2 + (z - self.defects['vac'][j,2])**2)
            
            if dist < cut_off:
                
                forbidden+=[j]
                index_list+=[j]
            
            j+=1
            while j in forbidden:
                j+=1
        
        if test < len(index_list):
            i=index_list[test]
            test+=1
        else:
            test=1
            
            # Iteration of the number of clusters
            if len(index_list) > 1:
                shape_list+=[len(index_list)]
                stor_index[nb_clusters][0:len(index_list)]=index_list
                nb_clusters+=1
            
            # Determination of the number of atoms in the biggest cluster to fit the shape of the following numpy array called clusters
            if len(index_list) > highest_length:
                highest_length=len(index_list)
            index_list=[]
            
            i=memory+1
            while i in forbidden:
                i+=1
        
    
    
    clusters=np.zeros((nb_clusters,highest_length,3))
    
    i=0
    while i < nb_clusters:
        inclustered_atoms = shape_list[i]
        j=0
        while j < inclustered_atoms:
            # Index of the numpy array that contains the coordinates of the vacancies, steming from the point_defects class
            index = round(stor_index[i,j])
            clusters[i,j,0:3] = self.defects['vac'][index,0:3]
            
            j+=1
        
        i+=1
    
    #print('forbidden''s length',len(forbidden),'forbidden',forbidden)    
    
    
    return {'clusters':clusters,'shape':shape_list}
    









def fcc_inter_clusters(self,NN):
    '''
    This is a subpart of the clusters() method of the point_defects class
    '''
    
    # recovery of the number of vacancies
    inter_count=self.defects['FP_count']    
    
    # recovery of the lattice parameter
    lat_param = self.defects['lat_param'] 
    
    # initialisation of a list to store the indexes of the clusters and the atoms of which they are constituted
    stor_index=np.zeros((inter_count,inter_count))
    
    # Initialisation of the shape list which contain one element per cluster and each element indicates the number of atoms included in the related cluster
    shape_list=[]
    
    # Distance cut-off between 2 atoms below which they are considered as neighbours
    if NN == 1:
        cut_off = ((math.sqrt(2)+1)/(2*math.sqrt(2))) * lat_param
    elif NN == 2:
        cut_off = ((math.sqrt(3)+math.sqrt(2))/(2*math.sqrt(2))) * lat_param
    elif NN == 3:
        cut_off = ((math.sqrt(3)+2)/(2*math.sqrt(2))) * lat_param
    elif NN == 4:
        cut_off = ((math.sqrt(5)+2)/(2*math.sqrt(2))) * lat_param
    else:
        print('Nearest neighbour variable NN is invalid in bcc_inter_clusters')
        
    
    # initialisation of two lists which will be useful in the following double loop
    forbidden=[]
    index_list=[]
    
    test=1
    memory=0
    nb_clusters=0
    highest_length=0
    i=0
    while i < inter_count:
        
        x = self.defects['inter'][i,0]
        y = self.defects['inter'][i,1]
        z = self.defects['inter'][i,2]
        
        if i not in forbidden:
            forbidden+=[i]
            index_list+=[i]        
        
        if test == 1:
            memory=i
        
        j=i+1
        while j in forbidden:
            j+=1
        
        while j < inter_count:
            
            dist = math.sqrt((x - self.defects['inter'][j,0])**2 + (y - self.defects['inter'][j,1])**2 + (z - self.defects['inter'][j,2])**2)
            
            if dist < cut_off:
                
                forbidden+=[j]
                index_list+=[j]
            
            j+=1
            while j in forbidden:
                j+=1
        
        if test < len(index_list):
            i=index_list[test]
            test+=1
        else:
            test=1
            
            # Iteration of the number of clusters
            if len(index_list) > 1:
                shape_list+=[len(index_list)]
                stor_index[nb_clusters][0:len(index_list)]=index_list
                nb_clusters+=1
            
            # Determination of the number of atoms in the biggest cluster to fit the shape of the following numpy array called clusters
            if len(index_list) > highest_length:
                highest_length=len(index_list)
            index_list=[]
            
            i=memory+1
            while i in forbidden:
                i+=1
        
    
    
    clusters=np.zeros((nb_clusters,highest_length,3))
    
    i=0
    while i < nb_clusters:
        inclustered_atoms = shape_list[i]
        j=0
        while j < inclustered_atoms:
            # Index of the numpy array that contains the coordinates of the vacancies, steming from the point_defects class
            index = round(stor_index[i,j])
            clusters[i,j,0:3] = self.defects['vac'][index,0:3]
            
            j+=1
        
        i+=1
        
    
    return {'clusters':clusters,'shape':shape_list}







    
    
def fcc_vac_clusters(self,NN):
    '''
    This is a subpart of the clusters() method of the point_defects class
    '''
    
    # recovery of the number of vacancies
    vac_count=self.defects['FP_count']    
    
    # recovery of the lattice parameter
    lat_param = self.defects['lat_param'] 
    
    # initialisation of a list to store the indexes of the clusters and the atoms of which they are constituted
    stor_index=np.zeros((vac_count,vac_count))
    
    # Initialisation of the shape list which contain one element per cluster and each element indicates the number of atoms included in the related cluster
    shape_list=[]
    
    # Distance cut-off between 2 atoms below which they are considered as neighbours
    if NN == 1:
        cut_off = ((math.sqrt(2)+1)/(2*math.sqrt(2))) * lat_param
    elif NN == 2:
        cut_off = ((math.sqrt(3)+math.sqrt(2))/(2*math.sqrt(2))) * lat_param
    elif NN == 3:
        cut_off = ((math.sqrt(3)+2)/(2*math.sqrt(2))) * lat_param
    elif NN == 4:
        cut_off = ((math.sqrt(5)+2)/(2*math.sqrt(2))) * lat_param
    else:
        print('Nearest neighbour variable NN is invalid in bcc_inter_clusters')
        
    
    # initialisation of two lists which will be useful in the following double loop
    forbidden=[]
    index_list=[]
    
    test=1
    memory=0
    nb_clusters=0
    highest_length=0
    i=0
    while i < vac_count:
        
        x = self.defects['vac'][i,0]
        y = self.defects['vac'][i,1]
        z = self.defects['vac'][i,2]
        
        if i not in forbidden:
            forbidden+=[i]
            index_list+=[i]        
        
        if test == 1:
            memory=i
        
        j=i+1
        while j in forbidden:
            j+=1
        
        while j < vac_count:
            
            dist = math.sqrt((x - self.defects['vac'][j,0])**2 + (y - self.defects['vac'][j,1])**2 + (z - self.defects['vac'][j,2])**2)
            
            # These following conditions in the if statement need to be verified
            if dist < cut_off:
                
                forbidden+=[j]
                index_list+=[j]
            
            j+=1
            while j in forbidden:
                j+=1
        
        if test < len(index_list):
            i=index_list[test]
            test+=1
        else:
            test=1
            
            # Iteration of the number of clusters
            if len(index_list) > 1:
                shape_list+=[len(index_list)]
                stor_index[nb_clusters][0:len(index_list)]=index_list
                nb_clusters+=1
            
            # Determination of the number of atoms in the biggest cluster to fit the shape of the following numpy array called clusters
            if len(index_list) > highest_length:
                highest_length=len(index_list)
            index_list=[]
            
            i=memory+1
            while i in forbidden:
                i+=1
        
    
    
    clusters=np.zeros((nb_clusters,highest_length,3))
    
    i=0
    while i < nb_clusters:
        inclustered_atoms = shape_list[i]
        j=0
        while j < inclustered_atoms:
            # Index of the numpy array that contains the coordinates of the vacancies, steming from the point_defects class
            index = round(stor_index[i,j])
            clusters[i,j,0:3] = self.defects['vac'][index,0:3]
            
            j+=1
        
        i+=1
        
    
    return {'clusters':clusters,'shape':shape_list}








